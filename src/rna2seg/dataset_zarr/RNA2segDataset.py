#%%

import json
import logging
from pathlib import Path
from time import time
from typing import Optional

import albumentations as A
import cv2
import numpy as np
import pandas as pd
import spatialdata as sd
import tifffile
import torch
from albumentations.core.transforms_interface import ImageOnlyTransform
from cellpose import transforms as tf_cp
from matplotlib import colors as mcolors
from scipy import ndimage as ndi
from sopa._constants import SopaKeys
from sopa.utils.utils import to_intrinsic
from torch.utils.data import Dataset
from tqdm import tqdm

from rna2seg._constant import RNA2segFiles
# from rna2seg.dataset_zarr import StainingTranscriptSegmentation
from rna2seg.dataset_zarr.background import get_background_mask
from rna2seg.dataset_zarr.data_augmentation import (cellbound_transform,
                                                    random_rotate_and_resize)
from rna2seg.utils import create_cell_contours

log = logging.getLogger(__name__)

import dask

dask.config.set(scheduler='synchronous')


def create_augmented_patch_df(patch_df, factor=2, std_dev=5):
    new_rows = []
    for idx, row in patch_df.iterrows():
        for i in range(factor - 1):
            x_new = row['x'] + np.random.normal(0, std_dev)
            y_new = row['y'] + np.random.normal(0, std_dev)
            new_row = row.copy()
            new_row['x'] = x_new
            new_row['y'] = y_new
            new_rows.append(new_row)

    new_df = pd.DataFrame(new_rows)
    return pd.concat([patch_df, new_df], ignore_index=True)


def pad_image3D_to_shape(image, target_shape=(1200, 1200)):
    """
    Pad the input image with zeros to match the target shape.

    Parameters:
    image (np.ndarray): Input 3D image.
    target_shape (tuple): Desired shape of the output image (height, width, channels).

    Returns:
    np.ndarray: Padded image.
    """
    assert len(image.shape) == 3, "Input image should be 3D"
    current_shape = image.shape
    padding = (
        (0, target_shape[0] - current_shape[1]),  # Padding for rows
        (0, target_shape[1] - current_shape[2])  # Padding for columns
    )
    new_image = np.zeros((image.shape[0], target_shape[0], target_shape[1]), dtype=image.dtype)
    for i in range(current_shape[0]):
        new_image[i] = pad_image2D_to_shape(image[i], target_shape=target_shape)
    return new_image


def pad_image2D_to_shape(image, target_shape=(1200, 1200)):
    """
    Pad the input image with zeros to match the target shape.

    Parameters:
    image (np.ndarray): Input 2D image.
    target_shape (tuple): Desired shape of the output image (height, width).

    Returns:
    np.ndarray: Padded image.
    """
    current_shape = image.shape
    padding = (
        (0, target_shape[0] - current_shape[0]),  # Padding for rows
        (0, target_shape[1] - current_shape[1])  # Padding for columns
    )
    padded_image = np.pad(image, padding, mode='constant', constant_values=0)
    return padded_image


class MaxFilter(ImageOnlyTransform):
    def __init__(self, max_filter_size=0, always_apply=False, p=1.0):
        super(MaxFilter, self).__init__(always_apply, p)
        self.max_filter_size = max_filter_size

    def apply(self, img, **params):
        if self.max_filter_size > 0:
            img = ndi.maximum_filter(img, size=(self.max_filter_size, self.max_filter_size, 1))
        return img


def compute_array_coord(bounds,
                        patch_df,
                        image_shape,
                        scaling_factor_coord):
    offset_x = bounds[0]
    offset_y = bounds[1]
    column_y = "y"
    column_x = "x"
    rna_df = patch_df.copy()
    rna_df[column_y] = rna_df[column_y].astype(int) - offset_y
    rna_df[column_x] = rna_df[column_x].astype(int) - offset_x
    ## drop the points outside the image
    str_assert = "rna coordinates max does not match image shape, check preprocessing :"
    assert (rna_df[column_y] >= 0).all(), str_assert
    assert (rna_df[column_y] < image_shape[0]).all(), str_assert
    assert (rna_df[column_x] >= 0).all, str_assert
    assert (rna_df[column_x] < image_shape[1]).all(), str_assert
    array_coord = np.array(rna_df[[column_y, column_x]].values)
    assert max(array_coord[:, 0]) < bounds[
        3], f"rna coordiante max doe not match image shape, check preprocessing : {max(array_coord[:, 0])} >= {max(bounds[3])}"
    assert max(array_coord[:, 1]) < bounds[
        2], f"rna coordiante max doe not match image shape, check preprocessing : {max(array_coord[:, 1])} >= {max(bounds[2])}"
    assert np.min(array_coord) >= 0, f" negative coordinate not allowed in patches"
    ## resize coordinate
    array_coord = (array_coord * scaling_factor_coord).astype(int)
    return array_coord


def rna2img(df_crop, dict_gene_value: dict | None,
            image_shape=(1200, 1200, 3),
            gene_column="genes",
            column_x="x",
            column_y="y",
            offset_x=0,
            offset_y=0,
            gaussian_kernel_size: float = 0,
            max_filter_size: float = 0,
            addition_mode=False):
    """
    Convert a dataframe of rna coordinates to an image
    Args:
        df_crop: dataframe with the rna coordinates
        image_shape: shape of the output image
        dict_gene_value: dictionary to convert gene to value
        gene_column: column name of the gene
        column_x: column name of the x coordinate
        column_y: column name of the y coordinate
    Returns:
    """

    df_crop = df_crop.copy()
    df_crop[column_y] = df_crop[column_y].astype(int) - offset_y
    df_crop[column_x] = df_crop[column_x].astype(int) - offset_x
    ## drop the points outside the image
    df_crop = df_crop[(df_crop[column_y] >= 0) & (df_crop[column_y] < image_shape[0])]
    df_crop = df_crop[(df_crop[column_x] >= 0) & (df_crop[column_x] < image_shape[1])]

    list_y = df_crop[column_y].values
    list_y = list_y.astype(int)
    list_x = df_crop[column_x].values
    list_x = list_x.astype(int)
    list_gene = df_crop[gene_column].values
    if dict_gene_value is not None:
        list_coord_value = [dict_gene_value[gene] for gene in list_gene]
    else:
        list_coord_value = [1 for gene in list_gene]

    img = np.zeros(image_shape, dtype=np.float32)
    if len(list_y) == 0:
        assert len(list_x) == 0
        assert len(list_coord_value) == 0
    else:
        if not addition_mode:
            # raise NotImplementedError("not using addition mode is deprecated")
            img[list_y, list_x] = list_coord_value
        else:
            for i, (y, x) in enumerate(zip(list_y, list_x)):
                img[y, x] += list_coord_value[i]

    # add max filter with scipy
    if max_filter_size > 0:
        img = ndi.maximum_filter(img, size=(max_filter_size, max_filter_size, 1))

    if gaussian_kernel_size > 0:  # it would be more optimise to apply it after resizing/downsampling
        img = ndi.gaussian_filter(img, sigma=(gaussian_kernel_size, gaussian_kernel_size, 0))
    return img


def remove_cell_in_background(agreement_segmentation,
                              background,
                              threshold_cell_in_bg=0.05,
                              threshold_nuclei_in_bg=0.5,
                              agreement_segmentation_without_nuclei=None,
                              segmentation_nuclei_not_in_cell=None,
                              remove_nucleus_seg_from_bg=True,
                              segmentation_nuclei=None,
                              ):
    """
    Remove the cells that are in the background
    Args:
        agreement_segmentation: mask of cell segmentation
        background: mask of background
        threshold_cell_in_bg:  maximum % of the cell in background to be accepted and taken into account during the backpropagation
    Returns:
    """
    t = time()

    def _remove_cell_in_background(agreement_segmentation, background, threshold_bg=0.5):
        overlap_bg_cell = ((agreement_segmentation > 0).astype(int) + (background > 0).astype(int)) == 2
        overlap_bg_cell = overlap_bg_cell * agreement_segmentation
        overlap_bg_cell_unique, overlap_bg_cell_count = np.unique(overlap_bg_cell, return_counts=True)
        if 0 in overlap_bg_cell_unique:
            assert overlap_bg_cell_unique[0] == 0
            overlap_bg_cell_unique = overlap_bg_cell_unique[1:]
            overlap_bg_cell_count = overlap_bg_cell_count[1:]
        ag_unique, ag_count = np.unique(agreement_segmentation, return_counts=True)
        if 0 in ag_unique:
            assert ag_unique[0] == 0
            ag_unique = list(ag_unique[1:])
            ag_count = ag_count[1:]
        new_ag_count = []
        for i in overlap_bg_cell_unique:
            new_ag_count.append(ag_count[ag_unique.index(i)])

        overlap_bg_cell = overlap_bg_cell_count / np.array(new_ag_count)

        cell_to_remove = overlap_bg_cell > threshold_bg
        overlap_bg_cell_unique = overlap_bg_cell_unique[cell_to_remove]

        corrected_agreement_segmentation = agreement_segmentation.copy()
        corrected_agreement_segmentation[np.isin(agreement_segmentation, overlap_bg_cell_unique)] = 0
        return corrected_agreement_segmentation

    if agreement_segmentation_without_nuclei is not None:
        corrected_agreement_segmentation_without_nuclei = _remove_cell_in_background(
            agreement_segmentation=agreement_segmentation_without_nuclei,
            background=background,
            threshold_bg=threshold_cell_in_bg)
    else:
        corrected_agreement_segmentation_without_nuclei = np.zeros_like(agreement_segmentation).astype(int)

    if segmentation_nuclei_not_in_cell is not None:
        corrected_segmentation_nuclei_not_in_cell = _remove_cell_in_background(
            agreement_segmentation=segmentation_nuclei_not_in_cell,
            background=background,
            threshold_bg=threshold_nuclei_in_bg)
        # remove the overlap between target
        corrected_segmentation_nuclei_not_in_cell[corrected_agreement_segmentation_without_nuclei > 0] = 0  #
        corrected_segmentation_nuclei_not_in_cell[agreement_segmentation > 0] = 0
    else:
        corrected_segmentation_nuclei_not_in_cell = np.zeros_like(agreement_segmentation).astype(int)

    agreement_segmentation = agreement_segmentation.astype(int)
    corrected_agreement_segmentation_without_nuclei = corrected_agreement_segmentation_without_nuclei.astype(int)
    corrected_segmentation_nuclei_not_in_cell = corrected_segmentation_nuclei_not_in_cell.astype(int)

    if remove_nucleus_seg_from_bg:
        assert segmentation_nuclei is not None, "segmentation_nuclei should be provided if remove_nucleus_seg_from_bg is True"
        correct_background = (background > 0).astype(int) - (
                corrected_agreement_segmentation_without_nuclei > 0).astype(int) \
                             - (agreement_segmentation > 0).astype(int) \
                             - (corrected_segmentation_nuclei_not_in_cell > 0).astype(int) - \
                             (segmentation_nuclei > 0).astype(int)
    else:
        correct_background = (background > 0).astype(int) - (
                corrected_agreement_segmentation_without_nuclei > 0).astype(int) \
                             - (agreement_segmentation > 0).astype(int) - (
                                     corrected_segmentation_nuclei_not_in_cell > 0).astype(int)
    correct_background = (correct_background > 0).astype(int)

    max_indice_cell = np.max(agreement_segmentation).astype(int)
    corrected_agreement_segmentation_without_nuclei[
        corrected_agreement_segmentation_without_nuclei > 0] += max_indice_cell
    max_indice_cell = np.max(corrected_agreement_segmentation_without_nuclei)
    corrected_segmentation_nuclei_not_in_cell[corrected_segmentation_nuclei_not_in_cell > 0] += max_indice_cell
    # max_indice_cell = np.max(segmentation_nuclei_not_in_cell)

    corrected_agreement_segmentation = agreement_segmentation + corrected_agreement_segmentation_without_nuclei + corrected_segmentation_nuclei_not_in_cell
    # mask_gradient = (corrected_agreement_segmentation > 0).astype(int) + (correct_background>0).astype(int)
    mask_gradient = ((corrected_agreement_segmentation > 0).astype(int) +
                     (corrected_agreement_segmentation_without_nuclei > 0).astype(int) +
                     (corrected_segmentation_nuclei_not_in_cell > 0).astype(int) * 3 +
                     (correct_background > 0).astype(int))

    # print(f"Time to remove cell in background: {time() - t:.6f}s")

    return mask_gradient, corrected_agreement_segmentation, correct_background  #, corrected_agreement_segmentation_without_nuclei


class RNA2segDataset(Dataset):

    def __init__(self,

                 sdata: sd.SpatialData,
                 channels_dapi: list[str],  # should be just STR
                 channels_cellbound: list[str] | None = None,

                 key_cell_consistent: str | None = None,
                 key_nucleus_consistent: str | None = None,
                 key_nuclei_segmentation: str | None = None,  # use to compute background

                 dict_gene_value: dict | None = None,

                 training_mode: bool = False,
                 evaluation_mode: bool = False,

                 patch_width: int = None,
                 patch_overlap: int = None,

                 list_patch_index: list[int] | None = None,
                 list_annotation_patches: list[int] | None = None,
                 gene_column="gene",
                 density_threshold: float | None = None,

                 kernel_size_background_density: None | float = 5,
                 kernel_size_rna2img: float = 0.5,
                 max_filter_size_rna2img: float = 2,
                 transform_resize: Optional[A.Compose] = None,
                 transform_dapi: Optional[A.Compose] = None,
                 augment_rna_density: bool = False,

                 min_nb_cell_per_patch=1,
                 remove_cell_in_background_threshold=0.05,
                 remove_nucleus_seg_from_bg=True,
                 addition_mode=True,
                 min_transcripts: int = 1,

                 # for rna_encoding
                 return_df=False,
                 gene2index=None,
                 augmentation_img=False,

                 # save flow
                 recompute_flow=True,

                 # for testing
                 test_return_background=False,

                 # for save cache

                 # path_cache
                 patch_dir: Path | str | None = None,  # Transcripts
                 experiment_name='input_target_rna2seg',  # for dev only
                 use_cache=False,
                 # IMG AUGMENTATION

                 # optional for testing
                 shape_patch_key=None,

                 ):

        """
        Initializes the dataset with the provided parameters.

        :param sdata: SpatialData object containing spatial transcriptomics data.
        :param channels_dapi: List of DAPI channels for nuclear staining.
        :param channels_cellbound: List of cell boundary channels, or None if not provided.
        :param key_cell_consistent: Key for consistent cell segmentation, or None.
        :param key_nucleus_consistent: Key for consistent nucleus segmentation, or None.
        :param key_nuclei_segmentation: Key for nuclei segmentation, or None.
        :param dict_gene_value: Dictionary containing gene encodings, or None.
        :param training_mode: Boolean flag to enable training mode.
        :param evaluation_mode: Boolean flag to enable evaluation mode.
        :param patch_width: Width of the patch for segmentation, or None.
        :param patch_overlap: Overlap of the patch for segmentation, or None.
        :param list_patch_index: List of patch indices, or None.
        :param list_annotation_patches: List of annotation patches to exclude, or None.
        :param gene_column: Column name for genes, default is "gene".
        :param density_threshold: Threshold for density calculation, or None.
        :param kernel_size_background_density: Kernel size for background density calculation, default is 5.
        :param kernel_size_rna2img: Gaussian kernel size for RNA image transformation, default is 0.5.
        :param max_filter_size_rna2img: Max filter size for RNA image transformation, default is 2.
        :param transform_resize: Resize transformation function for images, or None.
        :param transform_dapi: Transformation function for DAPI images, or None.
        :param augment_rna_density: Boolean flag to enable RNA density augmentation.
        :param min_nb_cell_per_patch: Minimum number of cells per patch for inclusion.
        :param remove_cell_in_background_threshold: Threshold for removing cells in the background.
        :param remove_nucleus_seg_from_bg: Boolean flag to remove nucleus segmentation from the background.
        :param addition_mode: Boolean flag to enable RNA spot value addition.
        :param return_df: Boolean flag to return DataFrame, default is False.
        :param gene2index: Dictionary mapping genes to indices, or None if `return_df` is False.
        :param augmentation_img: Boolean flag to enable image augmentation.
        :param recompute_flow: Boolean flag to recompute the flow field.
        :param test_return_background: Boolean flag to return background image for testing.
        :param patch_dir: Directory for patch storage, or None.
        :param experiment_name: Name of the experiment, default is 'input_target_rna2seg'.
        :param use_cache: Boolean flag to enable cache usage.
        :param shape_patch_key: Key for patch shape, or None.
        """

        self.training_mode = training_mode  # should be imporve
        self.evaluation_mode = evaluation_mode

        if self.training_mode:
            self.return_agg_segmentation = True
            self.return_flow = True

        else:
            self.return_agg_segmentation = False
            self.return_flow = False

        if self.evaluation_mode:
            self.return_agg_segmentation = True
            self.return_flow = True

        assert channels_cellbound is not None, "channels_cellbound should be provided"
        # path parameter
        if shape_patch_key is None:
            shape_patch_key = f"sopa_patches_rna2seg_{patch_width}_{patch_overlap}"
            print(f'default shape_patch_key set to {shape_patch_key}')
        assert shape_patch_key in sdata, f"shape_patch_key {shape_patch_key} not in sdata, set the correct shape_patch_key"

        # for compatibility with the sopa 2
        sdata[SopaKeys.PATCHES] = sdata[shape_patch_key]

        self.key_nuclei_segmentation = key_nuclei_segmentation

        from rna2seg.dataset_zarr import StainingTranscriptSegmentation

        st_segmentation = StainingTranscriptSegmentation(
            sdata=sdata,
            channels_dapi=channels_dapi,
            channels_cellbound=channels_cellbound,
            key_nuclei_segmentation=key_nuclei_segmentation,  # use to compute background
            key_cell_consistent=key_cell_consistent,
            key_nucleus_consistent=key_nucleus_consistent,
            density_threshold=density_threshold,
            patch_dir=patch_dir,
            shape_patch_key=shape_patch_key,
        )
        self.st_segmentation = st_segmentation
        self.st_segmentation.density_threshold = density_threshold # dummy to remove

        self.min_nb_cell_per_patch = min_nb_cell_per_patch
        self.list_patch_index = list_patch_index
        self.min_transcripts = min_transcripts

        self.gene_column = gene_column
        self.dict_gene_value = dict_gene_value
        self.kernel_size_background_density = kernel_size_background_density  # should be the same as the one used to compute density_threshold

        self.kernel_size_rna2img = kernel_size_rna2img
        self.max_filter_size_rna2img = max_filter_size_rna2img

        self.transform_resize = transform_resize
        if self.transform_resize is not None:
            assert str(type(self.transform_resize.__dict__["transforms"][
                                0])) == "<class 'albumentations.augmentations.geometric.resize.Resize'>", "The first transform should be a resize transform"
            assert len(self.transform_resize.__dict__["transforms"]) == 1, "Only one resize transform should be applied"
            self.resize = int(self.transform_resize.__dict__["transforms"][0].__dict__['width'])

        self.transform_dapi = transform_dapi
        self.augment_rna_density = augment_rna_density

        if self.dict_gene_value is not None:
            self.nb_channel_rna = len(self.dict_gene_value[list(self.dict_gene_value.keys())[0]])
        else:
            self.nb_channel_rna = 1

        self.test_return_background = test_return_background

        self.remove_cell_in_background_threshold = remove_cell_in_background_threshold

        # for reloding the dataset

        self.remove_nucleus_seg_from_bg = remove_nucleus_seg_from_bg

        self.addition_mode = addition_mode

        self.return_df = return_df

        if self.return_df:
            assert gene2index is not None, "gene2index should be provided if return_df is True"
            assert 0 not in gene2index.values(), "gene2index should not contain 0 as value, reserved for padding"
            self.nb_channel_rna = 3
            self.gene2index = gene2index

        self.augmentation_img = augmentation_img

        self.recompute_flow = recompute_flow

        self.use_cache = use_cache
        self.experiment_name = experiment_name
        self.patch_dir = patch_dir

        # data augmentation ######

        self.cellbound_transform = cellbound_transform
        self.transfrom_augment_resize = A.Compose([
            A.RandomScale(scale_limit=(-0.5, 0.5), interpolation=1, p=0.5, always_apply=None),
            A.PadIfNeeded(min_height=self.resize, min_width=self.resize, border_mode=cv2.BORDER_CONSTANT, value=0, p=1),
            A.CropNonEmptyMaskIfExists(height=self.resize,
                                       width=self.resize,
                                       p=1.0),
            # A.RandomResizedCrop(size=(self.resize , self.resize) , scale=(0.5, 1), ratio=(0.75, 1.33), p=0.5)
        ])

        self.training_mode = training_mode
        self.evaluation_mode = evaluation_mode

        if self.list_patch_index is None:

            nb_patches = len(self.st_segmentation.sdata[self.st_segmentation.shape_patch_key])
            self.list_patch_index = list(range(nb_patches))
            shape_segmentation_key = self.st_segmentation.key_cell_consistent_with_nuclei
            self.list_patch_index = self.st_segmentation.get_valid_patch(min_nb_cell=self.min_nb_cell_per_patch,
                                                                         list_path_index=self.list_patch_index,
                                                                         shape_segmentation_key=shape_segmentation_key,
                                                                         min_transcripts=self.min_transcripts)

            log.info(f"Number of valid patches: {len(self.list_patch_index)}")
            print(f"Number of valid patches: {len(self.list_patch_index)}")
        else:
            print(f"list_patch_index overwrite min_transcript and min_nb_cell_per_patch")

        if list_annotation_patches is not None:
            log.info(f"Number of valid patches before removing annotation: {len(self.list_patch_index)}")
            self.list_patch_index = list(set(list(self.list_patch_index)) - set(list(list_annotation_patches)))
            log.info(f"Number of valid patches after removing annotation: {len(self.list_patch_index)}")

        if self.training_mode:
            if patch_width is None:
                import re
                patch_width = re.search(r"sopa_patches_rna2seg_(\d+)_", shape_patch_key)
            if self.st_segmentation.density_threshold is None:
                self._set_threshold(
                    max_nb_crops=500,
                    kernel_size=9,
                    percentile_threshold=5,
                    shape=(patch_width, patch_width)
                )

    def __len__(self):
        return len(self.list_patch_index)

    def __getitem__(self, idx):
        """
        Retrieves a patch of spatial transcriptomics data based on the given index.

        :param idx: The index of the patch to retrieve.
        :type idx: int
        :return: A dictionary containing the following elements:

            - "img_cellbound": Tensor representing the cell boundary image.
            - "dapi": Tensor representing the DAPI-stained image for nuclear staining.
            - "rna_img": Tensor representing the spatial RNA expression image.
            - "mask_flow": Tensor representing the cellular flow field mask for segmentation.
            - "mask_gradient": Tensor representing the gradient mask for segmentation refinement.
            - "background" (optional): Tensor representing the background image if `test_return_background` is enabled.
            - "idx": Integer index of the patch.
            - "patch_index": Integer representing the patch identifier.
            - "bounds": List defining the spatial boundaries of the patch.
            - "segmentation_nuclei" (optional): Tensor representing the nuclear segmentation mask if available.
            - "list_gene" (optional): Tensor containing the list of detected genes if `return_df` is enabled.
            - "array_coord" (optional): Tensor containing spatial coordinates of detected transcripts
             if `return_df` is enabled.
        """

        patch_index = self.list_patch_index[idx]

        # check if cache exist ###################
        compute_all_patch = True

        if self.use_cache:
            folder_to_save = Path(self.patch_dir) / str(patch_index) / self.experiment_name
            try:
                img_cellbound = tifffile.imread(folder_to_save / RNA2segFiles.CELLBOUND)
                dapi = tifffile.imread(folder_to_save / RNA2segFiles.DAPI)
                if self.training_mode:
                    mask_flow = tifffile.imread(folder_to_save / RNA2segFiles.MASK_FLOW)
                    mask_gradient = tifffile.imread(folder_to_save / RNA2segFiles.GRADIENT)

                if self.return_flow:
                    background = tifffile.imread(folder_to_save / RNA2segFiles.BACKGROUND)
                if self.return_df:
                    list_gene = np.load(folder_to_save / RNA2segFiles.LIST_GENE)
                    array_coord = np.load(folder_to_save / RNA2segFiles.ARRAY_COORD)
                if self.key_nuclei_segmentation:
                    segmentation_nuclei = tifffile.imread(folder_to_save / RNA2segFiles.SEGMENTATION_NUCLEI)

                path_csv = Path(self.patch_dir) / str(patch_index)
                bounds = json.load(open(path_csv / RNA2segFiles.BOUNDS_FILE, "r"))
                patch_df = pd.read_csv(path_csv / RNA2segFiles.TRANSCRIPTS_FILE)

                bounds = bounds["bounds"]
                rna_img = rna2img(
                    df_crop=patch_df,
                    dict_gene_value=self.dict_gene_value,
                    image_shape=(dapi.shape[1], dapi.shape[2], self.nb_channel_rna),
                    gene_column=self.gene_column,
                    column_x="x", column_y="y", offset_x=bounds[0], offset_y=bounds[1],
                    gaussian_kernel_size=self.kernel_size_rna2img,
                    max_filter_size=self.max_filter_size_rna2img,
                    addition_mode=self.addition_mode
                ).transpose(2, 0, 1)

                compute_all_patch = False

            except Exception as e:

                print(f"No cache found for patch {patch_index} Recomputing the patch {patch_index}")

        if compute_all_patch:
            dapi, rna_img, img_cellbound, mask_flow, mask_gradient, background, list_gene, \
                array_coord, segmentation_nuclei, bounds = self._get_patch_input(
                    patch_index=patch_index,
                )

        # data augmentation ###################################

        if self.augmentation_img:
            # hard coded value
            self.prob_rotate_resize = 0.5
            self.prob_transfrom_augment_resize = 0.75

            dapi, rna_img, img_cellbound, mask_flow, mask_gradient, background = self._augment_input(
                dapi=dapi,
                rna_img=rna_img,
                img_cellbound=img_cellbound,
                mask_flow=mask_flow,
                mask_gradient=mask_gradient,
                background=background if self.test_return_background else None
            )

        ############### create dict_result ###################################

        dict_result = {}
        dict_result["img_cellbound"] = torch.tensor(img_cellbound)
        dict_result["dapi"] = torch.tensor(dapi).clone().detach()

        dict_result["mask_flow"] = torch.tensor(mask_flow.astype(np.float32))
        dict_result["mask_gradient"] = torch.tensor(mask_gradient.astype(np.float32))

        if self.test_return_background:
            transformed_background = self.transform_resize(image=background)
            background = transformed_background["image"]
            dict_result["background"] = torch.tensor(background)
        dict_result["idx"] = idx
        dict_result["patch_index"] = patch_index
        dict_result["bounds"] = bounds
        if segmentation_nuclei is not None:
            dict_result["segmentation_nuclei"] = torch.tensor(segmentation_nuclei.astype(np.float32))

        if self.return_df:
            dict_result["list_gene"] = torch.tensor(list_gene)
            dict_result["array_coord"] = torch.tensor(array_coord)
        else:
            dict_result["rna_img"] = torch.tensor(rna_img)

        return dict_result

    def _get_patch_input(self, patch_index):

        if self.training_mode:
            self.return_agg_segmentation = True
            self.return_flow = True

        else:
            self.return_agg_segmentation = False
            self.return_flow = False

        if self.evaluation_mode:
            self.return_agg_segmentation = True
            self.return_flow = True

        (dapi, patch_df, agreement_segmentation, agreement_segmentation_without_nuclei,
         segmentation_nuclei, segmentation_nuclei_not_in_cell, bounds) = self.st_segmentation.get_patch_input(
            patch_index=patch_index,
            return_agg_segmentation=self.training_mode,
        )

        width, height = dapi.shape[0], dapi.shape[1]
        target_width = bounds[2] - bounds[0]
        target_height = bounds[3] - bounds[1]
        assert target_height == target_width, (f"width of the patch does not match height {target_height} "
                                               f"!= {target_width}, not supported yet")
        if (width, height) != (target_width, target_height):
            dapi = pad_image2D_to_shape(dapi, target_shape=(target_width, target_height))
            if agreement_segmentation is not None:
                agreement_segmentation = pad_image2D_to_shape(agreement_segmentation,
                                                              target_shape=(target_width, target_height))
                log.info(f"patch {patch_index} is not square, adding padding")

        if agreement_segmentation is None:
            agreement_segmentation = np.zeros((target_width, target_height))

        if self.st_segmentation.channels_cellbound is not None:
            img_cellbound = self.st_segmentation.get_cellbound_staining(patch_index)
            img_cellbound = img_cellbound.astype('float32')
            # normalize
            if img_cellbound.ndim == 2:
                img_cellbound = tf_cp.normalize99(img_cellbound)
            else:
                assert img_cellbound.ndim == 3
                for i in range(len(img_cellbound)):
                    img_cellbound = tf_cp.normalize99(img_cellbound)
            # pad image if necessary
            if (img_cellbound.shape[-2], img_cellbound.shape[-1]) != (target_width, target_height):
                new_img_cellbound = np.zeros((len(img_cellbound), target_width, target_height))
                for clb in range(len(img_cellbound)):
                    new_img_cellbound[clb] = pad_image2D_to_shape(img_cellbound[clb],
                                                                  target_shape=(target_width, target_height))
                img_cellbound = new_img_cellbound
        else:
            img_cellbound = None

        assert len(dapi.shape) == 2, "only one channel (dapi) is allowed for now"
        assert (patch_df["y"] >= bounds[1]).all(), "y coordinate outside the image"
        assert (patch_df["y"] <= bounds[3]).all(), "y coordinate outside the image"
        assert (patch_df["x"] >= bounds[0]).all(), "x coordinate outside the image"
        assert (patch_df["x"] <= bounds[2]).all(), "x coordinate outside the image"

        if self.transform_dapi is not None:
            dapi = self.transform_dapi(image=dapi)["image"]
        dapi = tf_cp.normalize99(dapi)

        if len(patch_df) == 0:
            rna = np.zeros((target_width, target_height, self.nb_channel_rna))
        else:

            if self.return_df:

                list_gene = list(patch_df[self.gene_column].values)
                list_gene = [self.gene2index[gene] for gene in list_gene]

                image_shape = (dapi.shape[0], dapi.shape[1], self.nb_channel_rna)
                scaling_factor_coord = self.transform_resize.__dict__["transforms"][0].__dict__['width'] / dapi.shape[
                    -1]
                array_coord = compute_array_coord(bounds,
                                                  patch_df,
                                                  image_shape,
                                                  scaling_factor_coord)
                rna = np.zeros((target_width, target_height, self.nb_channel_rna))
            else:
                list_gene = None
                array_coord = None
                if self.augment_rna_density:
                    print("/!\ Adding the augmentation on RNA density")
                    patch_df = create_augmented_patch_df(patch_df, factor=2, std_dev=5)

                rna = rna2img(df_crop=patch_df, dict_gene_value=self.dict_gene_value,
                              image_shape=(dapi.shape[0], dapi.shape[1], self.nb_channel_rna),
                              gene_column=self.gene_column,
                              column_x="x", column_y="y", offset_x=bounds[0], offset_y=bounds[1],
                              gaussian_kernel_size=self.kernel_size_rna2img,
                              max_filter_size=self.max_filter_size_rna2img,
                              addition_mode=self.addition_mode)

        img_input = np.concatenate([dapi[:, :, None], rna], axis=2).astype(np.float32)

        if self.return_flow:
            assert self.st_segmentation.density_threshold is not None, \
                "density threshold should be provided or set with self._set_threshold() "
            background = get_background_mask(density_threshold=self.st_segmentation.density_threshold,
                                             df_crop=patch_df,
                                             shape=dapi.shape,
                                             kernel_size=self.kernel_size_background_density,
                                             column_y="y",
                                             column_x="x",
                                             offset_x=bounds[0],
                                             offset_y=bounds[1],
                                             )

            # check padding
            if (agreement_segmentation.shape[0], agreement_segmentation.shape[1]) != (target_width, target_height):
                agreement_segmentation = pad_image2D_to_shape(agreement_segmentation,
                                                              target_shape=(target_width, target_height))

            if segmentation_nuclei is not None:
                if (segmentation_nuclei.shape[0], segmentation_nuclei.shape[1]) != (target_width, target_height):
                    segmentation_nuclei = pad_image2D_to_shape(segmentation_nuclei,
                                                               target_shape=(target_width, target_height))

            (mask_gradient, agreement_segmentation, background) = remove_cell_in_background(
                agreement_segmentation=agreement_segmentation,
                background=background,
                threshold_cell_in_bg=self.remove_cell_in_background_threshold,
                agreement_segmentation_without_nuclei=agreement_segmentation_without_nuclei,
                segmentation_nuclei_not_in_cell=segmentation_nuclei_not_in_cell,
                remove_nucleus_seg_from_bg=self.remove_nucleus_seg_from_bg,
                segmentation_nuclei=segmentation_nuclei,

            )
        else:
            mask_gradient = np.zeros((target_width, target_height))
            background = np.zeros((target_width, target_height))

        if self.return_flow:
            from cellpose.dynamics import labels_to_flows
            target = labels_to_flows(
                [agreement_segmentation], files=None, redo_flows=False, device=torch.device("cpu")
            )[0]
            target[:, agreement_segmentation == 0] = 0
            target[:, background > 0] = 0
            label, pred, flow_x, flow_y = target
        else:
            label, pred, flow_x, flow_y = [agreement_segmentation,
                                           np.zeros(agreement_segmentation.shape),
                                           np.zeros(agreement_segmentation.shape),
                                           np.zeros(agreement_segmentation.shape),
                                           ]
        if self.evaluation_mode:
            target = [label, pred, flow_x, flow_y]
        else:
            target = [pred, flow_x, flow_y]

        target = np.array(target)
        # check padding for target
        if (target.shape[1], target.shape[2]) != (target_width, target_height):
            target = pad_image3D_to_shape(target, target_shape=(target_width, target_height))
        len_target = len(target)

        # apply transformation
        if self.transform_resize:
            # todo : problem mask_flow is not the same size depending of the mode should be
            #  simplified by returning the label outside mask_flow
            # do fonction resize
            mask_flow = np.array([*target, mask_gradient]).astype(np.float32)
            if img_cellbound is not None:
                img_cellbound = img_cellbound.astype(np.float32)
                mask_flow = np.concatenate([mask_flow, img_cellbound], axis=0)
            transformed = self.transform_resize(image=img_input, masks=mask_flow)
            img_input = transformed['image'].astype(np.float32)
            mask_flow = np.array(transformed['masks'][:len_target]).astype(np.float32)
            mask_gradient = (np.array(transformed['masks'][len_target])).astype(int)
            img_cellbound = np.array(transformed['masks'][len_target + 1:])
        else:
            raise NotImplementedError("transform should be provided")
        rna2seg_input = np.transpose(img_input, (2, 0, 1))

        # save in a cache for future use ##################
        dapi = rna2seg_input[:1]
        rna_img = rna2seg_input[1:]

        if self.use_cache:

            folder_to_save = Path(self.patch_dir) / str(patch_index) / self.experiment_name
            folder_to_save.mkdir(exist_ok=True, parents=False)

            tifffile.imwrite(folder_to_save / RNA2segFiles.CELLBOUND, img_cellbound)
            tifffile.imwrite(folder_to_save / RNA2segFiles.DAPI, dapi)
            tifffile.imwrite(folder_to_save / RNA2segFiles.RNA_img, rna_img)

            # save bounds as json
            with open(folder_to_save / RNA2segFiles.BOUNDS_FILE, 'w') as f:
                json.dump(bounds, f)

            if self.training_mode:
                if mask_flow.shape[0] == 4:
                    tifffile.imwrite(folder_to_save / RNA2segFiles.MASK_FLOW, mask_flow[1:].astype(float))
                elif mask_flow.shape[0] == 3:
                    tifffile.imwrite(folder_to_save / RNA2segFiles.MASK_FLOW, mask_flow[:].astype(float))

                tifffile.imwrite(folder_to_save / RNA2segFiles.GRADIENT, mask_gradient)
            tifffile.imwrite(folder_to_save / RNA2segFiles.CELLBOUND, img_cellbound)

            if self.return_flow:
                tifffile.imwrite(folder_to_save / RNA2segFiles.BACKGROUND, background)

            if self.return_df:
                np.save(folder_to_save / RNA2segFiles.LIST_GENE, list_gene)
                np.save(folder_to_save / RNA2segFiles.ARRAY_COORD, array_coord)
            if segmentation_nuclei is not None:
                tifffile.imwrite(folder_to_save / RNA2segFiles.SEGMENTATION_NUCLEI, segmentation_nuclei)

        return (dapi, rna_img, img_cellbound, mask_flow, mask_gradient, background, list_gene, array_coord,
                segmentation_nuclei, bounds)

    def _augment_input(self,
                       dapi,
                       rna_img,
                       img_cellbound,
                       mask_flow,
                       mask_gradient,
                       background=None
                       ):

        assert dapi.ndim == 3 and dapi.shape[0] == 1, f"dapi should have shape (1, H, W) but is {dapi.shape}"

        dapi = tf_cp.normalize99(dapi[0])
        dapi = self.cellbound_transform(image=dapi)['image']  # Augment
        dapi = dapi[None, :, :]

        img_cellbound = img_cellbound.astype('float32')
        if img_cellbound.ndim == 2:  # for compatibility
            img_cellbound = np.array([img_cellbound])
        assert img_cellbound.ndim == 3
        if np.max(img_cellbound) > 3:  # check that it is not already normailzed between 0 and 1
            for i in range(len(img_cellbound)):
                img_cellbound[i] = tf_cp.normalize99(img_cellbound[i])
        augmented_image = self.cellbound_transform(image=np.transpose(img_cellbound, (1, 2, 0)))['image']
        img_cellbound = np.transpose(augmented_image, (2, 0, 1))

        # rotation #
        if np.random.rand() < self.prob_rotate_resize:

            n_dapi = dapi.shape[0]
            n_rna_img = rna_img.shape[0]
            ncb = img_cellbound.shape[0]
            input_image = np.concatenate((dapi, rna_img, img_cellbound, [mask_gradient]), axis=0)
            if self.test_return_background:
                input_image = np.concatenate((input_image, [background]), axis=0)

            transformed_input_image, mask_flow, scale = random_rotate_and_resize(
                input_image=input_image, mask_flow=mask_flow)

            dapi = transformed_input_image[:n_dapi]
            rna_img = transformed_input_image[n_dapi:n_dapi + n_rna_img]
            img_cellbound = transformed_input_image[n_dapi + n_rna_img:n_dapi + n_rna_img + ncb]
            mask_gradient = transformed_input_image[n_dapi + n_rna_img + ncb]
            if self.test_return_background:
                background = transformed_input_image[-1]

        # augment resize ###

        # should I leave RNA in transfrom_augment_resize ?
        if np.random.rand() < self.prob_transfrom_augment_resize:

            n_dapi = dapi.shape[0]
            n_rna_img = rna_img.shape[0]
            ncb = img_cellbound.shape[0]
            input_image = np.concatenate((dapi,
                                          rna_img,
                                          img_cellbound,
                                          [mask_gradient]), axis=0)
            if self.test_return_background:
                input_image = np.concatenate((input_image, [background]), axis=0)
            original_shape = input_image.shape
            original_shape_mask = mask_flow.shape

            dict_aug_resize = self.transfrom_augment_resize(
                image=np.transpose(input_image, (1, 2, 0)),
                mask=np.transpose(mask_flow, (1, 2, 0))
            )

            transformed_input_image = dict_aug_resize['image']
            transformed_input_image = np.transpose(transformed_input_image, (2, 0, 1))
            assert transformed_input_image.shape == original_shape, \
                f"transformed_input_image.shape should be {original_shape} but is {transformed_input_image.shape}"

            mask_flow = dict_aug_resize['mask']
            mask_flow = np.transpose(mask_flow, (2, 0, 1))
            assert mask_flow.shape == original_shape_mask, \
                f"original_shape_mask.shape should be {original_shape_mask} but is {original_shape_mask.shape}"

            dapi = transformed_input_image[:n_dapi]
            rna_img = transformed_input_image[n_dapi:n_dapi + n_rna_img]
            img_cellbound = transformed_input_image[n_dapi + n_rna_img:n_dapi + n_rna_img + ncb]
            mask_gradient = transformed_input_image[n_dapi + n_rna_img + ncb]
        return dapi, rna_img, img_cellbound, mask_flow, mask_gradient, background

    def _set_threshold(self, max_nb_crops=2000,
                       kernel_size=9,
                       percentile_threshold=5,
                       shape=(1200, 1200), ):
        shape_segmentation_key = self.st_segmentation.key_cell_consistent_with_nuclei
        list_path_index_theshold = self.st_segmentation.get_valid_patch(min_nb_cell=1,
                                                                        shape_segmentation_key=shape_segmentation_key,
                                                                        min_transcripts=10)

        if len(list_path_index_theshold) > max_nb_crops:
            list_path_index_theshold = np.random.choice(list_path_index_theshold, max_nb_crops, replace=False)

        t = time()
        # compute background density threshold
        print("compute density threshold")
        density_threshold = self.st_segmentation.set_density_threshold_sequential(
            list_path_index=list_path_index_theshold,
            shape_segmentation_key=shape_segmentation_key,
            shape=shape,
            kernel_size=kernel_size,
            percentile_threshold=percentile_threshold,
        )
        print(f"Time to compute density threshold: {time() - t:.6f}s")
        print(f"Density threshold: {density_threshold}")

        self.st_segmentation.density_threshold = density_threshold

    def get_rna_img(
            self,
            bounds,
            key_transcripts="transcripts",
            dict_gene_value=None,
    ):
        """
        Generates an image of RNA from spatial transcriptomics data within the specified bounds.

        :param bounds: The bounding box coordinates for the extracted region (xmin, ymin, xmax, ymax).
        :type bounds: tuple[int, int, int, int]
        :param key_transcripts: The key to access transcriptomic data in the dataset. Defaults to "transcripts".
        :type key_transcripts: str
        :param dict_gene_value: Dictionary mapping gene names to encodings/colors. If None, all genes have a value of 1.
        :type dict_gene_value: dict[str, float] | None
        :return: An image representation of RNA expression in the selected region.
        :rtype: np.ndarray
        """
        print("Get RNA image ...")
        sdata = self.st_segmentation.sdata
        df = sdata[key_transcripts]
        img = sdata[self.st_segmentation.image_key]
        df = to_intrinsic(sdata, df, img)
        patch_df = df[
            (df["x"] >= bounds[0]) & (df["x"] <= bounds[2]) &
            (df["y"] >= bounds[1]) & (df["y"] <= bounds[3])
            ]
        patch_df = patch_df.compute()
        if dict_gene_value is None:
            dict_gene_value = self.dict_gene_value

        img = rna2img(
            patch_df,
            dict_gene_value=dict_gene_value,
            image_shape=(bounds[3] - bounds[1], bounds[2] - bounds[0], 3),
            offset_x=bounds[0], offset_y=bounds[1], gene_column=self.gene_column,
            max_filter_size=5,
            addition_mode=True
        )
        return img

    def get_staining_img(self, bounds):
        """
        Generates an image of dapi and cell boundaries stainings within the specified bounds. 

        :param bounds: The bounding box coordinates for the extracted region (xmin, ymin, xmax, ymax).
        :type bounds: tuple[int, int, int, int]
        :return: An image of the different staining (channels_dapi and channels_cellbound) in the selected region.
        :rtype: np.ndarray
        """
        print("Get image ...")
        xmin, ymin, xmax, ymax = bounds
        st_segmentation = self.st_segmentation
        stainings = [*st_segmentation.channels_dapi, *st_segmentation.channels_cellbound]
        image = st_segmentation.image.sel(
            c=stainings,
            x=slice(xmin, xmax),
            y=slice(ymin, ymax),
        ).values
        return image

    def get_segmentation_img(self, bounds, key_cell, image=None, color="red", size_line=5, min_size=5):
        """
        Generates an image of cell segmentation within the specified bounds.

        :param bounds: The bounding box coordinates for the extracted region (xmin, ymin, xmax, ymax).
        :type bounds: tuple[int, int, int, int]
        :param key_cell: The key to access cell segmentation data.
        :type key_cell: str
        :param image: The base image on which segmentation will be overlaid. If None, uses DAPI staining.
        :type image: np.ndarray | None
        :param color: The color used to outline cell contours. Defaults to "red".
        :type color: str
        :param size_line: The thickness of the segmentation contour lines. Defaults to 5.
        :type size_line: int
        :return: An image with cell segmentation overlaid.
        :rtype: np.ndarray
        """
        print("Get segmentation image ...")

        self.st_segmentation.sdata[key_cell] = to_intrinsic(
            self.st_segmentation.sdata, self.st_segmentation.sdata[key_cell], self.st_segmentation.sdata[self.st_segmentation.image_key])

        # Define image on which to plot the segmentation (dapi or cb)
        if image is None:
            image = self.get_staining_img(bounds)[0]  # Dapi staining if image is None
        image_with_mask = image / image.max()
        if image_with_mask.ndim == 2:
            image_with_mask = np.stack([image_with_mask] * 3, axis=-1)

        # Get segmentation
        segmentation = self.st_segmentation.get_segmentation_crop(
            bounds=bounds, shape=image_with_mask.shape[0:2], key_cell=key_cell)
        cell_contour_mask = create_cell_contours(segmentation, min_size=min_size, size_line=size_line)

        image_with_mask[cell_contour_mask > 0] = mcolors.hex2color(color)

        return image_with_mask


def custom_collate_fn(batch):
    batched_data = {}
    # Assuming all dictionaries have the same structure
    keys = batch[0].keys()

    max_cb_channel = 4
    for key in keys:
        items = [d[key] for d in batch]
        if key in ['img_cellbound']:
            items = [d[key].clone().detach() for d in batch]
            max_cb_channel = max([item.shape[0] for item in items])
            zero_img = torch.zeros_like(items[0][0:1])
            items = [torch.cat([item, zero_img.repeat(max_cb_channel - item.shape[0], 1, 1)], dim=0) for item in items]

        if key in ["list_gene", "array_coord"]:
            items_padded = torch.nn.utils.rnn.pad_sequence(items, batch_first=True, padding_value=0)
            batched_data[key] = items_padded
        else:
            # For fixed-size data or non-tensor data, simply stack or collect
            if isinstance(items[0], torch.Tensor):
                batched_data[key] = torch.stack(items)
            else:
                items = np.array(items)
                batched_data[key] = torch.tensor(items)

    return batched_data
