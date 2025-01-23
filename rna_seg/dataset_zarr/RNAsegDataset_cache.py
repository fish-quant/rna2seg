#%%


import tifffile
from pathlib import Path
from tqdm import tqdm
import pandas as pd
import numpy as np
import json
import torch
import sys
from rna_seg.dataset_zarr.RNAsegDataset import compute_array_coord, rna2img
from rna_seg.dataset_zarr.staining_transcript import StainingTranscriptSegmentation
from rna_seg.dataset_zarr.data_augmentation import (augment_embeddings,
                                                    random_rotate_and_resize,
                                                    cellbound_transform,
                                                    cellbound_transform2,
                                                    cellbound_transform3)
from rna_seg._constant import RNAsegFiles
from cellpose import transforms as tf_cp
import albumentations as A
from albumentations.core.transforms_interface import ImageOnlyTransform
#from albumentations import Compose, GaussNoise
from spatialdata import SpatialData
import spatialdata as sd
import logging
log = logging.getLogger(__name__)
from rna_seg.dataset_zarr.gene2color.utils import get_gene_random_vector
import cv2

from rna_seg._constant import RNAsegFiles


class RNAsegDatasetCache():

    def __init__(self,
                 path_to_cache: Path | str,
                 experiment_folder: str,
                 key_cell_consistent: str,
                 min_nb_cell_per_patch = 0,
                 list_path_index = None,
                 return_flow = True,
                 return_flow_omnipose : bool = False,
                 test_return_background = False,
                 test_mode = False,
                 list_annotation_patches = None,
                 return_df = False,
                 #todo param to recompute input from rna2img
                 recompute_rna2img = False,
                 gene_column = "gene",
                 gene2color = None,
                 original_image_shape = (1200, 1200),
                 kernel_size_rna2img = 0,
                 max_filter_size_rna2img = 0,
                 transform = None,
                 gene2index = None,
                 dapi_only_mode = False, # for ablative study
                augmentation_img = False,
                augmentation_rna = False,
                augmentation_version = 1,
                 prob_rotate_resize = 0.5,
                 use_random_color = False,
                 addition_mode = True,
                 ):


        # assert return_df == False or recompute_rna2img == False, "return_df and generate_rna2img cannot be True at the same time"
        self.path_to_cache = Path(path_to_cache)
        assert self.path_to_cache.exists(), f"{self.path_to_cache} does not exist"
        self.experiment_folder = experiment_folder
        self.key_cell_consistent = key_cell_consistent
        self.min_nb_cell_per_patch = min_nb_cell_per_patch
        if list_path_index is None:
            self.set_valid_indices()
        else:
            self.list_path_index = list_path_index
        if list_annotation_patches is not None:
            self.list_annotation_patches = [self.path_to_cache / str(idx) for idx in list_annotation_patches]
            ## remove the annotation patches from the list_path_index
            self.list_path_index = [idx for idx in self.list_path_index if idx not in self.list_annotation_patches]
            assert len(self.list_path_index) == len(set(self.list_path_index) - set(list_annotation_patches))
            print("annotation patches removed from the list_path_index")
            log.info("annotation patches removed from the list_path_index")
        self.return_flow = return_flow
        self.return_flow_omnipose = return_flow_omnipose
        self.test_return_background = test_return_background
        self.test_mode = test_mode
        ## todo :init the param to recompute input
        self.return_df = return_df
        # param to recompute input from rna2img
        self.recompute_rna2img = recompute_rna2img
        self.gene_column = gene_column
        self.gene2color = gene2color
        self.gene2index = gene2index
        if gene2index is not None and gene2color is not None:
            assert sorted(list(gene2index.keys())) == sorted(list(gene2index.keys())), "Genes in gene2index and gene2color should be the same."
        self.original_image_shape = original_image_shape
        self.kernel_size_rna2img = kernel_size_rna2img
        self.max_filter_size_rna2img = max_filter_size_rna2img
        self.transform = transform
        if self.transform is not None:
            assert str(type(self.transform.__dict__["transforms"][0])) == "<class 'albumentations.augmentations.geometric.resize.Resize'>",  "The first transform should be a resize transform"
            assert len(self.transform.__dict__["transforms"]) == 1, "Only the resize transform should is supported"
            self.resize = int(self.transform.__dict__["transforms"][0].__dict__['width'])

        self.dapi_only_mode = dapi_only_mode
        self.augmentation_img = augmentation_img
        self.augmentation_rna   = augmentation_rna
        self.prob_rotate_resize = prob_rotate_resize
        self.use_random_color = use_random_color
        assert augmentation_version in [1, 2, 3, 4, 5], "augmentation_version should be 1, 2, 3 or 4"
        self.augmentation_version = augmentation_version
        if self.augmentation_version == 2:
            self.cellbound_transform = cellbound_transform2
        elif self.augmentation_version == 1:
                self.cellbound_transform = cellbound_transform
        elif  self.augmentation_version in [3, 4, 5]:
            if self.augmentation_version in [3, 4]:
                self.cellbound_transform = cellbound_transform2
            elif self.augmentation_version == 5:
                self.cellbound_transform = cellbound_transform3

            self.transfrom_resize = A.Compose([
                A.RandomScale(scale_limit=(-0.5, 0.5), interpolation=1, p=0.5, always_apply=None),
                A.PadIfNeeded(min_height=self.resize , min_width=self.resize , border_mode=cv2.BORDER_CONSTANT, value=0, p=1),
                A.CropNonEmptyMaskIfExists(height=self.resize , width=self.resize , p=1.0),
                A.RandomResizedCrop(height=self.resize , width=self.resize , scale=(0.5, 1), ratio=(0.75, 1.33), p=0.5)
            ])

        self.addition_mode = addition_mode


    def __len__(self):
        return len(self.list_path_index)

    def __getitem__(self, idx):


        path_index = self.list_path_index[idx]

        rna_seg_input = self._get_rna_seg_input(path_index)
        mask_flow = self._get_mask_flow(path_index)
        img_cellbound = self._get_img_cellbound(path_index)
        mask_gradient = self._get_mask_gradient(path_index)
        if self.test_return_background:
            background = self._get_background(path_index)
        
        if self.augmentation_img:
            # Rotation and Resize on Input Image and Target
            if np.random.rand() < self.prob_rotate_resize:
                
                nchan = rna_seg_input.shape[0]
                ncb = img_cellbound.shape[0]
                input_image = np.concatenate((rna_seg_input, img_cellbound, [mask_gradient]), axis=0)
                if self.test_return_background:
                    input_image = np.concatenate((input_image, [background]), axis=0)

                transformed_input_image, mask_flow, scale = random_rotate_and_resize(
                    input_image=input_image, mask_flow=mask_flow)
                
                rna_seg_input = transformed_input_image[:nchan]
                img_cellbound = transformed_input_image[nchan:nchan+ncb]
                mask_gradient = transformed_input_image[nchan+ncb]
                if self.test_return_background:
                    background = transformed_input_image[-1]
        if self.augmentation_version in [3, 4, 5]:
            #print("augmentation version 3")

            nchan = rna_seg_input.shape[0]
            ncb = img_cellbound.shape[0]
            input_image = np.concatenate((rna_seg_input, img_cellbound, [mask_gradient]), axis=0)
            if self.test_return_background:
                input_image = np.concatenate((input_image, [background]), axis=0)
            original_shape = input_image.shape
            original_shape_mask = mask_flow.shape

            dict_aug_resize = self.transfrom_resize(
                    image=np.transpose(input_image, (1, 2, 0))#, mask=mask_flow)
                    ,mask=np.transpose(mask_flow, (1, 2, 0)))
            transformed_input_image = dict_aug_resize['image']
            transformed_input_image = np.transpose(transformed_input_image, (2, 0, 1))

            assert transformed_input_image.shape == original_shape,\
                f"transformed_input_image.shape should be {original_shape} but is {transformed_input_image.shape}"
            mask_flow = dict_aug_resize['mask']
            mask_flow = np.transpose(mask_flow, (2, 0, 1))
            assert mask_flow.shape == original_shape_mask,\
                f"original_shape_mask.shape should be {original_shape_mask} but is {original_shape_mask.shape}"
            rna_seg_input = transformed_input_image[:nchan]
            img_cellbound = transformed_input_image[nchan:nchan+ncb]
            mask_gradient = transformed_input_image[nchan+ncb]


        dict_result = {}
        if self.test_mode:
            dict_result["idx"] = idx
        dict_result["rna_seg_input"] = rna_seg_input # to remove
        dict_result["dapi"] = rna_seg_input[:1]
        dict_result["rna_img"] = rna_seg_input[1:]
        dict_result["mask_flow"] = mask_flow
        dict_result["mask_gradient"] = mask_gradient
        dict_result["img_cellbound"] = img_cellbound 
        if self.test_return_background:
            dict_result["background"] = background

        dict_result["patch_index"] = int(self.list_path_index[idx].stem)
        if self.return_df:
            try:
                list_gene = np.load(path_index / f"{self.experiment_folder}/{RNAsegFiles.LIST_GENE}")
                array_coord = np.load(path_index / f"{self.experiment_folder}/{RNAsegFiles.ARRAY_COORD}")
                dict_result["list_gene"] = torch.tensor(np.array(list_gene))
                dict_result["array_coord"] = torch.tensor(np.array(array_coord))
            except FileNotFoundError as e:
                #log.info(f"Error in loading list_gene and array_coord : {e}")
                #log.info('Load the df_crop and compute the list_gene and array_coord')
                patch_df = pd.read_csv(path_index / f"{RNAsegFiles.TRANSCRIPTS_FILE}")
                with open(path_index / f"{RNAsegFiles.BOUNDS_FILE}", "r") as f:
                    dict_bounds = json.load(f)
                list_gene = list(patch_df[self.gene_column].values)
                list_gene = [self.gene2index[gene] for gene in list_gene]
                ## create array of coordinates
                scaling_factor_coord = self.transform.__dict__["transforms"][0].__dict__['width'] / self.original_image_shape[0]
                array_coord = compute_array_coord(dict_bounds["bounds"],
                                                  patch_df=patch_df,
                                                  image_shape= self.original_image_shape,
                                                  scaling_factor_coord=scaling_factor_coord,)
                dict_result["list_gene"] = torch.tensor(np.array(list_gene))
                dict_result["array_coord"] = torch.tensor(np.array(array_coord))

        return dict_result

    def _get_rna_seg_input(self, path_index):

        rna_seg_input = tifffile.imread(path_index / f"{self.experiment_folder}/{RNAsegFiles.RNA_SEG_INPUT}")
        
        if self.recompute_rna2img:
            patch_df = pd.read_csv(path_index / f"{RNAsegFiles.TRANSCRIPTS_FILE}")
            patch_df = patch_df[~patch_df[self.gene_column].isna()].reset_index(drop=True)
            with open(path_index / f"{RNAsegFiles.BOUNDS_FILE}", "r") as f:
                dict_bounds = json.load(f)
                if self.augmentation_rna:
                    raise NotImplementedError("augmentation_rna is not implemented")
                    gene2color = augment_embeddings(self.gene2color,
                                                    noise_level=0.1,
                                                    min_scale_factor=0.9,
                                                    max_scale_factor=1.1,
                                                    min_shift=-1,
                                                    max_shift=1,
                                                    noise_p=0.5,
                                                    shift_p=0.5)
                else:
                    gene2color = self.gene2color

                if self.use_random_color:
                    gene2color = get_gene_random_vector(list(self.gene2color.keys()), nb_pc=3)

                rna = rna2img(df_crop=patch_df, dict_gene_value=gene2color, image_shape=(
                self.original_image_shape[0], self.original_image_shape[1], len(list(self.gene2color.values())[0])),
                              gene_column=self.gene_column, column_x="x", column_y="y",
                              offset_x=dict_bounds['bounds'][0], offset_y=dict_bounds['bounds'][1],
                              gaussian_kernel_size=self.kernel_size_rna2img,
                              max_filter_size=self.max_filter_size_rna2img, addition_mode=self.addition_mode)
                if self.dapi_only_mode:
                    rna[:] = 0
                rna = self.transform(image=rna)['image'].astype(np.float32)

                dapi = rna_seg_input[0]
                dapi_image = self.transform(image=dapi)['image']  # Resize
                if self.augmentation_img:
                    transformed_dapi  = self.cellbound_transform(image=dapi_image)['image'] # Augment
                else:
                    transformed_dapi = dapi

                rna_seg_input = np.concatenate([transformed_dapi[:,:,None], rna], axis=2).astype(np.float32)
                rna_seg_input = np.transpose(rna_seg_input, (2, 0, 1))

        return rna_seg_input

    def _get_mask_flow(self, path_index):
        if self.return_flow:
            if self.return_flow_omnipose:
                mask_flow = tifffile.imread(Path(path_index / f"{self.experiment_folder}/{RNAsegFiles.MASK_FLOW_OMNIPOSE}"))
                assert mask_flow.shape[0] == 6, f"mask_flow.shape[0] should be 6 but is {mask_flow.shape[0]}"
            else:
                mask_flow = tifffile.imread(path_index / f"{self.experiment_folder}/{RNAsegFiles.MASK_FLOW}")
                assert mask_flow.shape[0] == 3, f"mask_flow.shape[0] should be 3 but is {mask_flow.shape[0]}"
            
            transformed_mask_flow = self.transform(image=np.transpose(mask_flow, (1, 2, 0)))['image']
            mask_flow = np.transpose(transformed_mask_flow, (2, 0, 1)).astype(np.float32)
        
        else:
            mask_flow = np.zeros((3, self.original_image_shape[0], self.original_image_shape[1]))
    
        return mask_flow

    def _get_img_cellbound(self, path_index):
        img_cellbound = tifffile.imread(path_index / f"{self.experiment_folder}/{RNAsegFiles.CELLBOUND}")
        img_cellbound = img_cellbound.astype('float32')
        if img_cellbound.ndim == 2: # for compatibility
            img_cellbound = np.array([img_cellbound])
        assert img_cellbound.ndim  == 3
        if np.max(img_cellbound) > 3: # check that it is not already normailzed between 0 and 1
            for i in range(len(img_cellbound)):
                img_cellbound[i] = tf_cp.normalize99(img_cellbound[i])
        if self.augmentation_img:
            # Data Augmentation on Cellbound Images
            augmented_image = self.cellbound_transform(image=np.transpose(img_cellbound, (1, 2, 0)))['image']
            img_cellbound = np.transpose(augmented_image, (2, 0, 1))

        transformed_img_cellbound = self.transform(image=np.transpose(img_cellbound, (1, 2, 0)))['image'].astype(np.float32)
        img_cellbound = np.transpose(transformed_img_cellbound, (2, 0, 1)).astype(np.float32)

        return img_cellbound
    
    def _get_mask_gradient(self, path_index):
        mask_gradient = tifffile.imread(
            path_index / f"{self.experiment_folder}/{RNAsegFiles.GRADIENT}")
        return self.transform(image=mask_gradient)['image'].astype(np.float32)
    
    def _get_background(self, path_index):
            background = tifffile.imread(
                path_index / f"{self.experiment_folder}/{RNAsegFiles.BACKGROUND}"
            )
            return self.transform(image=background)['image']
    
    def set_valid_indices(self):

        if (self.path_to_cache / f"0/{self.experiment_folder}.npy").exists():
            list_path_index = np.load(self.path_to_cache / f"0/{self.experiment_folder}.npy",
                                        allow_pickle=True)
            log.info(f" valid index loaded from {self.path_to_cache / f'0/{self.experiment_folder}/.npy'}")
            print(f" valid index loaded from {self.path_to_cache / f'0/{self.experiment_folder}/.npy'}")
            print(f"valid len(list_path_index) : {len(list_path_index)}")

            ## checking valid index
            index_are_valid = True
            if not "open_vizgen" in str(self.path_to_cache):
                for path_idx in list_path_index:
                    if not Path(path_idx / f"{self.experiment_folder}/{RNAsegFiles.RNA_SEG_INPUT}").exists():
                        index_are_valid = False
                        break
            if len(list_path_index)>0 and index_are_valid:
                self.list_path_index = list_path_index
                print(f"path checked, using precomuted path for cache")
                return
        log.info(f"Compute the valid index")
        print(f"Compute the valid index")
        list_path_index = []

        for path_idx in tqdm(list(self.path_to_cache.glob("*")), file = sys.stdout, desc="valid index"):
            if not path_idx.is_dir():
                raise FileNotFoundError(f"{path_idx} is not a directory")
            if not (path_idx / f"{self.experiment_folder}").is_dir():
                log.info(f"{path_idx / f'{self.experiment_folder}'} does not exist check the preprocessing")
                print(f"{path_idx / f'{self.experiment_folder}'} does not exist check the preprocessing")
                #continue
                raise FileNotFoundError(f"{path_idx / f'{self.experiment_folder}'} does not exist")
            patch_df = pd.read_csv(path_idx / f"{RNAsegFiles.TRANSCRIPTS_FILE}")
            if patch_df.empty:
                continue

            if (path_idx / f"{self.experiment_folder}/{RNAsegFiles.NB_CELL_FILE}").exists():
                with open(path_idx / f"{self.experiment_folder}/{RNAsegFiles.NB_CELL_FILE}", "rb") as f:
                    nb_cells = int(f.readline())
                    ## read the first line
                    path = path_idx / f"{self.experiment_folder}/{RNAsegFiles.NB_CELL_FILE}"
                    #log.info(f"{path} exists")
                    #print(f"{path} exists")
            elif (path_idx / f"{self.key_cell_consistent}/{RNAsegFiles.NB_CELL_FILE}").exists():
                with open(path_idx / f"{self.key_cell_consistent}/{RNAsegFiles.NB_CELL_FILE}", "rb") as f:
                    nb_cells = int(f.readline())
                    path = path_idx / f"{self.key_cell_consistent}/{RNAsegFiles.NB_CELL_FILE}"

                    #log.info(f"{path} exists")
                    #print(f"{path} exists")
                    ## read the first line
            else:
                log.info(f"Compute the number of cells in patch :{path_idx}")
                print(f"Compute the number of cells in patch :{path_idx}")
                nb_cells = len(np.unique(np.load(path_idx / f"{self.key_cell_consistent}/{RNAsegFiles.LABEL_AND_MASK_FLOW}")[0])) - 1
                ## write the number of cells in the first line
                with open(path_idx / f"{self.key_cell_consistent}/{RNAsegFiles.NB_CELL_FILE}", "wb") as f:
                    f.write(str(nb_cells).encode())
            if nb_cells >= self.min_nb_cell_per_patch:
                list_path_index.append(path_idx)

        def sort_key(path):
            return int(path.name)
        # Sort the list using the custom key
        list_path_index = sorted(list_path_index, key=sort_key)

        self.list_path_index = list_path_index

        np.save(self.path_to_cache / f"0/{self.experiment_folder}", np.array(list_path_index))







