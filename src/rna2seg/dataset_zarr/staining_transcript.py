#### tile the images
from __future__ import annotations

import logging
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import rasterize
from scipy.ndimage import gaussian_filter
from shapely import Polygon
from skimage import exposure
from sopa._constants import SopaFiles
from sopa.segmentation._stainings import StainingSegmentation
from spatialdata import SpatialData
from tqdm import tqdm

from rna2seg.dataset_zarr.consistency import compute_polygon_intersection

log = logging.getLogger(__name__)

import dask

dask.config.set({'dataframe.query-planning': False})
dask.config.set(scheduler='synchronous')


class StainingTranscriptSegmentation(StainingSegmentation):
    def __init__(
            self,
            sdata: SpatialData,
            channels_dapi: list[str] | str,
            channels_cellbound: list[str] | str,
            shape_patch_key: str,
            image_key: str | None = None,
            key_cell_segmentation: str | None = None,
            key_nuclei_segmentation: str | None = None,
            key_cell_consistent: str | None = None,
            key_nucleus_consistent: str | None = None,
            density_threshold: float | None = None,
            patch_dir: str | Path | None = None,
            min_area: float = 0,
            clip_limit: float = 0,
            gaussian_sigma: float = 0,
    ):
        """
        :param sdata: SpatialData
        :type SpatialData
        :param channels_dapi: list of channel for the DAPI staining
        :type list
        :channels_cellbound: list of channel for the cellbound staining
        :type list
        :shape_patch_key: name of the shape containing the patch in the sdata
        :type str
        :param image_key: key of the image in the sdata
        :type str
        :param key_cell_segmentation: key of the cell segmentation in the sdata NOT USE TO REMOVE
        :type str
        :param key_cell_consistent: key of the consistent cell segmentation in the sdata
        :type str
        :param key_cell_consistent_without_nuclei: key of the consistent cell segmentation without nuclei in the sdata
         CAN BE DIRECTLY INCLUDED IN key_cell_consistent
        :type str
        :param key_nuclei_segmentation: key of the nuclei segmentation in the sdata
        :type str
        :param key_nuclei_consistent_not_in_cell: key of the nuclei segmentation not in cell in the sdata
        :type str
        """

        super().__init__(
            sdata=sdata,
            method=lambda x: x,  # dummy
            channels=channels_dapi,
            image_key=image_key,
            min_area=min_area,
            clip_limit=clip_limit,
            gaussian_sigma=gaussian_sigma
        )

        self.sdata = sdata
        del self.channels
        self.channels_dapi = channels_dapi
        self.channels_cellbound = channels_cellbound

        self.key_cell_segmentation = key_cell_segmentation  # to remove ?
        self.key_nuclei_segmentation = key_nuclei_segmentation

        self.key_cell_consistent = key_cell_consistent
        if self.key_cell_consistent is not None:
            self.key_cell_consistent_with_nuclei = f'{key_cell_consistent}_with_nuclei'
            self.key_cell_consistent_without_nuclei = f'{key_cell_consistent}_without_nuclei'
        else:
            self.key_cell_consistent_with_nuclei = None
            self.key_cell_consistent_without_nuclei = None

        self.key_nucleus_consistent = key_nucleus_consistent
        if self.key_nucleus_consistent is not None:
            self.key_nucleus_consistent_not_in_cell = f'{key_nucleus_consistent}_not_in_cell'
            self.key_nucleus_consistent_in_cell = f'{key_nucleus_consistent}_in_cell'
        else:
            self.key_nucleus_consistent_not_in_cell = None
            self.key_nucleus_consistent_in_cell = None

        for key in [key_cell_segmentation, self.key_cell_consistent_with_nuclei,
                    self.key_cell_consistent_without_nuclei,
                    self.key_nucleus_consistent_not_in_cell,
                    self.key_nuclei_segmentation]:
            if key is not None:
                assert key in self.sdata.shapes, f"{key} not in sdata.shapes {self.sdata.shapes}"

        if self.key_cell_consistent is not None:  # to remove
            self.cell_consistent_with_nuclei = self.sdata[self.key_cell_consistent_with_nuclei]
            self.cell_consistent_without_nuclei = self.sdata[self.key_cell_consistent_without_nuclei]

        self.density_threshold = density_threshold
        self.patch_dir = patch_dir
        if self.patch_dir is None:
            self.patch_dir = Path(self.sdata.path) / ".rna2seg"
        self.shape_patch_key = shape_patch_key

        self.y_max = sdata[self.image_key]["scale0"].dims['y']
        self.x_max = sdata[self.image_key]["scale0"].dims['x']

    def get_patch_input(self,
                        patch_index: int = None,
                        patch: Polygon | None = None,
                        return_agg_segmentation: bool = False):

        assert patch_index is not None or patch is not None, "either patch_index XOR patch should be provided"
        assert patch_index is None or patch is None, "either patch_index XOR patch should be provided"

        if patch_index is not None:
            patch = self.sdata[self.shape_patch_key].geometry[patch_index]

        bounds = [int(x) for x in patch.bounds]

        if self.channels_dapi:
            image = self.image.sel(
                c=self.channels_dapi,
                x=slice(bounds[0], bounds[2]),
                y=slice(bounds[1], bounds[3]),
            ).values

            assert image.shape[0] == 1, "only one channel (DAPI) is allowed for now"
            image = image[0]

            if self.gaussian_sigma > 0:
                image = gaussian_filter(image, sigma=self.gaussian_sigma)
            if self.clip_limit > 0:  # https://en.wikipedia.org/wiki/Adaptive_histogram_equalization#CLAHE
                image = exposure.equalize_adapthist(image, clip_limit=self.clip_limit)

        else:
            image = np.zeros((bounds[3] - bounds[1], bounds[2] - bounds[0]), dtype=np.uint16)

        patch_df = pd.read_csv(Path(self.patch_dir) / f'{patch_index}/{SopaFiles.TRANSCRIPTS_FILE}')

        if return_agg_segmentation:
            assert self.key_cell_consistent is not None, "aggrement segmentation is not loaded/ initialized"
            #t = time()
            consistent_seg_cell = self.get_segmentation_crop(
                bounds=bounds,
                shape=image.shape,
                key_cell=self.key_cell_consistent_with_nuclei)
            if self.key_cell_consistent_without_nuclei is not None:  # to remove
                consistent_seg_cell_without_nuclei = self.get_segmentation_crop(
                    bounds=bounds,
                    shape=image.shape,
                    key_cell=self.key_cell_consistent_without_nuclei)
            else:
                consistent_seg_cell_without_nuclei = None
        else:
            consistent_seg_cell = None
            consistent_seg_cell_without_nuclei = None

        if self.key_nuclei_segmentation is not None:
            segmentation_nuclei = self.get_segmentation_crop(
                bounds=bounds,
                shape=image.shape,
                key_cell=self.key_nuclei_segmentation)
        else:
            segmentation_nuclei = None

        if self.key_nucleus_consistent_not_in_cell is not None:  # to remove
            segmentation_nuclei_not_in_cell = self.get_segmentation_crop(
                bounds=bounds,
                shape=image.shape,
                key_cell=self.key_nucleus_consistent_not_in_cell)
        else:
            segmentation_nuclei_not_in_cell = None

        return (image,
                patch_df,
                consistent_seg_cell,
                consistent_seg_cell_without_nuclei,
                segmentation_nuclei,
                segmentation_nuclei_not_in_cell,
                bounds)

    def get_cellbound_staining(self,
                               patch_index: int, ):

        patch = self.sdata[self.shape_patch_key].geometry[patch_index]
        bounds = [int(x) for x in patch.bounds]
        #t = time()

        image = self.image.sel(
            c=self.channels_cellbound,
            x=slice(bounds[0], bounds[2]),
            y=slice(bounds[1], bounds[3]),
        ).values
        return image

    def set_density_threshold_sequential(self,
                                         list_path_index: list[int] | None = None,
                                         shape_segmentation_key: str = "cellpose_boundaries",
                                         shape: tuple = (1200, 1200),
                                         kernel_size=5,
                                         percentile_threshold=1,
                                         ):

        """
        Compute the density threshold for a list of patch
        :param list_path_index:
        :param shape_segmentation_key:
        :param patch_dir:
        :param shape:
        :param kernel_size: kernel size for the density mask estiamtion with gaussina filter
        :param percentile_threshold:  use as example :  np.percentile(all_list_density, percentile_threshold = 7)
                    which means that the 7 percentile of the density distribution of the instance object is used as threshold
        :return:
        """

        from rna2seg.dataset_zarr.background import (
            get_mean_density_per_polygon, get_rna_density)

        if list_path_index is None:
            list_path_index = list(range(len(self.sdata[self.shape_patch_key].geometry)))

        segmentation_shapes = self.sdata[shape_segmentation_key]
        list_density = []
        for patch_index in tqdm(list_path_index):
            patch = self.sdata[self.shape_patch_key].geometry[patch_index]
            bounds = [int(x) for x in patch.bounds]
            gdf_polygon_segmentation = segmentation_shapes.cx[bounds[0]:bounds[2], bounds[1]:bounds[3]]
            patch_df = pd.read_csv(Path(self.patch_dir) / f'{patch_index}/{SopaFiles.TRANSCRIPTS_FILE}')

            if len(patch_df) == 0 or len(gdf_polygon_segmentation) == 0:
                continue

            density_mask = get_rna_density(df_crop=patch_df,
                                           shape=shape,
                                           kernel_size=kernel_size,
                                           column_x="x",
                                           column_y="y",
                                           offset_y=bounds[1],
                                           offset_x=bounds[0],
                                           )

            list_density += get_mean_density_per_polygon(
                gdf_polygon=gdf_polygon_segmentation,
                density_mask=density_mask,
                shape=shape,
                x_trans=bounds[0],
                y_trans=bounds[1],
            )
        density_threshold = np.percentile(list_density, percentile_threshold)
        self.density_threshold = density_threshold
        self.list_density = list_density
        return density_threshold

    def get_valid_patch(self,
                        min_nb_cell=1,
                        list_path_index=None,
                        shape_segmentation_key=None,
                        check_border=True,
                        min_transcripts=1):

        if shape_segmentation_key is None:
            shape_segmentation_key = self.key_cell_consistent_with_nuclei

        if list_path_index is None:
            list_path_index = list(range(len(self.sdata[self.shape_patch_key].geometry)))

        non_empty_list_path_index = []
        # compute valid patches
        for patch_index in tqdm(list_path_index, file=sys.stdout, desc="Get valid patches"):
            patch = self.sdata[self.shape_patch_key].geometry[patch_index]
            bounds = [int(x) for x in patch.bounds]
            if check_border:
                if self.y_max < bounds[3] or self.x_max < bounds[2]:
                    continue

            if min_transcripts > 0:
                patch_df = pd.read_csv(Path(self.patch_dir) / f'{patch_index}/{SopaFiles.TRANSCRIPTS_FILE}')
                if len(patch_df) < min_transcripts:
                    continue

            if shape_segmentation_key is not None:
                segmentation_shapes = self.sdata[shape_segmentation_key]
                gdf_polygon_segmentation = segmentation_shapes.cx[bounds[0]:bounds[2], bounds[1]:bounds[3]]
                if len(gdf_polygon_segmentation) < min_nb_cell:
                    continue
            non_empty_list_path_index.append(patch_index)
        return non_empty_list_path_index

    def get_patch_annotation(self,
                             polygon_annotation: gpd.GeoDataFrame | list[Polygon]
                             ):
        """
        Get the patch index that intersect with the annotation
        :param polygon_annotation:
        :return: path index that intersect with the annotation and that should be removed from the training set
        """
        list_polygon_patch = list(
            self.sdata[self.shape_patch_key].geometry)  # can be optimise see read_patches_cells sopa function ?

        if isinstance(polygon_annotation, gpd.GeoDataFrame):
            polygon_annotation = list(polygon_annotation.geometry)
        else:
            assert isinstance(polygon_annotation,
                              list), "polygon_annotation should be a list of Polygon or a GeoDataFrame"

        intersection_p = compute_polygon_intersection(list_polygon_patch, polygon_annotation)
        return intersection_p[1]

    # function for pre-computing all the flow to train cellpose

    def get_segmentation_crop(self, bounds, shape, key_cell):

        assert key_cell in self.sdata.shapes, f"{key_cell} not in sdata.shapes {self.sdata.shapes}"
        cell_segmentation = self.sdata[key_cell].cx[bounds[0]:bounds[2], bounds[1]:bounds[3]]
        """
        ## can be optimised with R-tree index ? if needed ?
        sindex = self.image.sindex
        possible_matches_index = list(sindex.intersection(bounds))
        possible_matches = self.image.iloc[possible_matches_index]
        precise_matches = possible_matches.cx[bounds[0]:bounds[2], bounds[1]:bounds[3]]
        """
        try:
            polygons = list(cell_segmentation.geometry)
            ## rasterize the polygon
            # This is typically an affine transformation, but for simplicity, we'll use an identity transformation
            x_trans = bounds[0]
            y_trans = bounds[1]
            transform = rasterio.Affine(a=1, b=0, c=x_trans, d=0, e=1, f=y_trans)

            # Define the shape of the output numpy array
            cell_segmentation = rasterize(
                ((polygons[i], i + 1) for i in range(len(polygons))),
                out_shape=shape, transform=transform,
                fill=0, all_touched=True, dtype=np.uint16
            )
        except ValueError as e:
            #print(e)
            assert len(polygons) == 0, "no polygon in the consistent_seg_cell"
            return np.zeros(shape)
        return cell_segmentation
