
#### tile the images
from __future__ import annotations
import json
import shutil
import logging
import geopandas as gpd
import spatialdata as sd
from pathlib import Path
from functools import partial
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
from sopa._constants import SopaFiles, SopaKeys
from sopa._sdata import to_intrinsic
from sopa.segmentation import Patches2D
from sopa.patches.patches import TranscriptPatches

from rna2seg._constant import RNA2segFiles

import dask
dask.config.set({'dataframe.query-planning': False})



log = logging.getLogger(__name__)


def create_patch_rna2seg(sdata : sd.SpatialData,
                        image_key : str,
                        points_key : str,
                        patch_width : int, patch_overlap : int,
                        min_transcripts_per_patch : int,
                        folder_patch_rna2seg : Path | str | None = None,
                        overwrite : bool = True):

    """
    Creates patches from the spatial data to handle data in manageable sizes.
    Save the patches shapes into the zarr and precomputes transcript.csv files for each patch.

    :param sdata: SpatialData object containing the spatial dataset.
    :type SpatialData:
    :param image_key: Key identifying the image in sdata.
    :type str:
    :param points_key: Key identifying the points in sdata.
    :type str:
    :param patch_width: Width of each patch.
    :type int:
    :param patch_overlap: Overlap between adjacent patches.
    :type int:
    :param min_transcripts_per_patch: Minimum number of transcripts required for a patch.
    :type int:
    :param folder_patch_rna2seg: Directory where patches will be saved. If None, defaults to sdata.path/.rna2seg.
    :type Path | str | None:
    :param overwrite: Whether to overwrite the existing patches shape if it already exists in the zarr.
    :type bool:
    """

    if folder_patch_rna2seg is None:
        folder_patch_rna2seg = Path(sdata.path) / ".rna2seg"

    coordinate_system = f'_{image_key}_intrinsic'
    for element in sdata._gen_spatial_element_values():
        try:
            sd.transformations.operations.remove_transformation(element, coordinate_system)
        except KeyError:
            pass
    patches = Patches2D(sdata=sdata,
                        element_name=image_key,
                        patch_width=patch_width,
                        patch_overlap=patch_overlap)
    # save polygons patch in the sdata
    shape_patch_key = f"sopa_patches_rna2seg_{patch_width}_{patch_overlap}"
    if not overwrite:
        if not Path(sdata.path / f"shapes/{shape_patch_key}").exists():
            patches.write(shapes_key = f"sopa_patches_rna2seg_{patch_width}_{patch_overlap}",)
    else:
        if  Path(sdata.path / f"shapes/{shape_patch_key}").exists():
            shutil.rmtree(sdata.path / f"shapes/{shape_patch_key}")
        patches.write(shapes_key = f"sopa_patches_rna2seg_{patch_width}_{patch_overlap}",)


    if not overwrite:
        if Path(folder_patch_rna2seg).exists():
            raise ValueError(f"folder {folder_patch_rna2seg} already exists, set overwrite to True")
    csv_name = SopaFiles.TRANSCRIPTS_FILE
    # save a 'scaled' rna-csv  for each patch in the folder
    tp = TranscriptPatches_with_scale(
        sdata=sdata,
        patches_2d=patches,
        df=sdata[points_key],
        config_name="",
        csv_name=csv_name,
        min_transcripts_per_patch=min_transcripts_per_patch
    )

    tp.write_image_scale(
        temp_dir = folder_patch_rna2seg,
        cell_key = None,
        unassigned_value = None,
        image_key=image_key,
    )

class TranscriptPatches_with_scale(TranscriptPatches):
    """
    A class to handle the patching of transcript segmentation at the image scale
    """
    def __init__(
            self,
            sdata : sd.SpatialData,
            patches_2d: Patches2D | list,
            df: dd.DataFrame | gpd.GeoDataFrame,
            config_name: str,
            csv_name: str,
            min_transcripts_per_patch: int,
    ):
        self.patches_2d = patches_2d
        self.df = df
        self.min_transcripts_per_patch = min_transcripts_per_patch # to remove not use
        self.config_name = config_name
        self.csv_name = csv_name
        self.sdata = sdata #self.patches_2d.sdata if self.patches_2d is not None else None

    def write_image_scale(
        self,
        temp_dir: str,
        cell_key: str = None,
        unassigned_value: int | str = None,
        image_key : str = "image",
    ):
        """
        Write a sub-CSV for transcript segmentation for all patches at the images scale
        Args:
            temp_dir:
            cell_key:
            unassigned_value:
            use_prior:
            image_key:
            list_polygon: list of polygons for each patch to write arbitrary polygone localization

        Returns:

        """

        log.info("Writing sub-CSV for transcript segmentation")

        self.temp_dir = Path(temp_dir)

        if cell_key is None:
            cell_key = SopaKeys.DEFAULT_CELL_KEY

        if unassigned_value is not None and unassigned_value != 0:
            self.df[cell_key] = self.df[cell_key].replace(unassigned_value, 0)

        self.df = to_intrinsic(self.sdata, self.df, image_key)
        self._setup_patches_directory()
        if isinstance(self.patches_2d, list):
            patches_gdf = gpd.GeoDataFrame(geometry=self.patches_2d)
        else:
            assert isinstance(self.patches_2d, Patches2D), "patches_2d must be a Patches2D object or a list of polygons"
            patches_gdf = gpd.GeoDataFrame(geometry=self.patches_2d.polygons)
        if isinstance(self.df, dd.DataFrame):
            with ProgressBar():
                self.df.map_partitions(
                    partial(self._query_points_partition, patches_gdf), meta=()
                ).compute()
        else:
            self._write_points(patches_gdf, self.df)
        self._write_path_bound()
        return list(self.valid_indices())

    def _write_path_bound(self, temp_dir: str=None):

        if temp_dir is None:
            assert "temp_dir" in self.__dict__
        else:
            if 'temp_dir' not in self.__dict__:
                self.temp_dir = Path(temp_dir)
            else:
                assert self.temp_dir == Path(temp_dir), f"temp_dir is not the same as the one already set {self.temp_dir} != {temp_dir}"

        assert self.temp_dir.exists(), f"temp_dir {self.temp_dir} does not exist, save first csv file"

        patches_gdf = gpd.GeoDataFrame(geometry=self.patches_2d.polygons)
        for index, polygon in patches_gdf.iterrows() :


            dict2json = {"bounds" : list(polygon.geometry.bounds),
                         "bounds_min_x" : polygon.geometry.bounds[0],
                         "bounds_min_y" : polygon.geometry.bounds[1],
                         "bounds_max_x" : polygon.geometry.bounds[2],
                         "bounds_max_y" : polygon.geometry.bounds[3],
                         }
            path2save_json = self.temp_dir / f'{str(index)}/{RNA2segFiles.BOUNDS_FILE}'
            # assert not path2save_json.exists(), f"json file {path2save_json} already exists, overwriting is not allowed"
            with open(path2save_json, "w") as f:
                json.dump(dict2json, f)
