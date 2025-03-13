# tile the images
from __future__ import annotations

import json
import logging
import shutil
from functools import partial
from pathlib import Path

import dask
import dask.dataframe as dd
import geopandas as gpd
import numpy as np
import pandas as pd
import sopa
import spatialdata as sd
from dask.diagnostics import ProgressBar
from sopa._constants import SopaFiles
from sopa.patches._transcripts import OnDiskTranscriptPatches
from sopa.utils.utils import (get_feature_key, to_intrinsic)

from rna2seg._constant import RNA2segFiles

dask.config.set({'dataframe.query-planning': False})

log = logging.getLogger(__name__)


def create_patch_rna2seg(sdata: sd.SpatialData,
                         image_key: str,
                         points_key: str,
                         patch_width: int, patch_overlap: int,
                         min_points_per_patch: int,
                         folder_patch_rna2seg: Path | str | None = None,
                         overwrite: bool = False,
                         gene_column_name: str = "gene",
                         ):
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
    :param min_points_per_patch: Minimum number of transcripts required for a patch.
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

    shape_patch_key = f"sopa_patches_rna2seg_{patch_width}_{patch_overlap}"

    if Path(sdata.path / f"shapes/{shape_patch_key}").exists():
        if not overwrite:
            raise ValueError(f"folder {folder_patch_rna2seg} already exists, set overwrite to True")

    sopa.make_image_patches(sdata,
                            patch_width=patch_width,
                            patch_overlap=patch_overlap,
                            image_key=image_key,
                            key_added=shape_patch_key)

    csv_name = SopaFiles.TRANSCRIPTS_FILE
    # save a 'scaled' rna-csv  for each patch in the folder
    tp = TranscriptPatchesWithScale(
        sdata=sdata,
        points_key=points_key,
        patch_width=patch_width,
        patch_overlap=patch_overlap,
        shape_patch_key=shape_patch_key,
        df=sdata[points_key],
        config_name="",
        csv_name=csv_name,
        min_points_per_patch=min_points_per_patch,
        cache_dir=folder_patch_rna2seg,
        gene_column_name=gene_column_name,

    )
    tp.bboxes = np.array(sdata[shape_patch_key].bboxes)

    tp.write_image_scale(
        image_key,
        shape_patch_key,
    )

    return tp


class TranscriptPatchesWithScale(OnDiskTranscriptPatches):
    """
    A class to handle the patching of transcript segmentation at the image scale
    """

    def __init__(
            self,
            sdata: sd.SpatialData,
            points_key: str,
            patch_width: float,
            patch_overlap: float,
            shape_patch_key: str,
            df: dd.DataFrame | gpd.GeoDataFrame,
            config_name: str,
            csv_name: str,
            min_points_per_patch: int,
            cache_dir: str,
            gene_column_name: str = "gene",
    ):

        super().__init__(sdata=sdata,
                         points_key=points_key,
                         patch_width=patch_width,
                         patch_overlap=patch_overlap,
                         )

        self.df = df
        self.min_points_per_patch = min_points_per_patch  # to remove not use
        self.config_name = config_name
        self.csv_name = csv_name
        self.sdata = sdata  # self.patches_2d.sdata if self.patches_2d is not None else None
        self.shape_patch_key = shape_patch_key
        self.cache_dir = Path(cache_dir)
        self.gene_column_name = gene_column_name

    def write_image_scale(
            self,
            image_key: str,
            shape_patch_key: str,
            intrinsics: bool = True,
    ):
        """
        Write a sub-CSV for transcript segmentation for all patches at the images scale
        Args:
            cache_dir:

        Returns:

        """

        self.setup_patches_directory()

        # line to changed from sopa 2
        patches_geo_df = self.sdata[shape_patch_key]

        gene_column = get_feature_key(self.points)

        # line to changed from sopa 2
        if intrinsics:
            df_with_scale = to_intrinsic(self.sdata, self.points, image_key)
        else:
            df_with_scale = self.points
        # keep only the columns x, y and gene
        # df_with_scale = df_with_scale[["x", "y", gene_column]]

        with ProgressBar():
            df_with_scale.map_partitions(
                partial(self.query_points_partition, patches_geo_df, gene_column=gene_column), meta=()
            ).compute()

        with ProgressBar():
            self.points.map_partitions(
                partial(self.query_points_partition, patches_geo_df, gene_column=gene_column), meta=()
            ).compute()
        self._write_patch_bound()
        return list(self.valid_indices())

    def query_points_partition(
            self, patches_gdf: gpd.GeoDataFrame, df: pd.DataFrame, gene_column: str | None = None
    ) -> pd.DataFrame:

        points_gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df["x"], df["y"]))
        self.write_points(patches_gdf, points_gdf, mode="a")

    def write_points(
            self, patches_gdf: gpd.GeoDataFrame, points_gdf: gpd.GeoDataFrame, mode="w", csv_name: str | None = None
    ):
        patches_gdf.index.name = "index_right"  # to reuse the index name later
        df_merged: gpd.GeoDataFrame = points_gdf.sjoin(patches_gdf)

        for index, patch_df in df_merged.groupby("index_right"):
            patch_path = self.get_patch_path(index, csv_name=csv_name)
            patch_path.parent.mkdir(parents=True, exist_ok=True)
            patch_df = patch_df.drop(columns=["index_right", "geometry"])
            # set columns name x, y and gene

            assert "x" in patch_df.columns, f"column 'x' not found in {patch_df.columns}"
            assert "y" in patch_df.columns, f"column 'y' not found in {patch_df.columns}"
            assert self.gene_column_name in patch_df.columns, (f"column {self.gene_column_name} not found in "
                                                               f"{patch_df.columns}")

            list_x = patch_df["x"].tolist()
            list_y = patch_df["y"].tolist()
            list_gene = patch_df[self.gene_column_name].tolist()
            # create a new dataframe
            patch_df_new = pd.DataFrame(list(zip(list_x, list_y, list_gene)), columns=['x', 'y', self.gene_column_name])
            # print(f'patch_df_new {patch_df_new}')
            # set columns x, y and gene
            patch_df_new.to_csv(patch_path)
            # patch_df_new.to_csv(patch_path, mode=mode, header=mode == "w", index=False)

    def _write_patch_bound(self, cache_dir: str = None):

        if cache_dir is None:
            assert "cache_dir" in self.__dict__
        else:
            if 'cache_dir' not in self.__dict__:
                self.cache_dir = Path(cache_dir)
            else:
                assert self.cache_dir == Path(
                    cache_dir), f"cache_dir is not the same as the one already set {self.cache_dir} != {cache_dir}"

        assert self.cache_dir.exists(), f"cache_dir {self.cache_dir} does not exist, save first csv file"

        patches_gdf = gpd.GeoDataFrame(geometry=self.sdata[self.shape_patch_key].geometry)
        for index, polygon in patches_gdf.iterrows():
            dict2json = {"bounds": list(polygon.geometry.bounds),
                         "bounds_min_x": polygon.geometry.bounds[0],
                         "bounds_min_y": polygon.geometry.bounds[1],
                         "bounds_max_x": polygon.geometry.bounds[2],
                         "bounds_max_y": polygon.geometry.bounds[3],
                         }
            path2save_json = self.cache_dir / f'{str(index)}/{RNA2segFiles.BOUNDS_FILE}'
            with open(path2save_json, "w") as f:
                json.dump(dict2json, f)
