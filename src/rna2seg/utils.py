from pathlib import Path

import geopandas as gpd
import numpy as np
from scipy import ndimage as ndi
from sopa.segmentation.resolve import solve_conflicts
from spatialdata import SpatialData
from spatialdata._io import write_shapes
from spatialdata.models import ShapesModel
from tqdm import tqdm


def create_cell_contours(seg_mask, size_line=5, min_size=2500):
    contour_mask = np.zeros(seg_mask.shape)
    cell_to_contour = np.unique(seg_mask)

    nb_cell_in = 0
    for cell in tqdm(cell_to_contour):
        mask = (seg_mask == cell).astype(int)
        if cell == 0 or mask.sum() < min_size:
            continue
        nb_cell_in += 1
        contour_mask += (ndi.maximum_filter(mask, size=size_line) - ndi.minimum_filter(mask, size=size_line))

    return contour_mask


# fonction from sopa 1.0.14
def save_shapes(
        sdata: SpatialData,
        name: str,
        overwrite: bool = False,
) -> None:
    if not sdata.is_backed():
        return

    elem_group = sdata._init_add_element(name=name, element_type="shapes", overwrite=overwrite)

    write_shapes(
        shapes=sdata.shapes[name],
        group=elem_group,
        name=name,
    )


def load_segmentation2zarr(path_parquet_files):
    list_all_cells = []
    for path in tqdm(list(Path(path_parquet_files).glob("*.parquet"))):
        gdf = gpd.read_parquet(path)
        list_all_cells += list(gdf.geometry)
    return list_all_cells


def save_shapes2zarr(dataset, path_parquet_files, segmentation_key, overwrite=False):
    """

    Args:
        path_parquet_files:
        segmentation_key:

    Returns:

    """
    list_all_cells = load_segmentation2zarr(path_parquet_files=path_parquet_files)
    print(f"len(list_all_cells) {len(list_all_cells)}")
    unique_cells = solve_conflicts(
        cells=list_all_cells,
        threshold=0.25,
        patch_indices=None,
        return_indices=False,
    )

    sdata = dataset.st_segmentation.sdata
    geo_df = gpd.GeoDataFrame({"geometry": unique_cells.geometry})
    geo_df = ShapesModel.parse(geo_df)

    # segmentation_key = f"rna2seg_{segmentation_key}"
    sdata[segmentation_key] = geo_df
    sdata.write_element(segmentation_key, overwrite=overwrite)

    print(f"Added {len(geo_df)} cell boundaries in sdata['{segmentation_key}']")
