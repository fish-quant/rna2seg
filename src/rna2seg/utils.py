from tqdm import tqdm
from pathlib import Path
import geopandas as gpd
from sopa.segmentation import shapes
from spatialdata import SpatialData
from spatialdata._io import write_shapes
from spatialdata.models import ShapesModel


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


def save_shapes2zarr(dataset, segmentation_key):
    list_all_cells = load_segmentation2zarr(path_parquet_files=dataset.patch_dir)
    print(f"len(list_all_cells) {len(list_all_cells)}")
    unique_cells = shapes.solve_conflicts(
        cells = list_all_cells,
        threshold = 0.25,
        patch_indices = None,
        return_indices = False,
    )

    sdata = dataset.st_segmentation.sdata
    geo_df = gpd.GeoDataFrame({"geometry": unique_cells})
    geo_df = ShapesModel.parse(geo_df)

    segmentation_key = f"rna2seg_{segmentation_key}"
    sdata[segmentation_key] = geo_df
    sdata.write_element(segmentation_key)

    print(f"Added {len(geo_df)} cell boundaries in sdata['{segmentation_key}']")




