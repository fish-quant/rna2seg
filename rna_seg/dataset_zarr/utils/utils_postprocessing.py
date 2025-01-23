


import geopandas as gpd
import shapely
import spatialdata as sd
from sopa._sdata import get_spatial_image, save_shapes
from spatialdata.models import ShapesModel
from spatialdata.transformations import get_transformation
from tqdm import tqdm
import numpy as np
import tifffile

from  sopa.segmentation import shapes

# sopa.segmentation import shapes

path_zarr = "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/mouse_ileum/RNAseg_DATASET.zarr"
path_seg_cell = "/cluster/CBIO/data1/st_segmentation/data_release_baysor_merfish_gut/raw_dataset/data_analysis/cellpose/cell_boundaries/results/cellpose_membrane.tif"
path_seg_nuclei = "/cluster/CBIO/data1/st_segmentation/data_release_baysor_merfish_gut/raw_dataset/data_analysis/cellpose/cell_boundaries/results/cellpose_dapi.tif"



if __name__ == "__main__":

    ## mask to shape

    # Load the data
    sdata = sd.read_zarr(path_zarr)

    image_key = list(sdata.images.keys())[0]


    # Define the segmentation method
    seg_mask = tifffile.imread(path_seg_cell)
    seg_mask = seg_mask[3]
    cells = shapes.geometrize(seg_mask)
    geo_df = gpd.GeoDataFrame(geometry=cells)
    geo_df.index = image_key + geo_df.index.astype(str)


    image = get_spatial_image(sdata, image_key)
    geo_df = ShapesModel.parse(
        geo_df, transformations=get_transformation(image, get_all=True).copy()
    )
    shapes_key = "cellpose_cell_paper"
    sdata.shapes[shapes_key] = geo_df
    save_shapes(sdata, shapes_key, overwrite=True)


    seg_mask = tifffile.imread(path_seg_nuclei)
    for i in range(len(seg_mask)):
        seg_mask_z_i = seg_mask[3]
        cells = shapes.geometrize(seg_mask_z_i)
        geo_df = gpd.GeoDataFrame(geometry=cells)
        geo_df.index = image_key + geo_df.index.astype(str)

        image = get_spatial_image(sdata, image_key)
        geo_df = ShapesModel.parse(
            geo_df, transformations=get_transformation(image, get_all=True).copy()
        )
        shapes_key = f"z_{i}_cellpose_nuclei_paper_z"
        sdata.shapes[shapes_key] = geo_df
        save_shapes(sdata, shapes_key, overwrite=True)









