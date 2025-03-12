from __future__ import annotations

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.features import rasterize
from scipy.ndimage import gaussian_filter


def get_rna_density(df_crop,
                    shape,
                    kernel_size=5,
                    column_y="y",
                    column_x="x",
                    offset_x=0,
                    offset_y=0,
                    ):
    """
    Get the density of rna in the image
    :param df_crop:
    :param shape:
    :param kernel_size:
    :param column_x:
    :param column_y:
    :return:
    """

    df_crop = df_crop.copy()
    df_crop[column_y] = df_crop[column_y].astype(int) - offset_y
    df_crop[column_x] = df_crop[column_x].astype(int) - offset_x
    # drop the points outside the image
    df_crop = df_crop[(df_crop[column_y] >= 0) & (df_crop[column_y] < shape[0])]
    df_crop = df_crop[(df_crop[column_x] >= 0) & (df_crop[column_x] < shape[1])]
    list_y = df_crop[column_y].values
    list_x = df_crop[column_x].values
    coords = np.array(list(zip(list_y, list_x)))

    img_rna_count = np.zeros(shape, dtype=np.float32)
    img_rna_count[coords[:, 0], coords[:, 1]] += 1

    density_mask = gaussian_filter(img_rna_count, sigma=kernel_size)
    return density_mask


# compute the min density per instance object in a mask
def get_mean_density_per_polygon(gdf_polygon: gpd.GeoDataFrame,
                                 density_mask: np.array,
                                 shape: tuple,
                                 x_trans: float,
                                 y_trans: float,
                                 ):
    """
    Compute the list of density per instance object in a mask
    Args:
        mask_instance: mask with instance object
        density_mask: density mask

    Returns:
        min_density_per_instance: min density per instance object
    """
    polygons = list(gdf_polygon.geometry)
    # rasterize the polygon
    # This is typically an affine transformation, but for simplicity, we'll use an identity transformation
    transform = rasterio.Affine(a=1, b=0, c=x_trans, d=0, e=1, f=y_trans)
    # Define the shape of the output numpy array
    mask_instance = rasterize(((polygons[i], i + 1) for i in range(len(polygons))), out_shape=shape,
                              transform=transform, fill=0, all_touched=True, dtype=np.uint16)

    list_density = []
    unique_instance = np.unique(mask_instance)
    for instance in unique_instance:
        if instance == 0:
            continue
        density = density_mask[mask_instance == instance].mean()
        list_density.append(density)
    return list_density


def get_background_mask(density_threshold: float,
                        df_crop,
                        shape,
                        kernel_size=5,
                        column_y="y",
                        column_x="x",
                        offset_x=0,
                        offset_y=0,
                        ):
    """
    Get the background mask
    Args:
        density_mask: density mask
        density_threshold: density threshold

    Returns:
        background_mask: background mask
    """
    if len(df_crop) > 0:
        density_mask = get_rna_density(
            df_crop=df_crop,
            shape=shape,
            kernel_size=kernel_size,
            column_y=column_y,
            column_x=column_x,
            offset_x=offset_x,
            offset_y=offset_y,
        )
    else:
        density_mask = np.zeros(shape)

    background_mask = (density_mask < density_threshold).astype(np.uint8)
    return background_mask
