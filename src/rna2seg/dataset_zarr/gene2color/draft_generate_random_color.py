





from tifffile import imread
import tifffile
import pandas as pd

from pathlib import Path

import json
from skimage.measure import find_contours

import rasterio
from rasterio.features import rasterize
import numpy as np


from shapely.geometry import Polygon

import geopandas as gpd
from rna2seg.dataset_zarr.gene2color.utils import get_gene_random_vector



if __name__ == "__main__":

    gene2color_ref = np.load("/home/tom/share/st3/open_vizgen/RNAseg_DATASET/mouse_ileum/0710_dict_color_sample500/gene2vect_dist_mean_sigma_standard.npy",
                         allow_pickle=True).item()


    gene_list = list(gene2color_ref.keys())


    path_save_rd_dict = "/home/tom/share/st3/open_vizgen/RNAseg_DATASET/mouse_ileum/random 3color"

    for i_dict in range(30):
        gene2color = get_gene_random_vector(gene_list, nb_pc=3)
        np.save(path_save_rd_dict + f"/gene2color_rd_rgb_color_{i_dict}.npy", gene2color)

    gene2color_binary = {gene: np.array( [1]) for gene in gene2color}
    np.save(f"/home/tom/share/st3/open_vizgen/RNAseg_DATASET/mouse_ileum/gene2color_rd_binary_color_{i_dict}.npy", gene2color_binary)

