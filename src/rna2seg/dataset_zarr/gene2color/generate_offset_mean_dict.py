





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



    path_all_dict = "/home/tom/share/st4/open_vizgen/RNAseg_DATASET/all_dictionary"

    offset_dist = 2
    for dataset_folder in Path(path_all_dict).iterdir():

        dict_mean_only = {}
        dict_binary_and_mean = {}
        print(dataset_folder)

        ## load unorm dataset

        gene2color_ref = np.load(dataset_folder / "0710_dict_color_sample500/gene2vect_dist_unorm_mean_sigma_standard.npy",
                      allow_pickle=True).item()
        gene_list = list(gene2color_ref.keys())

        mean_list = [gene2color_ref[m][0] for m in gene2color_ref]
        print([np.mean(mean_list), np.min(mean_list), np.max(mean_list)])

        for gene in gene_list:
            dict_mean_only[gene] = np.array([gene2color_ref[gene][0] + offset_dist])
            dict_binary_and_mean[gene] = np.array([1, gene2color_ref[gene][0] + offset_dist])



        ## save result
        np.save(dataset_folder / f"0710_dict_color_sample500/dict_mean_only_offset{offset_dist}.npy", dict_mean_only)
        np.save(dataset_folder / f"0710_dict_color_sample500/dict_binary_and_mean_offset{offset_dist}.npy", dict_binary_and_mean)






    path_all_dict = "/home/tom/share/st5/open_vizgen/RNAseg_DATASET/all_dictionary"

    offset_dist = 2
    for dataset_folder in Path(path_all_dict).iterdir():

        dict_blank = {}
        print(dataset_folder)

        ## load unorm dataset

        gene2color_ref = np.load(dataset_folder / "0710_dict_color_sample500/gene2vect_dist_unorm_mean_sigma_standard.npy",
                                 allow_pickle=True).item()
        gene_list = list(gene2color_ref.keys())

        for gene in gene_list:
            dict_blank[gene] = np.array([0])
            #dict_binary_and_mean[gene] = np.array([1, gene2color_ref[gene][0] + offset_dist])

        np.save(dataset_folder / f"0710_dict_color_sample500/dict_blank.npy", dict_blank)





    path_all_dict = "/home/tom/share/st3/open_vizgen/RNAseg_DATASET/all_dictionary"

    for dataset_folder in Path(path_all_dict).iterdir():

        dict_binary = {}
        print(dataset_folder)

        ## load unorm dataset

        gene2color_ref = np.load(dataset_folder / "0710_dict_color_sample500/gene2vect_dist_unorm_mean_sigma_standard.npy",
                                 allow_pickle=True).item()
        gene_list = list(gene2color_ref.keys())

        for gene in gene_list:
            dict_binary[gene] = np.array([1])
            #dict_binary_and_mean[gene] = np.array([1, gene2color_ref[gene][0] + offset_dist])

        np.save(dataset_folder / f"0710_dict_color_sample500/dict_binary.npy", dict_binary)