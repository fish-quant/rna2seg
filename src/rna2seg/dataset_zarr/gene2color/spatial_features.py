



#%%




import time
import numpy as np
import os
import pandas
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import tifffile

from tqdm import tqdm
from pathlib import Path
import json
from CNN_gene.rna_seg._constant import RNAsegFiles
import argparse
## import ndi
from scipy import ndimage as ndi
from datetime import datetime
import random


def compute_rna2nuclei_distance(nuclei_mask,
                                df_spots,
                                x_column_name="global_x",
                                y_column_name="global_x",
                                use_max_threshold = True,
                                max_threshold = 20,
                                pixel_um_size = 0.108,
                                ):
    """
    Args:
        nuclei_mask:  np.array of shape (z, y, x) containing the mask of the nuclei for each z
        df_spots: dataframe containing the spots with the columns x, y, z, gene
        plot_distribution:
    Returns:
        compute the "relative" mean and std distance of eahc RNA species to the nuclei boundary
        relative mean here negative if inside, positive if outside
        distance2nuclei_boundary_mean is a pandas series with the gene as index and the mean distance as value
        distance2nuclei_boundary_std is a pandas series with the gene as index and the std distance as value
    """
    assert nuclei_mask.ndim == 2
    mask_contour = ndi.maximum_filter(nuclei_mask, size=2) - ndi.minimum_filter(
        nuclei_mask, size=2)
    mask_contour_bool = mask_contour != 0
    inverted_mask = np.ones(mask_contour_bool.shape, dtype = np.uint8) - (mask_contour_bool != 0).astype(np.uint8)
    distance = ndi.distance_transform_edt(inverted_mask,
                                          sampling=pixel_um_size)
    if not use_max_threshold:
        max_threshold = np.inf
    f_sign = lambda x: -1 if x > 0 else 1
    list_distance2nuclei_boundary = []
    for index, spots in df_spots.iterrows():
        x = int(round(spots[x_column_name]))
        y = int(round(spots[y_column_name]))
        try:
            distance2nuclei_boundary = f_sign(nuclei_mask[y, x]) * distance[y, x]
            if distance2nuclei_boundary >= max_threshold and f_sign(nuclei_mask[y, x]) == 1:
                list_distance2nuclei_boundary.append(np.nan)
            else:
                list_distance2nuclei_boundary.append(distance2nuclei_boundary)
        except IndexError as e:
            assert x == nuclei_mask.shape[1] or y == nuclei_mask.shape[0]
            list_distance2nuclei_boundary.append(np.nan)
    df_spots["distance2nuclei_boundary"] = list_distance2nuclei_boundary

    return df_spots



def get_gene_nuclei_dist_vector(
        list_path_nuclei,
        list_path_df_spots,
        norm = "standard",
        x_column_name = "global_x",
        y_column_name = "global_x",
        #z_column_name = "global_z",
        dict_scale={"x": 1, "y": 1, "z": 1},
        use_max_threshold=True,
        max_threshold=15,
        gene_column = "gene",
        pixel_um_size = 1,
):

    """
    Args:
        nuclei_all_z:  np.array of shape (z, y, x) containing the mask of the nuclei for each z
        df_spots: dataframe containing the spots with the columns x, y, z, gene
        path_to_save: folder to save results (distance_all_z and gene2vect).
        norm: normalization to apply to the distance vectors. Can be min max or standardization.
        x_column_name: name of the column containing the x coordinate of the spots.
        y_column_name: name of the column containing the y coordinate of the spots.
        z_column_name: name of the column containing the z coordinate of the spots.
        dict_scale: dictionary containing the scaling factors in µm for the x, y, and z coordinates.
        use_max_threshold: boolean to use a max threshold for the distance to the nuclei boundary.
        max_threshold: value of the max threshold in µm for the distance to the nuclei boundary.
    Returns:
        gene2vect: Dictionary mapping gene names to nuclei distance vectors.
    """

    assert len(list_path_nuclei) == len(list_path_df_spots)
    assert len(list_path_nuclei) > 0, "No nuclei mask found"
    list_df_spots = []
    for index_crop in tqdm(range(len(list_path_nuclei))):
        nuclei_mask = tifffile.imread(list_path_nuclei[index_crop])
        df_spots = pd.read_csv(list_path_df_spots[index_crop])


        ## load json for scaling
        with open(Path(list_path_df_spots[index_crop]).parent / "bounds.json", "r") as f:
            dict_bounds = json.load(f)

        df_spots[x_column_name] = df_spots[x_column_name] - dict_bounds["bounds"][0]
        df_spots[y_column_name] = df_spots[y_column_name] - dict_bounds["bounds"][1]

        size_bounds = dict_bounds["bounds"][2] -  dict_bounds["bounds"][0]


        if int(size_bounds) != int(nuclei_mask.shape[0]): ## do not take non square input
            continue

        ### scaling ###
        df_spots[x_column_name] = df_spots[x_column_name] / dict_scale["x"]
        df_spots[y_column_name] = df_spots[y_column_name] / dict_scale["y"]
        #df_spots[z_column_name] = df_spots[z_column_name] / dict_scale["z"]

        assert df_spots[x_column_name].max() <= nuclei_mask.shape[0]
        assert df_spots[y_column_name].max() <= nuclei_mask.shape[0]

        df_spots = compute_rna2nuclei_distance(
            nuclei_mask= nuclei_mask,
            df_spots = df_spots,
            x_column_name = x_column_name,
            y_column_name = y_column_name,
            #z_column_name = z_column_name,
            use_max_threshold = use_max_threshold,
            max_threshold = max_threshold,
            pixel_um_size = pixel_um_size
        )
        list_df_spots.append(df_spots)

    assert len(list_df_spots)!=0, "problem with the nuclei mask, check that the nuclei mask is not resized"
    df_spots = pd.concat(list_df_spots)
    distance2nuclei_boundary_mean = df_spots.groupby(gene_column)["distance2nuclei_boundary"].mean()
    distance2nuclei_boundary_std = df_spots.groupby(gene_column)["distance2nuclei_boundary"].std()


    ### generate a color vector for each gene
    assert list(distance2nuclei_boundary_std.index)  == list(distance2nuclei_boundary_mean.index)

    if norm == "min_max":
        _min = np.array(distance2nuclei_boundary_mean).min()
        _max = np.array(distance2nuclei_boundary_mean).max()
        norm_distance2nuclei_boundary_mean = (np.array(distance2nuclei_boundary_mean) - _min) / (_max - _min)
    elif norm == "standard":
        _mean = np.array(distance2nuclei_boundary_mean).mean()
        _std = np.array(distance2nuclei_boundary_std).std()
        norm_distance2nuclei_boundary_mean = (np.array(distance2nuclei_boundary_mean) - _mean) / _std
    else:
        raise NotImplementedError
    gene2vect_mean = {}
    gene2vect_mean_sigma = {}
    gene2vect_unorm_mean_sigma = {}
    gene2vect_rgb = {}
    for gene in distance2nuclei_boundary_mean.index:
        gene_index = list(list(distance2nuclei_boundary_std.index)).index(gene)
        gene2vect_rgb[gene] = [1 - norm_distance2nuclei_boundary_mean[gene_index],
                               norm_distance2nuclei_boundary_mean[gene_index],
                               1 - norm_distance2nuclei_boundary_mean[gene_index]]

        gene2vect_mean[gene] = norm_distance2nuclei_boundary_mean[gene_index]
        gene2vect_mean_sigma[gene] = [norm_distance2nuclei_boundary_mean[gene_index], distance2nuclei_boundary_std[gene]]
        gene2vect_unorm_mean_sigma[gene] = [distance2nuclei_boundary_mean[gene_index], distance2nuclei_boundary_std[gene]]

    return gene2vect_mean, gene2vect_mean_sigma, gene2vect_unorm_mean_sigma, gene2vect_rgb, df_spots


#%%
if __name__ == "__main__":


    parser = argparse.ArgumentParser(description='RNA Spatial Embedding Generation')
    parser.add_argument('--path_to_save', type=str,
                        default="/home/tom/share/st3/open_vizgen/RNAseg_DATASET/mouse_ileum/RNAseg_DATASET.zarr/test_dict",
                        help='Path to save the results')
    parser.add_argument('--path_to_cache_path_folder',  nargs='+',
                        default=["/home/tom/share/st3/open_vizgen/RNAseg_DATASET/mouse_ileum/RNAseg_DATASET.zarr/.sopa_cache/rna_seg_1200"],
                        help='Path to the nuclei tif file')
    parser.add_argument('--experiement_name',
                        type=str,
                        default="0709_0309_cell_consistent_01_Cellbound1_dapi_cyto3_0.95_tii0.05_pt5_kernel10_set01_nuc_not_in_cell_nuc_bg",
                        help='Path to the spots CSV file')
    parser.add_argument('--dapi_seg', type=str, default="DAPI",
                        help='Path to the spots CSV file')
    parser.add_argument('--sample_crop', type=int, default=15,
                        help='')
    # parser gene_column
    parser.add_argument('--gene_column', type=str, default="gene",
                        help='Name of the column containing the gene name')

    parser.add_argument('--use_max_threshold', type=int, default=1,
                        help='Name of the column containing the gene name')

    parser.add_argument('--max_threshold', type=int, default=20,
                        help='Name of the column containing the gene name')


    parser.add_argument('--pixel_um_size', type=int, default=0.108,
                        help='Name of the column containing the gene name')

    parser.add_argument('--save_df', type=int, default=1,
                        help='Name of the column containing the gene name')

    parser.add_argument("--norm_gene_nuclei_dist", default='standard',
                        help="Normalization to apply to the distance vectors. choose in ['min_max', 'standard']")

    parser.add_argument("--port", default=3950)
    parser.add_argument("--mode", default='client')
    parser.add_argument("--host", default='127.0.0.2')
    args = parser.parse_args()
    assert isinstance(args.path_to_cache_path_folder, list), "path_to_cache_path_folder should be a list"




    ### set save model path
    args.e = datetime.now()
    args.date_str = f"{args.e.month}_d{args.e.day}_h{args.e.hour}_min{args.e.minute}_s{args.e.second}_r" + str(random.randint(0, 500))

    Path(args.path_to_save).mkdir(parents=True, exist_ok=True)
    np.save(Path(args.path_to_save) / f"script_parameter_{args.date_str}", args)
    with open(Path(args.path_to_save) / f"script_parameter_{args.date_str}.txt", "w") as f:
        for k, v in args.__dict__.items():
            f.write(f"{k} : {v}\n")
            print(f"{k} : {v}")


    list_path_df_spots = []
    list_path_nuclei = []
    list_path_index = []

    for path in args.path_to_cache_path_folder:
        print(path)
        list_path_index += list(Path(path).glob("*"))

    print('len(list_path_index)', len(list_path_index))


    # random shuffle the list of index
    np.random.shuffle(list_path_index)
    for path_index in list_path_index:
        path_df = path_index  / RNAsegFiles.TRANSCRIPTS_FILE
        df = pd.read_csv(path_df)
        if len(df) > 0:
            list_path_df_spots.append(path_df)
            assert (path_index / args.experiement_name / f"seg_{args.dapi_seg}.tif").is_file(), \
                f"{path_index / args.experiement_name / f'seg_{args.dapi_seg}'} is not a file cehck the segmentation path to nuclei segmentation"
            list_path_nuclei.append(str(path_index / args.experiement_name / f"seg_{args.dapi_seg}.tif"))
            print(f"Found {len(list_path_nuclei)}")

        if len(list_path_nuclei) == args.sample_crop:
            break




    gene2vect_dist_mean, gene2vect_dist_mean_sigma, gene2vect_dist_unorm_mean_sigma, gene2vect_dist_rgb, df_spots = get_gene_nuclei_dist_vector(
        list_path_nuclei= list_path_nuclei,
        list_path_df_spots=list_path_df_spots,
        norm = args.norm_gene_nuclei_dist,
        x_column_name="x",
        y_column_name="y",
        dict_scale={"x": 1, "y": 1, "z": 1},  # leave it to one if the df is already scaled in pixel
        gene_column = args.gene_column,
        use_max_threshold=args.use_max_threshold,
        max_threshold=args.max_threshold,
        pixel_um_size = args.pixel_um_size
    )


    np.save(Path(args.path_to_save) / f"gene2vect_dist_mean_{args.norm_gene_nuclei_dist}.npy", gene2vect_dist_mean)
    np.save(Path(args.path_to_save) / f"gene2vect_dist_mean_sigma_{args.norm_gene_nuclei_dist}.npy", gene2vect_dist_mean_sigma)
    np.save(Path(args.path_to_save) / f"gene2vect_dist_unorm_mean_sigma_{args.norm_gene_nuclei_dist}.npy", gene2vect_dist_unorm_mean_sigma)
    np.save(Path(args.path_to_save) / f"gene2vect_dist_rgb_{args.norm_gene_nuclei_dist}.npy", gene2vect_dist_rgb)

    print("saving df")

    if args.save_df:
        ## sample the df
        df_spots = df_spots.sample(min(2952117, len(df_spots)))

        df_spots.to_csv(Path(args.path_to_save) / f"df_spots_sampled_{args.norm_gene_nuclei_dist}.csv")
        print("done")







