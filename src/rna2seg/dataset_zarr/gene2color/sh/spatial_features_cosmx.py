








import os
import pandas as pd
import numpy as np
from pathlib import Path
import argparse
from tqdm import tqdm
import tifffile




if __name__ == "__main__":

    path_out = "/cluster/CBIO/data1/data3/tdefard/st_seg_code/st_seg/CNN_gene/rna_seg/dataset_zarr/sh/cosmx/out"
    path_sh = "/cluster/CBIO/data1/data3/tdefard/st_seg_code/st_seg/CNN_gene/rna_seg/dataset_zarr/sh/cosmx/dataset_preprocessing.sh"
    path_lung_folder = "/cluster/CBIO/data1/st_segmentation/cosmx"
    list_valid_zarr = []
    list_valid_zarr_folder = []
    folder_to_use = ["Lung5_Rep1", "Lung5_Rep2", "Lung5_Rep3", "Lung6"]
    for path_lung in Path(path_lung_folder).glob("Lung*"):

        if path_lung.stem not in folder_to_use:
            continue

        path_folder_zarr = path_lung / "fov_zarr_0809"

        for path_zarr in path_folder_zarr.glob("*.zarr"):
            path_zarr_shape = path_zarr / "shapes"
            print(path_zarr_shape)
            present_shape = [path_shape.name for path_shape in path_zarr_shape.glob("*")]

            required_shape = [
                "DAPI_d70",
                "Membrane_d60",
                "CD45_d60",
                "CD3_d60",
                "PanCK_d60",
            ]

            if  len(set(required_shape) - set(present_shape))==0:
                list_valid_zarr.append(str(path_zarr) + "/.sopa_cache/rna_seg_1200")

        path_to_cache_path_folder = ""
        for pt in list_valid_zarr:
            path_to_cache_path_folder += f'{pt} '



    print(f'runing {path_zarr.stem}')
    print(f'runing {len(path_to_cache_path_folder)} zarr')
    # SBATCH --nodelist=node008
    os.system(f'sbatch --output spatial_cosmx.out '
              f'-J spatial_cosmx spatial_features_cosmx.sh "{path_to_cache_path_folder[:-1]}"')





















