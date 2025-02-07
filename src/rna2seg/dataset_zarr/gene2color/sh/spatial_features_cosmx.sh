#!/bin/bash

#SBATCH -N 1

#nombre threads max sur GPU 48acti

#SBATCH -n 1

#SBATCH --cpus-per-task=1

#SBATCH -J 500spatial_features

#SBATCH --output="500spatial_features.out"


#SBATCH --partition=cbio-cpu

#SBATCH --mem-per-cpu=20000




#randint=$RANDOM

module load cuda/11.3

cd ..

echo "start from slurmscript"

path_to_cache_path_folder=$1

echo "path_to_cache_path_folder: $path_to_cache_path_folder"

python spatial_features.py \
--path_to_save "/cluster/CBIO/data1/st_segmentation/cosmx/dict_color" \
--path_to_cache_path_folder $path_to_cache_path_folder \
--experiement_name "1209_0309_cell_consistent_01_CD3_d60_0.95_tii0.05_pt5_kernel10_" \
--norm_gene_nuclei_dist standard \
--dapi_seg "DAPI_d70" \
--sample_crop 2000 \
--gene_column 'target' \

echo "done from slurmscript"