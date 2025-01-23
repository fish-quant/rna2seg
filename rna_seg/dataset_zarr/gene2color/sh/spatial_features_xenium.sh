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

path_to_cache_path_folder=/cluster/CBIO/data1/st_segmentation/xenium_dataset/xenium_dataset/Xenium_Prime_Human_Skin_FFPE_outs.zarr/.sopa_cache/rna_seg_1200

echo "path_to_cache_path_folder: $path_to_cache_path_folder"

python spatial_features.py \
--path_to_save "/cluster/CBIO/data1/st_segmentation/xenium_dataset/xenium_dataset/dict_color" \
--path_to_cache_path_folder $path_to_cache_path_folder \
--experiement_name "2409_cell_consistent_cell_boundaries10_p5_2008" \
--norm_gene_nuclei_dist standard \
--dapi_seg "nucleus_boundaries" \
--sample_crop 2000 \
--gene_column 'feature_name' \

echo "done from slurmscript"