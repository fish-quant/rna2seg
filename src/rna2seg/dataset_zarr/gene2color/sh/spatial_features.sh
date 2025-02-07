#!/bin/bash

#SBATCH -N 1

#nombre threads max sur GPU 48acti

#SBATCH -n 1

#SBATCH --cpus-per-task=1

#SBATCH -J 500spatial_features

#SBATCH --output="500spatial_features.out"


#SBATCH --partition=cbio-cpu

#SBATCH --mem-per-cpu=40000




#randint=$RANDOM

module load cuda/11.3

cd ..

organ=$1


# Ovarian  Prostate,  Uterine



if [ $organ == "Breast" ];
then
echo "start from slurmscript"
echo "organ: $organ"
python spatial_features.py \
--path_to_save "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/breast/0710_dict_color_sample500" \
--path_to_cache_path_folder "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/breast/entire_dataset_z3_cb1_dapi.zarr/.sopa_cache/rna_seg_1200" \
--experiement_name "0709_0309_cell_consistent_01_Cellbound1_dapi_cyto3_d60_2008_tic0.95_tii0.05_pt5_kernel10_set01_nuc_not_in_cell_nuc_bg" \
--norm_gene_nuclei_dist standard \
--dapi_seg "DAPI" \
--sample_crop 1500 \

fi




if [ $organ == "Lung" ];
then
echo "start from slurmscript"
echo "organ: $organ"
python spatial_features.py \
--path_to_save "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/lung/0710_dict_color_sample500" \
--path_to_cache_path_folder "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/lung/lung_patient1.zarr/.sopa_cache/rna_seg_1200" \
--experiement_name "0709_0309_cell_consistent_01_Cellbound1_dapi_cyto3_d60_2008_0.95_tii0.05_pt5_kernel10_set01_nuc_not_in_cell_nuc_bg" \
--norm_gene_nuclei_dist standard \
--dapi_seg "DAPI" \
--sample_crop 1500 \

fi




if [ $organ == "Melanoma" ];
then
echo "start from slurmscript"
echo "organ: $organ"
python spatial_features.py \
--path_to_save "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/Melanoma/0710_dict_color_sample500" \
--path_to_cache_path_folder "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/Melanoma/HumanMelanomaPatient1.zarr/.sopa_cache/rna_seg_1200" \
--experiement_name "0709_0309_cell_consistent_01_Cellbound1_dapi_cyto3_d60_2008__0.95_tii0.05_pt5_kernel10_set01_nuc_not_in_cell_nuc_bg" \
--norm_gene_nuclei_dist standard \
--dapi_seg "DAPI_d90" \
--sample_crop 1500 \

fi


if [ $organ == "Colon" ];
then
echo "start from slurmscript"
echo "organ: $organ"
python spatial_features.py \
--path_to_save "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/colon/0710_dict_color_sample500" \
--path_to_cache_path_folder "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/colon/HumanColonCancerPatient2.zarr/.sopa_cache/rna_seg_1200" \
--experiement_name "0709_0309_cell_consistent_01_Cellbound1_dapi_cyto3_d60_2008_0.95_tii0.05_pt5_kernel10_set01_nuc_not_in_cell_nuc_bg" \
--norm_gene_nuclei_dist standard \
--dapi_seg "DAPI" \
--sample_crop 1500 \

fi



if [ $organ == "Ovarian" ];
then
echo "start from slurmscript"
echo "organ: $organ"
python spatial_features.py \
--path_to_save "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/Ovarian/0710_dict_color_sample500" \
--path_to_cache_path_folder "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/Ovarian/HumanOvarianCancerPatient2Slice1.zarr/.sopa_cache/rna_seg_1200" \
--experiement_name "0709_0309_cell_consistent_01_Cellbound1_dapi_cyto3_d60_2008_0.95_tii0.05_pt5_kernel10_set01_nuc_not_in_cell_nuc_bg" \
--norm_gene_nuclei_dist standard \
--dapi_seg "DAPI_d90" \
--sample_crop 1500 \

fi


if  [ $organ == "Prostate" ];
then
echo "start from slurmscript"
echo "organ: $organ"
python spatial_features.py \
--path_to_save "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/Prostate/0710_dict_color_sample500" \
--path_to_cache_path_folder "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/Prostate/HumanProstateCancerPatient1.zarr/.sopa_cache/rna_seg_1200" \
--experiement_name "0709_0309_cell_consistent_01_Cellbound3_dapi_cyto3_d60_2008_0.95_tii0.05_pt5_kernel10_set01_nuc_not_in_cell_nuc_bg" \
--norm_gene_nuclei_dist standard \
--dapi_seg "DAPI_d90" \
--sample_crop 1500 \

fi

if  [ $organ == "Uterine" ];
then
echo "start from slurmscript"
echo "organ: $organ"
python spatial_features.py \
--path_to_save "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/Uterine/0710_dict_color_sample500" \
--path_to_cache_path_folder "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/Uterine/HumanUterineCancerPatient1.zarr/.sopa_cache/rna_seg_1200" \
--experiement_name "0709_0309_cell_consistent_01_Cellbound1_dapi_cyto3_d60_2008_0.95_tii0.05_pt5_kernel10_set01_nuc_not_in_cell_nuc_bg" \
--norm_gene_nuclei_dist standard \
--dapi_seg "DAPI_d90" \
--sample_crop 1500 \

fi

if  [ $organ == "mouse_ileum" ];
then
echo "start from slurmscript"
echo "organ: $organ"
python spatial_features.py \
--path_to_save "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/mouse_ileum/0710_dict_color_sample500" \
--path_to_cache_path_folder "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/mouse_ileum/RNAseg_DATASET.zarr/.sopa_cache/rna_seg_1200" \
--experiement_name "0709_0309_cell_consistent_01_Cellbound1_dapi_cyto3_0.95_tii0.05_pt5_kernel10_set01_nuc_not_in_cell_nuc_bg" \
--norm_gene_nuclei_dist standard \
--dapi_seg "DAPI" \
--sample_crop 1000 \

fi




echo  "done from slurmscript"