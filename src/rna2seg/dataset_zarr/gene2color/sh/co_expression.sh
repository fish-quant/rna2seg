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

cd /cluster/CBIO/data1/data3/tdefard/st_seg_code/st_seg/CNN_gene/rna_seg/dataset_zarr/gene2color || exit

organ=$1
# Ovarian  Prostate,  Uterine



if [ $organ == "Breast" ];
then
echo "start from slurmscript"
echo "organ: $organ"
python coexpression.py \


fi




if [ $organ == "Lung" ];
then
echo "start from slurmscript"
echo "organ: $organ"
python coexpression.py \

fi




if [ $organ == "Melanoma" ];
then
echo "start from slurmscript"
echo "organ: $organ"
python coexpression.py \


fi


if [ $organ == "Colon" ];
then
echo "start from slurmscript"
echo "organ: $organ"
python coexpression.py \

fi



if [ $organ == "Ovarian" ];
then
echo "start from slurmscript"
echo "organ: $organ"
python coexpression.py \
--path_to_save "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/all_dictionary/ovarian/" \
--path_to_cache "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/Ovarian/HumanOvarianCancerPatient2Slice1.zarr/.sopa_cache/rna_seg_1200" \
--sampling_size 500000 \
--nb_crop_sample 250 \

fi


if  [ $organ == "Prostate" ];
then
echo "start from slurmscript"
echo "organ: $organ"
python coexpression.py \
--path_to_save "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/all_dictionary/prostate/" \
--path_to_cache "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/Prostate/HumanProstateCancerPatient1.zarr/.sopa_cache/rna_seg_1200" \
--sampling_size 500000 \
--nb_crop_sample 250 \


fi

if  [ $organ == "Uterine" ];
then
echo "start from slurmscript"
echo "organ: $organ"
python coexpression.py \
--path_to_save "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/all_dictionary/uterine/" \
--path_to_cache "/cluster/CBIO/data1/st_segmentation/open_vizgen/RNAseg_DATASET/Uterine/HumanUterineCancerPatient1.zarr/.sopa_cache/rna_seg_1200" \
--sampling_size 500000 \
--nb_crop_sample 250 \

fi

if  [ $organ == "mouse_ileum" ];
then
echo "start from slurmscript"
echo "organ: $organ"
python coexpression.py \


fi




echo  "done from slurmscript"