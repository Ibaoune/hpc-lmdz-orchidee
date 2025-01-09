#!/bin/bash
#Author: M. El Aabaribaoune (@um6p)

path_to_lmdz_simus="/home/${USER}/lustre/climat-um6p-st-iwri-7ksifkvwkuy/users/${USER}/lmdz_simus/"
if [[ ! -e $path_to_lmdz_simus ]]; then
    mkdir $path_to_lmdz_simus
fi

path_to_similar_dirs_jz="${path_to_lmdz_simus}/similar_dirs_jz_tbkl"
if [[ ! -e $path_to_similar_dirs_jz ]]; then
    mkdir $path_to_similar_dirs_jz
fi

working_path=$(pwd)

cd $path_to_similar_dirs_jz
mkdir STORE_TBKL STORAGE_TBKL SCRATCH_TBKL ARCHIVE_TBKL
mkdir gpfswork_tbkl gpfsstore_tbkl gpfsscratch_tbkl 
mkdir gpfsscratch_tbkl/rech  gpfswork_tbkl/rech

# Go back to the model path
cd $working_path
