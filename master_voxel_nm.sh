# Voxel-wise, tissue-specific normative modelling of microstructure in UK Biobank

#################################################################################
############################ Make micro matrices ################################
#################################################################################

#region Matrices of BISON labels

# Make BISON file lists in chunks
lists=('Label')

for i in ${!lists[@]}
do
    echo ${lists[i]}
    ls ../WMH_micro_spatial/maps_UKB_space/sub-1*ses-2*${lists[i]}*.mnc > ./micro_file_lists/sub1_ses2_${lists[i]}.txt
    ls ../WMH_micro_spatial/maps_UKB_space/sub-2*ses-2*${lists[i]}*.mnc > ./micro_file_lists/sub2_ses2_${lists[i]}.txt
    ls ../WMH_micro_spatial/maps_UKB_space/sub-3*ses-2*${lists[i]}*.mnc > ./micro_file_lists/sub3_ses2_${lists[i]}.txt
    ls ../WMH_micro_spatial/maps_UKB_space/sub-4*ses-2*${lists[i]}*.mnc > ./micro_file_lists/sub4_ses2_${lists[i]}.txt
    ls ../WMH_micro_spatial/maps_UKB_space/sub-5*ses-2*${lists[i]}*.mnc > ./micro_file_lists/sub5_ses2_${lists[i]}.txt
    ls ../WMH_micro_spatial/maps_UKB_space/sub-6*ses-2*${lists[i]}*.mnc > ./micro_file_lists/sub6_ses2_${lists[i]}.txt

    ls ../WMH_micro_spatial/maps_UKB_space/*ses-3*${lists[i]}*.mnc > ./micro_file_lists/ses3_${lists[i]}.txt
done

# Make BISON matrices

for file in ./micro_file_lists/*Label*
do
    echo Rscript ./scripts/bison_matrices.R $file ../UKB/temporary_template/Mask_2mm_dil2.mnc ./bison_matrices/$(basename $file .txt)_whole_brain.tsv
done > joblist_bison_matrices

qbatch -c 1 -w 1:00:00 joblist_bison_matrices

# Concatenate together for ses2 (and delete chunks)
cat ./bison_matrices/sub*_ses2_Label_whole_brain.tsv > ./bison_matrices/ses2_Label_whole_brain.tsv
cat ./bison_matrices/ids_sub*_ses2_Label_whole_brain.txt > ./bison_matrices/ids_ses2_Label_whole_brain.txt

rm ./bison_matrices/sub*_ses2_Label_whole_brain.tsv
rm ./bison_matrices/ids_sub*_ses2_Label_whole_brain.txt

#endregion

#region Tissue prevalence maps

# Make template mask dilated by 2 voxels (to make sure to get all brain voxels)
mincmorph -clobber -3D26 -successive DD ../UKB/temporary_template/Mask_2mm.mnc ../UKB/temporary_template/Mask_2mm_dil2.mnc

# Calculate tissue prevalence, mask of prevalence >0 and >=100

Rscript ./scripts/tissue_prevalence.R

#endregion

#region Tissue-specific matrices of microstructure

# Create file lists

micro=('dti_FA' 'dti_MD' 'NODDI_ICVF' 'NODDI_ECVF' 'NODDI_ISOVF' 'NODDI_OD' 'T2star' 'QSM')

for i in ${!micro[@]}
do
    echo ${micro[i]}
    ls ../WMH_micro_spatial/maps_UKB_space/sub-1*ses-2*${micro[i]}*.mnc > ./micro_file_micro/sub1_ses2_${micro[i]}.txt
    ls ../WMH_micro_spatial/maps_UKB_space/sub-2*ses-2*${micro[i]}*.mnc > ./micro_file_lists/sub2_ses2_${micro[i]}.txt
    ls ../WMH_micro_spatial/maps_UKB_space/sub-3*ses-2*${micro[i]}*.mnc > ./micro_file_lists/sub3_ses2_${micro[i]}.txt
    ls ../WMH_micro_spatial/maps_UKB_space/sub-4*ses-2*${micro[i]}*.mnc > ./micro_file_lists/sub4_ses2_${micro[i]}.txt
    ls ../WMH_micro_spatial/maps_UKB_space/sub-5*ses-2*${micro[i]}*.mnc > ./micro_file_lists/sub5_ses2_${micro[i]}.txt
    ls ../WMH_micro_spatial/maps_UKB_space/sub-6*ses-2*${micro[i]}*.mnc > ./micro_file_lists/sub6_ses2_${micro[i]}.txt

    ls ../WMH_micro_spatial/maps_UKB_space/*ses-3*${micro[i]}*.mnc > ./micro_file_lists/ses3_${micro[i]}.txt
done

# Make matrices

micro=('dti_FA' 'dti_MD' 'NODDI_ICVF' 'NODDI_ECVF' 'NODDI_ISOVF' 'NODDI_OD' 'T2star' 'QSM')
new_micro=('FA' 'MD' 'ICVF' 'ECVF' 'ISOVF' 'OD' 'T2star' 'QSM')

for i in ${!micro[@]}
do
    echo Rscript ./scripts/micro_matrices.R ./micro_file_lists/sub1_ses2_${micro[i]}.txt ../UKB/temporary_template/Mask_2mm_dil2.mnc sub1_ses2_${new_micro[i]}
    echo Rscript ./scripts/micro_matrices.R ./micro_file_lists/sub2_ses2_${micro[i]}.txt ../UKB/temporary_template/Mask_2mm_dil2.mnc sub2_ses2_${new_micro[i]}
    echo Rscript ./scripts/micro_matrices.R ./micro_file_lists/sub3_ses2_${micro[i]}.txt ../UKB/temporary_template/Mask_2mm_dil2.mnc sub3_ses2_${new_micro[i]}
    echo Rscript ./scripts/micro_matrices.R ./micro_file_lists/sub4_ses2_${micro[i]}.txt ../UKB/temporary_template/Mask_2mm_dil2.mnc sub4_ses2_${new_micro[i]}
    echo Rscript ./scripts/micro_matrices.R ./micro_file_lists/sub5_ses2_${micro[i]}.txt ../UKB/temporary_template/Mask_2mm_dil2.mnc sub5_ses2_${new_micro[i]}
    echo Rscript ./scripts/micro_matrices.R ./micro_file_lists/sub6_ses2_${micro[i]}.txt ../UKB/temporary_template/Mask_2mm_dil2.mnc sub6_ses2_${new_micro[i]}

    echo Rscript ./scripts/micro_matrices.R ./micro_file_lists/ses3_${micro[i]}.txt ../UKB/temporary_template/Mask_2mm_dil2.mnc ses3_${new_micro[i]}
done > joblist_make_micro_matrices

qbatch -c 1 -w 2:00:00 joblist_make_micro_matrices

# Concatenate together for ses2 (and delete chunks)
for i in ${!new_micro[@]}
do
    echo ${new_micro[i]}
    cat ./micro_matrices/sub*_ses2_${new_micro[i]}.tsv > ./micro_matrices/ses2_${new_micro[i]}.tsv
    cat ./micro_matrices/ids_sub*_ses2_${new_micro[i]}.txt > ./micro_matrices/ids_ses2_${new_micro[i]}.txt
    rm ./micro_matrices/sub*_ses2_${new_micro[i]}.tsv
    rm ./micro_matrices/ids_sub*_ses2_${new_micro[i]}.txt
done

#endregion

#################################################################################
############################## Normative models #################################
#################################################################################

#region Bayesian Linear Regression with age and sex

module load cobralab
module load python/3.9.8

cd ~/.virtualenvs
virtualenv --system-site-packages ~/.virtualenvs/PCN_env
source ~/.virtualenvs/PCN_env/bin/activate

git clone https://github.com/predictive-clinical-neuroscience/PCNtoolkit-demo.git

pip install -r PCNtoolkit-demo/tutorials/BLR_protocol/requirements.txt
pip install --upgrade numpy
pip install --upgrade pytensor
pip install --upgrade numba
pip install plotnine
pip install pyarrow
pip install dask
pip install --upgrade pandas "dask[complete]"
pip install polars

# Run
module load cobralab
module load python/3.9.8
source ~/.virtualenvs/PCN_env/bin/activate

for micro in $(seq 0 6)
do
    for chunk in $(seq 0 30) # 31 chunks of 10,000 voxels
    do
        echo python ./norm_models_BLR_1_run.py ${micro} ${chunk} 
    done
done > joblist_norm_models_BLR

qbatch -c 1 -w 0:30:00 joblist_test
qbatch -c 1 -w 0:30:00 joblist_norm_models_BLR

#endregion

