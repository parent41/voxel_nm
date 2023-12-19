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

#region Denoised micro maps in UKB space (dx only)

mkdir maps_UKB_space_anlm_dx

# Denoise maps

while IFS= read -r id; do
    echo ./scripts/denoise_micro.sh sub-${id}_ses-2 ../WMH_micro_spatial/maps_UKB_space ./maps_UKB_space_anlm_dx
done < "../WMH_micro_spatial/QC/inclusions_only_dx_new.txt" > joblist_anlm_micro_maps_dx

qbatch -c 200 joblist_anlm_micro_maps_dx

# Make matrices

micro_long=('dti_FA' 'dti_MD' 'NODDI_ICVF' 'NODDI_ISOVF' 'NODDI_OD' 'T2star' 'QSM')
micro=('FA' 'MD' 'ICVF' 'ISOVF' 'OD' 'T2star' 'QSM')

for i in ${!micro[@]}
do
    echo Rscript ./scripts/micro_matrices_anlm_dx.R ${micro_long[i]} ${micro[i]}
done > joblist_micro_matrices_anlm_dx

qbatch -c 1 joblist_micro_matrices_anlm_dx

#endregion

#region Denoised micro maps in UKB space (all subjects)

# Denoise micro maps

for file in ../WMH_micro_spatial/maps_UKB_space/*_dti_FA_UKB.mnc
do
    echo ./scripts/denoise_micro.sh $(basename $file _dti_FA_UKB.mnc) ../WMH_micro_spatial/maps_UKB_space ./maps_UKB_space_anlm_all
done > joblist_anlm_micro_maps_all

head joblist_anlm_micro_maps_all -n 10000 > joblist_anlm_micro_maps_all_10000
head joblist_anlm_micro_maps_all -n 20000 | tail -n 10000 > joblist_anlm_micro_maps_all_20000
head joblist_anlm_micro_maps_all -n 30000 | tail -n 10000 > joblist_anlm_micro_maps_all_30000
head joblist_anlm_micro_maps_all -n 40000 | tail -n 10000 > joblist_anlm_micro_maps_all_40000
tail joblist_anlm_micro_maps_all -n 1889 > joblist_anlm_micro_maps_all_41889

qbatch -c 500 joblist_anlm_micro_maps_all_10000
qbatch -c 500 joblist_anlm_micro_maps_all_20000
qbatch -c 500 joblist_anlm_micro_maps_all_30000
qbatch -c 500 joblist_anlm_micro_maps_all_40000
qbatch -c 500 joblist_anlm_micro_maps_all_41889

# Create file lists

micro=('dti_FA' 'dti_MD' 'NODDI_ICVF' 'NODDI_ISOVF' 'NODDI_OD' 'T2star' 'QSM')

for i in ${!micro[@]}
do
    echo ${micro[i]}
    ls ./maps_UKB_space_anlm_all/sub-1*ses-2*${micro[i]}*_UKB_anlm.mnc > ./micro_file_lists/sub1_ses2_${micro[i]}_anlm.txt
    ls ./maps_UKB_space_anlm_all/sub-2*ses-2*${micro[i]}*_UKB_anlm.mnc > ./micro_file_lists/sub2_ses2_${micro[i]}_anlm.txt
    ls ./maps_UKB_space_anlm_all/sub-3*ses-2*${micro[i]}*_UKB_anlm.mnc > ./micro_file_lists/sub3_ses2_${micro[i]}_anlm.txt
    ls ./maps_UKB_space_anlm_all/sub-4*ses-2*${micro[i]}*_UKB_anlm.mnc > ./micro_file_lists/sub4_ses2_${micro[i]}_anlm.txt
    ls ./maps_UKB_space_anlm_all/sub-5*ses-2*${micro[i]}*_UKB_anlm.mnc > ./micro_file_lists/sub5_ses2_${micro[i]}_anlm.txt
    ls ./maps_UKB_space_anlm_all/sub-6*ses-2*${micro[i]}*_UKB_anlm.mnc > ./micro_file_lists/sub6_ses2_${micro[i]}_anlm.txt

    ls ./maps_UKB_space_anlm_all/*ses-3*${micro[i]}*_UKB_anlm.mnc > ./micro_file_lists/ses3_${micro[i]}_anlm.txt
done

# Make matrices

micro=('dti_FA' 'dti_MD' 'NODDI_ICVF' 'NODDI_ISOVF' 'NODDI_OD' 'T2star' 'QSM')
new_micro=('FA' 'MD' 'ICVF' 'ISOVF' 'OD' 'T2star' 'QSM')

for i in ${!micro[@]}
do
    echo Rscript ./scripts/micro_matrices.R ./micro_file_lists/sub1_ses2_${micro[i]}_anlm.txt ../UKB/temporary_template/Mask_2mm_dil2.mnc sub1_ses2_${new_micro[i]}_anlm
    echo Rscript ./scripts/micro_matrices.R ./micro_file_lists/sub2_ses2_${micro[i]}_anlm.txt ../UKB/temporary_template/Mask_2mm_dil2.mnc sub2_ses2_${new_micro[i]}_anlm
    echo Rscript ./scripts/micro_matrices.R ./micro_file_lists/sub3_ses2_${micro[i]}_anlm.txt ../UKB/temporary_template/Mask_2mm_dil2.mnc sub3_ses2_${new_micro[i]}_anlm
    echo Rscript ./scripts/micro_matrices.R ./micro_file_lists/sub4_ses2_${micro[i]}_anlm.txt ../UKB/temporary_template/Mask_2mm_dil2.mnc sub4_ses2_${new_micro[i]}_anlm
    echo Rscript ./scripts/micro_matrices.R ./micro_file_lists/sub5_ses2_${micro[i]}_anlm.txt ../UKB/temporary_template/Mask_2mm_dil2.mnc sub5_ses2_${new_micro[i]}_anlm
    echo Rscript ./scripts/micro_matrices.R ./micro_file_lists/sub6_ses2_${micro[i]}_anlm.txt ../UKB/temporary_template/Mask_2mm_dil2.mnc sub6_ses2_${new_micro[i]}_anlm

    echo Rscript ./scripts/micro_matrices.R ./micro_file_lists/ses3_${micro[i]}_anlm.txt ../UKB/temporary_template/Mask_2mm_dil2.mnc ses3_${new_micro[i]}_anlm
done > joblist_make_micro_matrices_anlm

qbatch -c 1 -w 2:00:00 joblist_make_micro_matrices_anlm

# Concatenate together for ses2 (and delete chunks)
for i in ${!new_micro[@]}
do
    echo ${new_micro[i]}
    # cat ./micro_matrices/sub*_ses2_${new_micro[i]}_anlm.tsv > ./micro_matrices/ses2_${new_micro[i]}_anlm.tsv
    # cat ./micro_matrices/ids_sub*_ses2_${new_micro[i]}_anlm.txt > ./micro_matrices/ids_ses2_${new_micro[i]}_anlm.txt
    rm ./micro_matrices/sub*_ses2_${new_micro[i]}_anlm.tsv
    rm ./micro_matrices/ids_sub*_ses2_${new_micro[i]}_anlm.txt
done


#endregion

#################################################################################
################################## Analyses #####################################
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

#region Voxel-wises z-score for subjects with dx in common space

cd Analyses/subj_zscores_common_dx

# Make tsv of micro and label for dx only

module load cobralab

micro=('FA' 'MD' 'ICVF' 'ISOVF' 'OD' 'T2star' 'QSM')

for i in ${!micro[@]}
do
    echo Rscript subj_zscores_common_dx_1_select_rows.R ${micro[i]}
done > joblist_subj_zscores_common_dx_1_select_rows

qbatch -c 1 -w 1:00:00 joblist_subj_zscores_common_dx_1_select_rows

# Z-scores for participants with dx

module load cobralab
module load python/3.9.8
source ~/.virtualenvs/PCN_env/bin/activate

for micro in $(seq 0 6)
do
    echo python ./subj_zscores_common_dx_2_zscore.py ${micro}
done > joblist_subj_zscores_common_dx_2_zscore

qbatch -c 1 -w 1:30:00 joblist_subj_zscores_common_dx_2_zscore

# Make mnc files

micro=('FA' 'MD' 'ICVF' 'ISOVF' 'OD' 'T2star' 'QSM')

for i in ${!micro[@]}
do
    echo Rscript subj_zscores_common_dx_3_mnc.R ${micro[i]}
done > joblist_subj_zscores_common_dx_3_mnc

qbatch -c 1 -w 1:30:00 joblist_subj_zscores_common_dx_3_mnc

# Visualize

echo Rscript subj_zscores_common_dx_4_viz.R > joblist_subj_zscores_common_dx_4_viz

qbatch -c 1 -w 3:00:00 joblist_subj_zscores_common_dx_4_viz

# Add histograms

echo Rscript ./subj_zscores_common_dx_5_viz_hist.R > joblist_subj_zscores_common_dx_5_viz_hist

qbatch -c 1 -w 12:00:00 joblist_subj_zscores_common_dx_5_viz_hist

#endregion

#region Voxel-wises z-score for healthy subjects in common space

# Z-scores

module load cobralab
module load python/3.9.8
source ~/.virtualenvs/PCN_env/bin/activate

micro=('FA' 'MD' 'ICVF' 'ISOVF' 'OD' 'T2star' 'QSM')

for m in ${!micro[@]}
do
    for i in $(seq 0 32)
    do
        echo python ./subj_zscores_common_hc_1_zscore.py ${m} ${i}
    done
done > joblist_subj_zscores_common_hc_1_zscore

qbatch -c 2 joblist_subj_zscores_common_hc_1_zscore
qbatch -c 2 -w 3:20:00 joblist_subj_zscores_common_hc_1_zscore

# Make mnc files

micro=('FA' 'MD' 'ICVF' 'ISOVF' 'OD' 'T2star' 'QSM')

for i in ${!micro[@]}
do
    echo Rscript subj_zscores_common_hc_2_mnc.R ${micro[i]}
done > joblist_subj_zscores_common_hc_2_mnc

qbatch -c 1 -w 3:00:00 joblist_subj_zscores_common_hc_2_mnc

# Visualize

echo Rscript ./subj_zscores_common_hc_3_viz.R > joblist_subj_zscores_common_hc_3_viz

qbatch -c 1 -w 3:00:00 joblist_subj_zscores_common_hc_3_viz

# Add histograms

echo Rscript ./subj_zscores_common_hc_4_viz_hist.R > joblist_subj_zscores_common_hc_4_viz_hist

qbatch -c 1 -w 12:00:00 joblist_subj_zscores_common_hc_4_viz_hist

#endregion

#region Percentage of abnormal voxels per tissue type

micro=('FA' 'MD' 'ICVF' 'ISOVF' 'OD' 'T2star' 'QSM')

# args[1] = "FA"
# args[2] = "../subj_zscores_common_dx/results/ses2_Label_whole_brain_dx.tsv"
# args[3] = "../subj_zscores_common_dx/results/zscores_anlm_FA.tsv"
# args[4] = "../../../WMH_micro_spatial/QC/inclusions_only_dx_new.txt"
# args[5] = "./results/perc_abnormal_dx_FA_anlm.tsv"

# first for dx
for m in ${!micro[@]}
do
    echo Rscript ./perc_abnormal_vox.R ${micro[m]} \
        ../subj_zscores_common_dx/results/ses2_Label_whole_brain_dx.tsv \
        ../subj_zscores_common_dx/results/zscores_anlm_${micro[m]}.tsv \
        ../../../WMH_micro_spatial/QC/inclusions_only_dx_new.txt \
        ./results/raw/perc_abnormal_dx_${micro[m]}_anlm.tsv
done > joblist_perc_abnormal_vox

# then for hc

for m in ${!micro[@]}
do
    for i in $(seq 0 32)
    do
        echo Rscript ./perc_abnormal_vox.R ${micro[m]} \
                ../subj_zscores_common_hc/tmp/label_c${i}_FA.tsv \
                ../subj_zscores_common_hc/results/zscores_c${i}_${micro[m]}_anlm.tsv \
                ../subj_zscores_common_hc/tmp/ids_label_c${i}_${micro[m]}.txt \
                ./results/perc_abnormal_hc_c${i}_${micro[m]}_anlm.tsv
    done >> joblist_perc_abnormal_vox
done

qbatch -c 5 -w 2:45:00 joblist_perc_abnormal_vox



#endregion

#region Z-score averages by dx

for i in $(seq 1 21)
do
    echo Rscript dx_average.R ${i}
done > joblist_dx_average

qbatch -c 1 -w 2:00:00 joblist_dx_average


#endregion
