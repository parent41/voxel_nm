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
    # ls ../WMH_micro_spatial/maps_UKB_space/sub-1*ses-2*${lists[i]}*.mnc > ./micro_file_lists/sub1_ses2_${lists[i]}.txt
    # ls ../WMH_micro_spatial/maps_UKB_space/sub-2*ses-2*${lists[i]}*.mnc > ./micro_file_lists/sub2_ses2_${lists[i]}.txt
    # ls ../WMH_micro_spatial/maps_UKB_space/sub-3*ses-2*${lists[i]}*.mnc > ./micro_file_lists/sub3_ses2_${lists[i]}.txt
    # ls ../WMH_micro_spatial/maps_UKB_space/sub-4*ses-2*${lists[i]}*.mnc > ./micro_file_lists/sub4_ses2_${lists[i]}.txt
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


#endregion

#region Tissue prevalence maps

# Make template mask dilated by 2 voxels (to make sure to get all brain voxels)
mincmorph -clobber -3D26 -successive DD ../UKB/temporary_template/Mask_2mm.mnc ../UKB/temporary_template/Mask_2mm_dil2.mnc

# Calculate tissue prevalence, mask of prevalence >0 and >=100

Rscript ./scripts/tissue_prevalence.R

#endregion

#region Matrices of micro values

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


#endregion















