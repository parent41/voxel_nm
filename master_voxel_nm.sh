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
    echo Rscript ./scripts/make_label_matrix.R $file ../UKB/temporary_template/Mask_2mm.mnc ./micro_matrices/$(basename $file .txt)_whole_brain.tsv
done > joblist_bison_matrices

qbatch -c 1 -w 2:00:00 joblist_bison_matrices


#endregion

#region Tissue prevalence maps

# Make template mask dilated by 2 voxels (to make sure to get all brain voxels)
mincmorph -clobber -3D26 -successive DD ../UKB/temporary_template/Mask_2mm.mnc ../UKB/temporary_template/Mask_2mm_dil2.mnc



#endregion

















