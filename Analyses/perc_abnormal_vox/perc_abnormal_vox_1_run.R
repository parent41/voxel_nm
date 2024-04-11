
library(RMINC)
library(data.table)
library(matrixStats)

# Arg 1: Micro name
args = commandArgs(trailingOnly=TRUE)

# args = c()
# args[1] = "FA"
# args[2] = "../subj_zscores_common_all/tmp/label_c2_FA.tsv"
# args[3] = "../subj_zscores_common_all/results/zscores_c2_FA_anlm.tsv"
# args[4] = "../../../WMH_micro_spatial/QC/inclusions_without_excluding_dx_new.txt"
# args[5] = "../subj_zscores_common_all/tmp/ids_label_c2_FA.txt"
# args[6] = "./results/test.tsv"

tissue_nm = c('Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM')
tissue_all=c('Ventricules', 'CSF', 'Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH')

ids = as.data.frame(fread(args[5]))
ids = ids$V1

# Function to count values above zscore thresholds per label
perc_values_above_thresholds <- function(id, micro_id, label_id) {
  thresholds = c(1,2,3)

  result = expand.grid(label_value = 3:9, threshold = thresholds)
  result = cbind(rep(id, nrow(result)), result)
  colnames(result)[1] = "ID"

  result$count_vox_label = NA
  result$count_vox_above_thresh = NA
  result$perc_vox_above_thresh = NA

  for (i in 1:nrow(result)) {
    result$count_vox_label[i] = length(micro_id[which(label_id == result$label_value[i])])
    result$count_vox_above_thresh[i] = sum(abs(micro_id[which(label_id == result$label_value[i])]) > result$threshold[i], na.rm=TRUE)
    result$perc_vox_above_thresh[i] = result$count_vox_above_thresh[i] / result$count_vox_label[i]
  }

  return(result)
}

# Calculate % abnormal voxels per tissue type
print("Load voxel data")

label = as.data.frame(fread(args[2]))
micro = as.data.frame(fread(args[3], header=TRUE))

print(dim(label))
print(dim(micro))
print(length(ids))

print("Run for each ID")


results_list = list()
for (id in 1:nrow(micro)) {
  print(ids[id])
  results_list[[id]] = perc_values_above_thresholds(ids[id], as.numeric(micro[id,]), as.numeric(label[id,]))
}

results = do.call(rbind, results_list)

fwrite(as.data.frame(results), args[6], row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
