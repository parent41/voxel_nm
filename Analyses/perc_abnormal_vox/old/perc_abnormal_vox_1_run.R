
library(RMINC)
library(data.table)
library(matrixStats)

# Arg 1: Micro name
args = commandArgs(trailingOnly=TRUE)

# args = c()
# args[1] = "FA"
# args[2] = "../subj_zscores_common_dx/results/ses2_Label_whole_brain_dx.tsv"
# args[3] = "../subj_zscores_common_dx/results/zscores_anlm_FA.tsv"
# args[4] = "../../../WMH_micro_spatial/QC/inclusions_only_dx_new.txt"
# args[5] = "./results/raw/perc_abnormal_dx_FA_anlm.tsv"

tissue_nm = c('Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM')
tissue_all=c('Ventricules', 'CSF', 'Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH')

# names = c("FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM")

# Load tab data
print("Load tab data")
inclusions = as.data.frame(fread(args[4]))
colnames(inclusions) = c("ID")

demo = as.data.frame(fread("../../../UKB/tabular/df_demo/UKBB_demo_wider.tsv"))
demo = subset(demo, InstanceID == 2, select=c('SubjectID', 'Sex_31_0', 'Age_when_attended_assessment_centre_21003_0'))
colnames(demo) = c("ID", "Sex", "Age")
demo = merge(inclusions, demo, by="ID", all.x=TRUE)

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
print(dim(inclusions))

print("Run for each ID")

results = perc_values_above_thresholds(demo$ID[1], as.numeric(micro[1,]), as.numeric(label[1,]))

for (id in 2:nrow(micro)) {
  print(demo$ID[id])
  
  id_result = perc_values_above_thresholds(demo$ID[id], as.numeric(micro[id,]), as.numeric(label[id,]))

  results = rbind(results, id_result)
}

fwrite(as.data.frame(results), args[5], row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")




# dx = as.data.frame(fread("../../../UKB/QC/exclusion_lists/dx_single.csv"))
# colnames(dx) = c("Stroke", "TIA", "Subdural_H", "Subarachnoid_H", "Head_trauma", "Psychiatric", "Infection_nerv", "Abscess", 
#                 "Encephalitis", "Meningitis", "Chronic_neurol", "Motor_neuron_disease", "MS", "PD", "Cog_imp", "Epilepsy",
#                 "Head_injury", "Alcoholism", "Opioid_dep", "Other_dep", "Other_neurol", "Haemorrhage", "Ischaemic_stroke", "Fracture")

# # Max number of dx for one person
# max_occurrences <- function(row) {
#   max(table(row))
# }
# max_counts <- apply(dx, 1, max_occurrences)
# max_dx <- max(max_counts)

# # Merge demo and dx
# demo$dx_1 = NA
# demo$dx_2 = NA

# for (c in 1:ncol(dx)) {
#   # print(colnames(dx)[c])
#   for (id in 1:nrow(dx)) {
#     # print(dx[id,c])
#     row_demo = which(demo$ID %in% dx[id,c])
#     if (length(row_demo) > 0) {
#       if (is.na(demo$dx_1[row_demo])==TRUE) {
#         demo$dx_1[row_demo] = colnames(dx)[c]
#       } else {
#         demo$dx_2[row_demo] = colnames(dx)[c]
#       }
#     }
#   }
# }