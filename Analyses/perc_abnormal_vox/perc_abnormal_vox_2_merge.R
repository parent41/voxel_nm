
library(RMINC)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(data.table)
library(ggpubr)


tissue_nm = c('Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM')
tissue_all=c('Ventricules', 'CSF', 'Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH')

names = c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")

# Load tab data
print("Load tab data")
inclusions = as.data.frame(fread("../../../WMH_micro_spatial/QC/inclusions_excluding_dx_new.txt"))
colnames(inclusions) = "ID"
inclusions_dx = as.data.frame(fread("../../../WMH_micro_spatial/QC/inclusions_only_dx_new.txt"))
colnames(inclusions_dx) = "ID"

demo = as.data.frame(fread("../../../UKB/tabular/df_demo/UKBB_demo_wider.tsv"))
demo = subset(demo, InstanceID == 2, select=c('SubjectID', 'Sex_31_0', 'Age_when_attended_assessment_centre_21003_0'))
colnames(demo) = c("ID", "Sex", "Age")
demo = merge(inclusions, demo, by="ID", all.x=TRUE)

demo$dx_1 = NA
demo$dx_2 = NA

# Load and merge dx data

demo_dx = as.data.frame(fread("../../../UKB/tabular/df_demo/UKBB_demo_wider.tsv"))
demo_dx = subset(demo_dx, InstanceID == 2, select=c('SubjectID', 'Sex_31_0', 'Age_when_attended_assessment_centre_21003_0'))
colnames(demo_dx) = c("ID", "Sex", "Age")
demo_dx = merge(inclusions_dx, demo_dx, by="ID", all.x=TRUE)

dx = as.data.frame(fread("../../../UKB/QC/exclusion_lists/dx_single.csv"))
colnames(dx) = c("Stroke", "TIA", "Subdural_H", "Subarachnoid_H", "Head_trauma", "Psychiatric", "Infection_nerv", "Abscess", 
                "Encephalitis", "Meningitis", "Chronic_neurol", "Motor_neuron_disease", "MS", "PD", "Cog_imp", "Epilepsy",
                "Head_injury", "Alcoholism", "Opioid_dep", "Other_dep", "Other_neurol", "Haemorrhage", "Ischaemic_stroke", "Fracture")

# Max number of dx for one person
max_occurrences <- function(row) {
  max(table(row))
}
max_counts <- apply(dx, 1, max_occurrences)
max_dx <- max(max_counts)

# Merge demo and dx
demo_dx$dx_1 = NA
demo_dx$dx_2 = NA

for (c in 1:ncol(dx)) {
  # print(colnames(dx)[c])
  for (id in 1:nrow(dx)) {
    # print(dx[id,c])
    row_demo = which(demo_dx$ID %in% dx[id,c])
    if (length(row_demo) > 0) {
      if (is.na(demo_dx$dx_1[row_demo])==TRUE) {
        demo_dx$dx_1[row_demo] = colnames(dx)[c]
      } else {
        demo_dx$dx_2[row_demo] = colnames(dx)[c]
      }
    }
  }
}

# Merge hc and dx

demo_all = rbind(demo, demo_dx)

demo_all[,c(2,4,5)] = as.data.frame(lapply(demo_all[,c(2,4,5)], as.factor))

# Load perc abnormal vox results

results_fl = list.files("./results", pattern="*", full.names = TRUE)

results = as.data.frame(fread(results_fl[1]))
micro_res = sub(".*_(.*)_anlm.tsv", "\\1", results_fl[1])
results$micro = micro_res

for (f in 2:length(results_fl)) {
    print(results_fl[f])
    results_tmp = as.data.frame(fread(results_fl[f]))
    micro_res = sub(".*_(.*)_anlm.tsv", "\\1", results_fl[f])
    results_tmp$micro = micro_res
    results = rbind(results, results_tmp)
}

results = results[order(results$ID),]

# Remove dx with prevalence less than 30

dx_prev = as.data.frame(table(demo_dx$dx_1))
dx_prev_order = dx_prev[order(-dx_prev$Freq),]
dx_prev_order = subset(dx_prev_order, Freq >= 30, select="Var1")
dx_prev_order = as.character(dx_prev_order$Var1)
dx_toremove = as.character(dx_prev[which(dx_prev$Freq < 30), 1])

demo_dx_tobind = subset(demo_dx, is.na(dx_2) == FALSE)
demo_dx_tobind$dx_1 = demo_dx_tobind$dx_2
demo_dx_tobind = demo_dx_tobind[,c(1,2,3,4)]
colnames(demo_dx_tobind)[ncol(demo_dx_tobind)] = "dx"
demo_dx_tobind$dx_num = 2

demo_dx_prev = demo_dx
demo_dx_prev = demo_dx_prev[,c(1,2,3,4)]
colnames(demo_dx_prev)[ncol(demo_dx_prev)] = "dx"
demo_dx_prev$dx_num = 1

demo_dx_prev = rbind(demo_dx_prev, demo_dx_tobind)
demo_dx_prev$dx[which(demo_dx_prev$dx %in% dx_toremove)] <- NA
demo_dx_prev = subset(demo_dx_prev, is.na(dx) == FALSE)

ids_toremove = inclusions_dx$ID[-which(inclusions_dx$ID %in% demo_dx_prev$ID)]

# Merge results with dx

results$dx = NA
results = results[-which(results$ID %in% ids_toremove),]

results_dx = results
for (i in 1:nrow(demo_dx_prev)) {
  print(demo_dx_prev$ID[i])
  print(demo_dx_prev$dx_num[i])
  print(demo_dx_prev$dx[i])

  if (demo_dx_prev$dx_num[i] == 1) {
    results_dx[which(results_dx$ID == demo_dx_prev$ID[i]),]$dx <- demo_dx_prev$dx[i]
  } else if (demo_dx_prev$dx_num[i] == 2) {
    # Important: subjects with two neurological dx appear twice
    results_tmp = results_dx[which(results_dx$ID == demo_dx_prev$ID[i]),]
    results_tmp$dx = demo_dx_prev$dx[i]
    results_dx = rbind(results_dx, results_tmp)
  }
}

results_dx$dx[is.na(results_dx$dx)] <- "HC"

fwrite(results_dx, "./results/final_results_dx.tsv", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

