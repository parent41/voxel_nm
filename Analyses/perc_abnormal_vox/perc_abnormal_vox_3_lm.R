
library(RMINC)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(data.table)
library(ggpubr)
library(scales)
# install.packages(c("fmsb"))
# install.packages("devtools")
# devtools::install_github("ricardo-bion/ggradar", dependencies = TRUE, lib="/gpfs/fs1/home/m/mchakrav/parent41/R/x86_64-pc-linux-gnu-library/4.1")
# library(fmsb)
# library(ggradar)

tissue_nm = c('Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM')
tissue_all=c('Ventricules', 'CSF', 'Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH', 'Cerebral_WM')
tissue_abn = c('Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH', 'Cerebral_WM')

names = c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM", "jacobians_rel", "jacobians_abs")

# Load data
results = as.data.frame(fread("./results/results.tsv"))
inclusions = as.data.frame(fread("../../../WMH_micro_spatial/QC/inclusions_without_excluding_dx_new.txt"))
colnames(inclusions) = "ID"

demo = as.data.frame(fread("../../../UKB/tabular/df_demo/UKBB_demo_wider.tsv"))
demo = subset(demo, InstanceID == 2, select=c('SubjectID', 'Sex_31_0', 'Age_when_attended_assessment_centre_21003_0'))
colnames(demo) = c("ID", "Sex", "Age")
demo = merge(inclusions, demo, by="ID", all.x=TRUE)

df_dx = as.data.frame(fread("../../../WMH_micro_spatial/Analyses_nm/predict_firstocc_categ/results/firstocc_categ.tsv"))

no_dx_ids = as.data.frame(inclusions$ID[!(inclusions$ID %in% df_dx$ID)])
colnames(no_dx_ids) = "ID"
df_dx = merge(no_dx_ids, df_dx, by="ID", all=TRUE)

df = merge(demo, df_dx, all.y=TRUE, by="ID")
df[,c(2,5,6,9,10)] = as.data.frame(lapply(df[,c(2,5,6,9,10)], as.factor))

results[,c(2,3,7)] = as.data.frame(lapply(results[,c(2,3,7)], as.factor))

controls = as.data.frame(fread("../../QC/inclusions_nodx_inputNM.txt"))
controls = controls$V1

# Function to calculate LM between dx and controls for % of abnormal voxels
lm_perc_abn = function(dx_icd, dx_name) {
  dx_ids = df$ID[which(df$icd_code %in% dx_icd)]

  # Check if n(dx) is higher than 0
  if (length(dx_ids) > 0) {

    # Prep dataframe for linear model
    dx_tolm = subset(results, ID %in% dx_ids)
    dx_tolm = results[which(results$ID %in% dx_ids),]
    hc_tolm = results[which(results$ID %in% controls),]

    dx_tolm$group = "DX"
    hc_tolm$group = "HC"

    df_tolm = rbind(dx_tolm, hc_tolm)
    df_tolm = merge(demo, df_tolm, all.y=TRUE, by="ID")

    df_tolm[,c(2,10)] = as.data.frame(lapply(df_tolm[,c(2,10)], as.factor))
    df_tolm$group <- relevel(df_tolm$group, ref = "HC")

    # Run linear models for each combination of label, threshold, and micro
    final_res_list = list()

    i=1
    for (l in levels(df_tolm$label_value)) {
    # print(l)
      for (t in levels(df_tolm$threshold)) {
        # print(t)
        for (m in levels(df_tolm$micro)) {
          # print(m)

          lm_df = subset(df_tolm, label_value == l & threshold == t & micro == m)

          lm_res = lm(scale(perc_vox_above_thresh) ~ group + scale(Age) + Sex, data=lm_df)
          sum_lm = summary(lm_res)
          coefs = as.data.frame(t(unlist(as.data.frame(sum_lm$coefficients), use.names=TRUE)))
          colnames(coefs) = c("Estimate_intercept", "Estimate_group", "Estimate_Age", "Estimate_Sex",
                          "stderr_intercept", "stderr_group", "stderr_Age", "stderr_Sex",
                          "tval_intercept", "tval_group", "tval_Age", "tval_Sex",
                          "pval_intercept", "pval_group", "pval_Age", "pval_Sex")

          final_res_list[[i]] = data.frame(dx=dx_name, Label=l, Threshold=t, Micro=m)
          final_res_list[[i]] = cbind(final_res_list[[i]], coefs)

          i=i+1
        }
      }
    }

    # Return all coefficients from linear models
    final_res = do.call(rbind, final_res_list)
    
  } else {
    final_res = data.frame(dx=dx_name, Label=NA, Threshold=NA, Micro=NA, Estimate_intercept=NA, Estimate_group=NA, Estimate_Age=NA, Estimate_Sex=NA,
                          stderr_intercept=NA, stderr_group=NA, stderr_Age=NA, stderr_Sex=NA,
                          tval_intercept=NA, tval_group=NA, tval_Age=NA, tval_Sex=NA,
                          pval_intercept=NA, pval_group=NA, pval_Age=NA, pval_Sex=NA)
  }
  
  return(final_res)
}

# Results for all categories of diagnoses

icd_codes_list <- readRDS("../../../WMH_micro_spatial/Analyses_nm/predict_firstocc_categ/results/icd_codes_list.rds")

res_dx_list = list()
for (dx in 2:length(icd_codes_list)) {
  print(names(icd_codes_list)[dx])
  res_dx_list[[dx-1]] = lm_perc_abn(icd_codes_list[[dx]],names(icd_codes_list)[dx])
}

res_dx = do.call(rbind, res_dx_list)
fwrite(res_dx, "./results/lm_results.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

