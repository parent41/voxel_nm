
library(data.table)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)
library(Hmisc)
library(corrplot)

names = c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM", "jacobians_rel", "jacobians_abs")
regions = c('Whole_brain', 'Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH', 'Cerebral_WM')

# Load data

inclusions = as.data.frame(fread(paste0("../../../WMH_micro_spatial/QC/inclusions_without_excluding_dx_new.txt")))
colnames(inclusions) = "ID"
ids = inclusions$ID

demo = as.data.frame(fread("../../../UKB/tabular/df_demo/UKBB_demo_wider.tsv"))
demo = subset(demo, InstanceID == 2, select=c('SubjectID', 'Sex_31_0', 'Age_when_attended_assessment_centre_21003_0'))
colnames(demo) = c("ID", "Sex", "Age")
demo = merge(inclusions, demo, by="ID", all=FALSE)

anatVol = mincArray(mincGetVolume("../../../UKB/temporary_template/avg.020_2mm.mnc"))
mask = mincGetVolume("../../../UKB/temporary_template/Mask_2mm.mnc")
mask_path = "../../../UKB/temporary_template/Mask_2mm.mnc"

dx = as.data.frame(fread("../../../WMH_micro_spatial/Analyses_nm/predict_firstocc_categ/results/firstocc_categ.tsv"))
dx = dx[which(dx$ID %in% ids),]
dx = merge(demo, dx, by="ID", all=TRUE)
dx[,c(2,5,6,9,10)] = as.data.frame(lapply(dx[,c(2,5,6,9,10)], as.factor))

hc_ids = as.data.frame(fread("../../QC/inclusions_nodx_inputNM.txt"))
hc_ids = hc_ids$V1

corr_r = as.data.frame(fread("./results/indiv_concat_raw.tsv"))
corr_z = as.data.frame(fread("./results/indiv_concat_zscores.tsv"))

corr_r[,c(1,2,3,6)] = as.data.frame(lapply(corr_r[,c(1,2,3,6)], as.factor))
corr_z[,c(1,2,3,6)] = as.data.frame(lapply(corr_z[,c(1,2,3,6)], as.factor))

# Mass univariate differences between dx and hc

for (dx in 1:length(levels())) {
    
}
