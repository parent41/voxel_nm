
library(data.table)
library(dplyr)
library(tidyverse)

# Load data

inc = as.data.frame(fread("../../WMH_micro_spatial/QC/inclusions_without_excluding_dx_new.txt"))
inc = as.numeric(inc$V1)

demo = as.data.frame(fread("../../UKB/tabular/df_demo/UKBB_demo_wider.tsv"))
demo = subset(demo, InstanceID == 2 & is.na(Age_when_attended_assessment_centre_21003_0) == FALSE, select=c("SubjectID", "Sex_31_0", "Age_when_attended_assessment_centre_21003_0"))
colnames(demo) = c("ID", "sex", "age")

dx_categ = as.data.frame(fread("../../WMH_micro_spatial/Analyses_nm/predict_firstocc_categ/results/ICD0_codes_categ.csv"))

firstocc = as.data.frame(fread("../../WMH_micro_spatial/Analyses_nm/predict_firstocc_categ/results/firstocc_categ.tsv"))

# Add ppl without any dx
no_dx_ids = as.data.frame(inc[!(inc %in% firstocc$ID)])
colnames(no_dx_ids) = "ID"
firstocc = merge(no_dx_ids, firstocc, by="ID", all=TRUE)
firstocc = merge(demo, firstocc, by="ID")

hc = firstocc %>% filter(is.na(icd_code) == FALSE)
dx = firstocc %>% filter(is.na(icd_code) == FALSE)

# Only dx to include in NMs are hemmoroids, depressive episode, other anxiety disorders, migraine
dx_not_excluded = c("I84", "F32", "F41", "G43")

firstocc_final = firstocc %>%
    mutate(icd_code = ifelse(icd_code %in% dx_not_excluded | is.na(icd_code) == TRUE, NA, icd_code)) %>%
    mutate(any_dx = factor(ifelse(is.na(icd_code) == TRUE, "HC", "DX"))) %>%
    glimpse()

ggplot(firstocc_final, aes(x=age, color=any_dx)) + 
    geom_density() +
    # scale_color_discrete(labels=c("HC", "DX")) + 
    theme_classic() +
    theme(text = element_text(size=20), plot.title = element_text(hjust=0.5))
ggsave("age_dist.png")

# Write list of inclusions for running NM
inc_ids = firstocc_final %>% filter(is.na(icd_code) == TRUE) %>% select("ID")
inc_ids = as.data.frame(unique(inc_ids$ID))

fwrite(inc_ids, "inclusions_nodx_inputNM.txt", row.names=FALSE, col.names=FALSE)
