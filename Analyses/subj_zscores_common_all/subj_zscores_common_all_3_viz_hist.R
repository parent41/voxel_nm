
library(data.table)
library(RMINC)
library(splitstackshape)
library(ggridges)
library(ggplot2)
library(scales)
library(tidyverse)

# Arg 1: Micro name
args = commandArgs(trailingOnly=TRUE)
num_chunk = as.numeric(args[2])

# args = c()
# num_chunk = 1

names = c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM", "jacobians_abs", "jacobians_rel")
names_toplot = c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM", "J(abs)", "J(rel)")

# Load data

inclusions = as.data.frame(fread("../../../WMH_micro_spatial/QC/inclusions_without_excluding_dx_new.txt"))
colnames(inclusions) = "ID"

demo = as.data.frame(fread("../../../UKB/tabular/df_demo/UKBB_demo_wider.tsv"))
demo = subset(demo, InstanceID == 2, select=c('SubjectID', 'Sex_31_0', 'Age_when_attended_assessment_centre_21003_0'))
colnames(demo) = c("ID", "Sex", "Age")
demo = merge(inclusions, demo, by="ID", all=FALSE)

dx = as.data.frame(fread("../../../WMH_micro_spatial/Analyses_nm/predict_firstocc_categ/results/firstocc_categ.tsv"))
dx = merge(demo, dx, by="ID", all=TRUE)
dx[,c(2,5,6,9,10)] = as.data.frame(lapply(dx[,c(2,5,6,9,10)], as.factor))

# Load MRI data

ids_label = as.data.frame(fread(paste0("./tmp/ids_label_c",num_chunk,"_FA.txt")))
colnames(ids_label) = "ID"
ids_micro_all = list()

for (n in 1:length(names)) {
    print(names[n])
    ids_micro_all[[n]] = as.data.frame(fread(paste0("./tmp/ids_micro_c",num_chunk,"_",names[n],"_anlm.txt")))
}

# Make density plots
ids = as.numeric(ids_micro_all[[1]]$V1)

tissue=c('Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH')
color_scale = c(MD = "#04319E", ISOVF = "#8298CF", FA = "#448312", ICVF="#2B520B", OD="#A2C189", T2star="#D36108", QSM="#E9B084", jacobians_abs = "#00A7B9", jacobians_rel = "#47EDFF")
mask = mincGetVolume("../../../UKB/temporary_template/Mask_2mm_dil2.mnc")

for (i in 1:length(ids)) {
    print(ids[i])
    demo_id = demo[which(demo$ID == ids[i]),]
    dx_id = dx[which(dx$ID %in% ids[i]),]

    id_string = paste0(ids[i], "_", dx_id$Age[1], ifelse(dx_id$Sex[1] == "Female", "F", "M"))
    if(is.na(dx_id$icd_code[1]) == FALSE) {
        for (d in 1:length(dx_id$icd_code)) {
            id_string = paste0(id_string, "_", dx_id$icd_code[d])
        }
    }

    vis_dir=paste0("./visualization/",id_string)
    dir.create(vis_dir, showWarnings=FALSE)

    label_id = mincGetVolume(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-",ids[i],"_ses-2_Label_UKB.mnc"))
    df = list()
    for (n in 1:length(names)) {
        print(names[n])
        micro_id = mincGetVolume(paste0("./results/mnc/",ids[i],"_zscore_anlm_",names[n],".mnc"))
        df[[n]] = as.data.frame(cbind(as.numeric(micro_id[][which(mask[]>0.5)]), round(as.numeric(label_id[][which(mask[]>0.5)]),0)))
        colnames(df[[n]]) = c("Values", "Label")
        df[[n]]$Label = as.factor(df[[n]]$Label)
        df[[n]] = df[[n]][which(df[[n]]$Label %in% c(3,4,5,6,7,8,9)),]
        df[[n]] = df[[n]][complete.cases(df[[n]]),]
        df[[n]]$Micro = names[n]
    }
    df_id = do.call(rbind, df)
    df_id$Micro = as.factor(df_id$Micro)

    ggplot(df_id, aes(x=Label, y=Values, fill=factor(Micro, levels=names))) +
        geom_boxplot(outlier.shape = NA) + 
        geom_hline(yintercept = 0) + 
        scale_fill_manual(name="", values=color_scale, labels=names_toplot) +
        scale_x_discrete(labels = tissue, name="") +
        scale_y_continuous(limits = c(-2, 2), name=paste0("\nZ-scores"), breaks = seq(-2,3)) +
        theme_classic() + 
        theme(text=element_text(size=20), axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
    ggsave(paste0(vis_dir, "/", id_string, "_zscore_anlm_hist.png"), width=2500, height=1500, dpi=300, units = "px")
    print(paste0(vis_dir, "/", id_string, "_zscore_anlm_hist.png"))

    # command = paste0("convert -gravity East ",vis_dir,"/*_zscore_anlm_max3.png ",vis_dir,"/*_zscore_anlm_hist.png -append ",vis_dir, "/", demo_id$Age,ifelse(demo_id$Sex == "Male", "M", "F"),"_",ids[i], "_zscore_anlm_all.png")
    # system(command)
    # print(paste0(vis_dir, "/", demo_id$Age,ifelse(demo_id$Sex == "Male", "M", "F"),"_",ids[i], "_zscore_anlm_all.png"))
}



