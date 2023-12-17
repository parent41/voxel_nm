
library(data.table)
library(RMINC)
library(splitstackshape)
library(ggridges)
library(ggplot2)
library(scales)
library(tidyverse)

# Arg 1: Micro name
# args = commandArgs(trailingOnly=TRUE)

# args = c()
# args[1] = "FA"

names = c("FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM")
names = c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")

# Load data

inclusions = as.data.frame(fread("../../../WMH_micro_spatial/QC/inclusions_excluding_dx_new.txt"))
colnames(inclusions) = "ID"

demo = as.data.frame(fread("../../../UKB/tabular/df_demo/UKBB_demo_wider.tsv"))
demo = subset(demo, InstanceID == 2, select=c('SubjectID', 'Sex_31_0', 'Age_when_attended_assessment_centre_21003_0'))
colnames(demo) = c("ID", "Sex", "Age")
demo = merge(inclusions, demo, by="ID", all=FALSE)

# Load MRI data
nchunks = 32

label=as.data.frame(fread(paste0("./tmp/sampled_label_",names[1],".tsv")))
ids_label = as.data.frame(fread(paste0("./tmp/ids_sampled_label_",names[1],".txt")))

zscores_anlm_all = list()
ids_micro_all = list()

for (n in 1:length(names)) {
    print(names[n])
    zscores_anlm_list = list()
    ids_micro_list = list()

    for(i in 0:nchunks) {
        print(i)
        # Load matrices
        zscores_anlm_list[[i+1]] = as.data.frame(fread(paste0("./tmp/sampled_label_c",i,"_",names[n],"_anlm.tsv")))
        ids_micro_list[[i+1]] = as.data.frame(fread(paste0("./tmp/ids_sampled_micro_c",i,"_",names[n],".txt")))
    }
    
    zscores_anlm_all[[n]] <- do.call(rbind, zscores_anlm_list)
    ids_micro_all[[n]] = do.call(rbind, ids_micro_list)
}


# Make density plots

tissue=c('Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH')

ids = as.numeric(ids_micro_all[[n]]$V1)

color_scale = c(MD = "#04319E", ISOVF = "#8298CF", FA = "#448312", ICVF="#2B520B", OD="#A2C189", T2star="#D36108", QSM="#E9B084")

for (i in 1:length(ids)) {
    print(ids[i])
    demo_id = demo[which(demo$ID == ids[i]),]
    vis_dir=paste0("./visualization/",demo_id$Age,ifelse(demo_id$Sex == "Male", "M", "F"),"_",ids[i])

    df = list()
    for (n in 1:length(names)) {
        print(names[n])
        df[[n]] = as.data.frame(cbind(as.numeric(zscores_anlm_all[[n]][which(demo$ID == ids[i]),]), as.numeric(label[which(demo$ID == ids[i]),])))
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
        scale_fill_manual(name="", values=color_scale) +
        scale_x_discrete(labels = tissue, name="") +
        scale_y_continuous(limits = c(-3, 3), name=paste0("Z-scores"), breaks = seq(-3,3)) +
        theme_classic() + 
        theme(text=element_text(size=20), axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
    ggsave(paste0(vis_dir, "/", demo_id$Age,ifelse(demo_id$Sex == "Male", "M", "F"),"_",ids[i], "_zscore_anlm_hist.png"), width=2500, height=1500, dpi=300, units = "px")

    command = paste0("convert -gravity East ",vis_dir,"/*_zscore_anlm_max3.png ",vis_dir,"/*_zscore_anlm_hist.png -append ",vis_dir, "/", demo_id$Age,ifelse(demo_id$Sex == "Male", "M", "F"),"_",ids[i], "_zscore_anlm_all.png")
    system(command)
}



