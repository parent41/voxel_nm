
library(data.table)
library(RMINC)
library(splitstackshape)
library(ggridges)
library(ggplot2)
library(scales)
library(tidyverse)

names = c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")

# Load data

inclusions = as.data.frame(fread("../../../WMH_micro_spatial/QC/inclusions_only_dx_new.txt"))
colnames(inclusions) = "ID"

demo = as.data.frame(fread("../../../UKB/tabular/df_demo/UKBB_demo_wider.tsv"))
demo = subset(demo, InstanceID == 2, select=c('SubjectID', 'Sex_31_0', 'Age_when_attended_assessment_centre_21003_0'))
colnames(demo) = c("ID", "Sex", "Age")
demo = merge(inclusions, demo, by="ID", all=FALSE)

dx = as.data.frame(fread("../../../UKB/QC/exclusion_lists/dx_single.csv"))
colnames(dx) = c("Stroke", "TIA", "Subdural_H", "Subarachnoid_H", "Head_trauma", "Psychiatric", "Infection_nerv", "Abscess", 
                "Encephalitis", "Meningitis", "Chronic_neurol", "Motor_neuron_disease", "MS", "PD", "Cog_imp", "Epilepsy",
                "Head_injury", "Alcoholism", "Opioid_dep", "Other_dep", "Other_neurol", "Haemorrhage", "Ischaemic_stroke", "Fracture")

dx <- dx[, colSums(!is.na(dx)) > 0]

# Load MRI data
nchunks = 32

label=as.data.frame(fread(paste0("./results/ses2_Label_whole_brain_dx.tsv")))

zscores_anlm_all = list()

for (n in 1:length(names)) {
    print(names[n])
    zscores_anlm_all[[n]] = as.data.frame(fread(paste0("./results/zscores_anlm_",names[n],".tsv")))
}


# Make density plots

tissue=c('Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH')

ids = as.numeric(inclusions$ID)

color_scale = c(MD = "#04319E", ISOVF = "#8298CF", FA = "#448312", ICVF="#2B520B", OD="#A2C189", T2star="#D36108", QSM="#E9B084")

for (d in 1:ncol(dx)) {
    cat(paste0("\nDiagnosis = ", colnames(dx)[d]))
    dx_ids = dx[,d]
    dx_ids = na.omit(dx_ids)

    dir.create(paste0("./visualization/",colnames(dx)[d]), showWarnings=FALSE)

    for (i in 1:length(dx_ids)) {
        print(dx_ids[i])
        demo_id = demo[which(demo$ID == dx_ids[i]),]
        vis_dir=paste0("./visualization/",colnames(dx)[d], "/", dx_ids[i])

        df = list()
        for (n in 1:length(names)) {
            print(names[n])
            df[[n]] = as.data.frame(cbind(as.numeric(zscores_anlm_all[[n]][which(demo$ID == dx_ids[i]),]), as.numeric(label[which(demo$ID == dx_ids[i]),])))
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
            scale_y_continuous(limits = c(-2, 2), name=paste0("Z-scores"), breaks = seq(-2,3)) +
            theme_classic() + 
            theme(text=element_text(size=20), axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
        ggsave(paste0(vis_dir, "/", dx_ids[i], "_zscore_anlm_hist.png"), width=2500, height=1500, dpi=300, units = "px")
        print(paste0(vis_dir, "/", dx_ids[i], "_zscore_anlm_hist.png"))

        command = paste0("convert -gravity East ",vis_dir,"/*_anlm_zscore_max3.png ",vis_dir,"/*_zscore_anlm_hist.png -append ",vis_dir, "/", dx_ids[i], "_zscore_anlm_all.png")
        system(command)
        print(paste0(vis_dir, "/", dx_ids[i], "_zscore_anlm_all.png"))
    }
}

