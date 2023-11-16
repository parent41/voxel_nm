
library(data.table)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reticulate)
library(RMINC)
library(MRIcrotome)
library(viridis)
library(RColorBrewer)
library(foreach)
library(doParallel)

# Load data

names = c("FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM")
# names = c("FA")

tissues = c('Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM')
tissues_all = c('Ventricules', 'CSF', 'Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH')

nm = list()

mask = mincGetVolume("../../../UKB/temporary_template/Mask_2mm_dil2.mnc")

# Load NM
dir.create("./results/mnc", showWarnings=FALSE)
for (n in 1:length(names)) {
    print(names[n])
    dir.create(paste0("./results/mnc/",names[n]), showWarnings=FALSE)

    nm[[n]] = list()

    for (t in 1:length(tissues)) {
        print(tissues[t])
        dir.create(paste0("./results/mnc/",names[n],"/",tissues[t]), showWarnings=FALSE)

        nm[[n]][[t]] = as.data.frame(fread(paste0("./results/nm_",names[n],"_",tissues[t],".tsv")))
    }
}


# Find interesting voxels to visualize (ones in-between tissue types)

vox_tissue = list()

for (t in 1:length(tissues)) {
    vox_tissue[[t]] = subset(nm[[1]][[t]], N>10000, select=c("ROI", "N"))
}

vox_nawm_corticalGM = merge(vox_tissue[[5]], vox_tissue[[6]], by="ROI", all=FALSE)
vox_nawm_subcorticalGM = merge(vox_tissue[[4]], vox_tissue[[6]], by="ROI", all=FALSE)

set.seed(123)
random_samples = sort(c(sample(vox_nawm_corticalGM$ROI, 10), sample(vox_nawm_subcorticalGM$ROI, 10)))
fwrite(as.data.frame(random_samples), "./visualization/random_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Then run python script save_vox_micro_label.py to save matrices of those voxels only...

inclusions = as.data.frame(fread("../../../WMH_micro_spatial/QC/inclusions_excluding_dx_new.txt"))
colnames(inclusions) = "ID"

micro = list()

for (n in 1:length(names)) {
    print(names[n])

    micro[[n]] = as.data.frame(fread(paste0("./visualization/vox/",names[n],"_vox_to_viz_micro.tsv")))
    ids_micro = as.data.frame(fread(paste0("../../micro_matrices/ids_ses2_",names[n],".txt")))
    micro[[n]] = cbind(ids_micro, micro[[n]])
    colnames(micro[[n]])[1] = "ID"

    micro[[n]] = merge(inclusions, micro[[n]], by = "ID", all=FALSE)
}

label = as.data.frame(fread(paste0("./visualization/vox/vox_to_viz_label.tsv")))
ids_label = as.data.frame(fread("../../bison_matrices/ids_ses2_Label_whole_brain.txt"))
label = cbind(ids_label, label)
colnames(label)[1] = "ID"
label = merge(inclusions, label, by="ID", all=FALSE)

# Make plots of NM for certain voxels

demo = as.data.frame(fread("../../../UKB/tabular/df_demo/UKBB_demo_wider.tsv"))
demo = subset(demo, InstanceID == 2, select=c('SubjectID', 'Sex_31_0', 'Age_when_attended_assessment_centre_21003_0'))
colnames(demo) = c("ID", "Sex", "Age")
demo = merge(inclusions, demo, by="ID", all=FALSE)

anatVol = mincArray(mincGetVolume("../../../UKB/temporary_template/avg.020_2mm.mnc"))
mask = mincGetVolume("../../../UKB/temporary_template/Mask_2mm_dil2.mnc")

# Visualize a few voxels where there is tissue overlap

dir.create("./visualization/vox/mnc", showWarnings=FALSE)
for (i in 1:length(random_samples)) {
    plots=list()
    dir.create(paste0("./visualization/vox/vox_",random_samples[i]), showWarnings=FALSE)

    # Visualize micro and NM
    for (n in 1:length(names)) {
        print(names[n])

        df = as.data.frame(cbind(demo, micro[[n]][,i+1], label[,i+1]))
        colnames(df)[c(4,5)] = c("Value", "Label")
        df$Label = ifelse(is.na(df$Label), NA, tissues_all[df$Label])
        df$Label = as.factor(df$Label)

        nm_df = list()
        for (t in levels(df$Label)) {
            tissue = which(tissues == t)
            tissue_all = which(tissues_all == t)

            if (t %in% tissues) {
                nm_df_1 = nm[[n]][[tissue]][random_samples[i],]
                male_avg = as.numeric(nm_df_1[, grep("Male_mean", names(nm_df_1), value = TRUE)])
                male_sd = as.numeric(nm_df_1[, grep("Male_sd", names(nm_df_1), value = TRUE)])
                female_avg = as.numeric(nm_df_1[, grep("Female_mean", names(nm_df_1), value = TRUE)])
                female_sd = as.numeric(nm_df_1[, grep("Female_sd", names(nm_df_1), value = TRUE)])

                nm_df_2 = data.frame("Age" = seq(45,81), "Male_avg" = male_avg, "Male_SD" = male_sd, "Female_avg" = female_avg, "Female_SD" = female_sd)
                nm_df_2 = nm_df_2 %>%
                    pivot_longer(cols = -Age, 
                    names_to = c("Sex", ".value"), 
                    names_pattern = "(Male|Female)_(avg|SD)")

                nm_df[[as.character(tissues_all[tissue_all])]] = as.data.frame(nm_df_2)
            }
        }

        combined_nm_df = do.call(rbind, Map(cbind, nm_df, df_name = names(nm_df)))

        plots[[n]] = ggplot(df, aes(x=Age, y=Value, color=Label, shape=Sex)) +
                        geom_point(alpha=0.1, size=1) + 
                        geom_line(data = combined_nm_df, aes(x=Age, y=avg, linetype=Sex,color=df_name)) + 
                        ggtitle(names[n]) +
                        theme_classic() + 
                        theme(text=element_text(size=20), plot.title=element_text(hjust=0.5, size=30))
        ggsave(paste0("./visualization/vox/vox_",random_samples[i],"/",names[n],"_",random_samples[i],".png"), width=8, height=5)
        print(paste0("./visualization/vox/vox_",random_samples[i],"/",names[n],"_",random_samples[i],".png"))
    }

    wrap_plots(plots, ncol=3)
    ggsave(paste0("./visualization/vox/vox_",random_samples[i],"/all_",random_samples[i],".png"), width=24, height=15)

    # Visualize where voxels are located
    outvol <- mincGetVolume("../../../UKB/temporary_template/avg.020_2mm.mnc")
    outvol[] <- 0
    vox = c(rep(0, random_samples[i] - 1), 1, rep(0, (nrow(nm[[1]][[1]]) - (random_samples[i]))))
    outvol[mask > 0.5] <- vox
    mincWriteVolume(outvol, paste0("./visualization/vox/vox_",random_samples[i],"/vox_",random_samples[i],".mnc"), clobber=TRUE)
    command=paste0("mincmorph -clobber -3D26 -successive DD ./visualization/vox/vox_",random_samples[i],"/vox_",random_samples[i],".mnc ./visualization/vox/vox_",random_samples[i],"/vox_",random_samples[i],"_dil.mnc")
    system(command)

    vox_mask = mincArray(mincGetVolume(paste0("./visualization/vox/vox_",random_samples[i],"/vox_",random_samples[i],".mnc")))
    slice = which(vox_mask[] == 1, arr.ind = TRUE)[3]
    png(file=paste0("./visualization/vox/vox_",random_samples[i],"/vox_",random_samples[i],".png"), width=3000, height=3000, pointsize = 80)
    sliceSeries(nrow=1, ncol=1, begin=slice, end=slice, dimension=3) %>%
        #addtitle("On average") %>%
        anatomy(anatVol, low=10, high=140) %>%
        overlay(mincArray(mincGetVolume(paste0("./visualization/vox/vox_",random_samples[i],"/vox_",random_samples[i],".mnc"))),
                low=0, high=1, col="red") %>%
        draw(layout="row")
    dev.off()
    print(paste0("./visualization/vox/vox_",random_samples[i],"/vox_",random_samples[i],".png"))

    png(file=paste0("./visualization/vox/vox_",random_samples[i],"/vox_",random_samples[i],"_dil.png"), width=3000, height=3000, pointsize = 80)
    sliceSeries(nrow=1, ncol=1, begin=slice, end=slice, dimension=3) %>%
        #addtitle("On average") %>%
        anatomy(anatVol, low=10, high=140) %>%
        overlay(mincArray(mincGetVolume(paste0("./visualization/vox/vox_",random_samples[i],"/vox_",random_samples[i],"_dil.mnc"))),
                low=0, high=1, col="red") %>%
        draw(layout="row")
    dev.off()
    print(paste0("./visualization/vox/vox_",random_samples[i],"/vox_",random_samples[i],"_dil.png"))

    command = paste0("convert -gravity center ./visualization/vox/vox_",random_samples[i],"/vox_",random_samples[i],".png ./visualization/vox/vox_",random_samples[i],"/all_",random_samples[i],".png +append ./visualization/vox/vox_",random_samples[i],"/vox_all_",random_samples[i],".png")
    system(command)
}

