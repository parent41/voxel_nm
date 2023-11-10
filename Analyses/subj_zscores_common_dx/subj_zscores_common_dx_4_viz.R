
library(data.table)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(RMINC)
library(MRIcrotome)
library(viridis)
library(RColorBrewer)
library(foreach)
library(doParallel)

names = c("FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM")
names_long = c("dti_FA", "dti_MD", "NODDI_ICVF", "NODDI_ISOVF", "NODDI_OD", "T2star", "QSM")

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

anatVol = mincArray(mincGetVolume("../../../UKB/temporary_template/avg.020_2mm.mnc"))
mask = mincGetVolume("../../../UKB/temporary_template/Mask_2mm_dil2.mnc")

# Make mnc and png files for each subject grouped by dx

zscore_png = function(input, output, up_thresh, down_thresh) {
    input_id = sub(".*/", "", input)

    color_scale_div_2 = colorRampPalette(brewer.pal(9,"Blues"))(255)
    color_scale_div_1 = colorRampPalette(brewer.pal(9,"Reds"))(255)

    anatVol = mincArray(mincGetVolume(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-",input_id,"_ses-2_flair_UKB.mnc")))
    micro_subj = list()
    for (n in 1:length(names)) {
        # print(names[n])
        micro_subj[[n]] = mincArray(mincGetVolume(paste0(input, "_",names[n],".mnc")))
    }

    png(file=paste0(output, "_all_zscore_min",down_thresh,"_max",up_thresh,".png"), width=8500, height=4000, pointsize = 150)
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle("FLAIR") %>%
        anatomy(anatVol, low=10, high=200) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[1]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(micro_subj[[1]], low=down_thresh, high=up_thresh,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[2]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(micro_subj[[2]], low=down_thresh, high=up_thresh,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[3]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(micro_subj[[3]], low=down_thresh, high=up_thresh,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[4]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(micro_subj[[4]], low=down_thresh, high=up_thresh,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[5]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(micro_subj[[5]], low=down_thresh, high=up_thresh,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[6]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(micro_subj[[6]], low=down_thresh, high=up_thresh,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[7]) %>%
        anatomy(anatVol, low=10, high=200) %>%
        overlay(micro_subj[[7]], low=down_thresh, high=up_thresh,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
    legend("Z-values") %>%
    draw()
    dev.off()
}

raw_png = function(input, output) {
    input_id = sub(".*/", "", input)

    color_scale_div_2 = colorRampPalette(brewer.pal(9,"Blues"))(255)
    color_scale_div_1 = colorRampPalette(brewer.pal(9,"Reds"))(255)

    down_lims = c(0.05,0.05,0.05,0.05,0.05,0.7,0.001)
    up_lims = c(0.95,0.95,0.95,0.95,0.95,0.99,0.999)

    anatVol = mincArray(mincGetVolume(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-",input_id,"_ses-2_flair_UKB.mnc")))
    micro_subj = list()
    for (n in 1:length(names)) {
        # print(names[n])
        micro_subj[[n]] = mincArray(mincGetVolume(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-",input_id,"_ses-2_",names_long[n],"_UKB.mnc")))
    }

    png(file=paste0(output, "_all_raw.png"), width=8500, height=4000, pointsize = 150)
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle("FLAIR") %>%
        anatomy(anatVol, low=10, high=200) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[1]) %>%
        anatomy(micro_subj[[1]], low=quantile(micro_subj[[1]][][micro_subj[[1]][]!=0],down_lims[1]), high=quantile(micro_subj[[1]][][micro_subj[[1]][]!=0], up_lims[1])) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[2]) %>%
        anatomy(micro_subj[[2]], low=quantile(micro_subj[[2]][][micro_subj[[2]][]!=0],down_lims[2]), high=quantile(micro_subj[[2]][][micro_subj[[2]][]!=0], up_lims[2])) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[3]) %>%
        anatomy(micro_subj[[3]], low=quantile(micro_subj[[3]][][micro_subj[[3]][]!=0],down_lims[3]), high=quantile(micro_subj[[3]][][micro_subj[[3]][]!=0], up_lims[3])) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[4]) %>%
        anatomy(micro_subj[[4]], low=quantile(micro_subj[[4]][][micro_subj[[4]][]!=0],down_lims[4]), high=quantile(micro_subj[[4]][][micro_subj[[4]][]!=0], up_lims[4])) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[5]) %>%
        anatomy(micro_subj[[5]], low=quantile(micro_subj[[5]][][micro_subj[[5]][]!=0],down_lims[5]), high=quantile(micro_subj[[5]][][micro_subj[[5]][]!=0], up_lims[5])) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[6]) %>%
        anatomy(micro_subj[[6]], low=quantile(micro_subj[[6]][][micro_subj[[6]][]!=0],down_lims[6]), high=quantile(micro_subj[[6]][][micro_subj[[6]][]!=0], up_lims[6])) %>%
    sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
        addtitle(names[7]) %>%
        anatomy(micro_subj[[7]], low=quantile(micro_subj[[7]][][micro_subj[[7]][]!=0],down_lims[7]), high=quantile(micro_subj[[7]][][micro_subj[[7]][]!=0], up_lims[7])) %>%
    legend("Micro") %>%
    draw()
    dev.off()
}

# Run
cores = detectCores()
parallelCluster <- makeCluster(cores-5, type = "SOCK", methods = FALSE)
setDefaultCluster(parallelCluster)
registerDoParallel(parallelCluster)

down_thresholds <<- c(0.015, 1, 2, 3)
up_thresholds <<- c(3,3,4,5)

for (d in 1:ncol(dx)) {
    cat(paste0("\nDiagnosis = ", colnames(dx)[d]))
    dx_ids = dx[,d]
    dx_ids = na.omit(dx_ids)

    dir.create(paste0("./visualization/",colnames(dx)[d]), showWarnings=FALSE)

    foreach(i=1:length(dx_ids), .packages = c('RMINC', 'MRIcrotome', 'RColorBrewer', 'tidyverse')) %dopar% {
    # for (i in 1:length(dx_ids)) {
        cat(paste0("\n\tID = ", dx_ids[i], "\t"))
        id_row = which(inclusions$ID == dx_ids[i])

        # If dx ID is in inclusions
        if(length(id_row)>0) {

            res_dir = paste0("./results/",colnames(dx)[d], "/", dx_ids[i])
            viz_dir = paste0("./visualization/",colnames(dx)[d], "/", dx_ids[i])
            dir.create(viz_dir, showWarnings=FALSE)
            
            # Make PNG
            input = paste0(res_dir, "/", dx_ids[i])
            output = paste0(viz_dir, "/", dx_ids[i])
            for (t in 1:length(down_thresholds)) {
                zscore_png(input, output, up_thresholds[t], down_thresholds[t])
                cat(paste0("\n\tZ-score PNG = ", output, "_all_zscore_min",down_thresholds[t],"_max",up_thresholds[t],".png"))
            }
            raw_png(input, output)
            cat(paste0("\n\tRaw PNG = ", output, "_all_raw.png"))

        }
    }
}

stopCluster(parallelCluster)




