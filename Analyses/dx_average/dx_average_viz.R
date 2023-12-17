
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
library(ggnewscale)
library(scales)
library(stringr)

# Source yohan's ggplot for MRI extension
source("../../yohan_ggplot/plotting_functions.R")
source("../../yohan_ggplot/plotting_functions_labels.R")

dx = as.data.frame(fread("../../../UKB/QC/exclusion_lists/dx_single.csv"))
colnames(dx) = c("Stroke", "TIA", "Subdural_H", "Subarachnoid_H", "Head_trauma", "Psychiatric", "Infection_nerv", "Abscess", 
                "Encephalitis", "Meningitis", "Chronic_neurol", "Motor_neuron_disease", "MS", "PD", "Cog_imp", "Epilepsy",
                "Head_injury", "Alcoholism", "Opioid_dep", "Other_dep", "Other_neurol", "Haemorrhage", "Ischaemic_stroke", "Fracture")

dx <- dx[, colSums(!is.na(dx)) > 0]

anatVol = mincArray(mincGetVolume("../../../UKB/temporary_template/avg.020_2mm.mnc"))
mask = mincGetVolume("../../../UKB/temporary_template/Mask_2mm.mnc")
mask_path = "../../../UKB/temporary_template/Mask_2mm.mnc"

names = c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")

# Make png files for each dx average

zscore_png = function(input, output, up_thresh) {
    # Load data for slices across different axes
    slices = list()
    slices[[1]] = c(0,20) # Slices for z-axis
    slices[[2]] = c(-5, -30) # Slices for x-axis
    slices[[3]] = c(-20) # Slices for y-axis

    slices_axes = c("z", "x", "y")

    n_slices = 0

    micro_final = tibble(matrix(nrow = 0, ncol = 12))
    label_ukb_final = tibble(matrix(nrow = 0, ncol = 11))

    for (a in 1:length(slices)) {
        # print(slices_axes[a])
        slices_num = seq(n_slices + 1, n_slices + length(slices[[a]]))

        # Load image data

        micro = prepare_masked_anatomy(paste0(input, "_",names[1],"_anlm.mnc"), mask_path, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]])
        ifelse(slices_axes[a] == "x", micro$anatomy_df$x_toplot <- micro$anatomy_df$y, ifelse(slices_axes[a] == "y", micro$anatomy_df$x_toplot <- micro$anatomy_df$x, ifelse(slices_axes[a] == "z", micro$anatomy_df$x_toplot <- micro$anatomy_df$x, NA)))
        ifelse(slices_axes[a] == "x", micro$anatomy_df$y_toplot <- micro$anatomy_df$z, ifelse(slices_axes[a] == "y", micro$anatomy_df$y_toplot <- micro$anatomy_df$z, ifelse(slices_axes[a] == "z", micro$anatomy_df$y_toplot <- micro$anatomy_df$y, NA)))
        micro$anatomy_df$slice_index = micro$anatomy_df$slice_index + n_slices

        micro$anatomy_df$Micro = names[1]
        for (n in 2:length(names)) {
            micro_tmp = prepare_masked_anatomy(paste0(input, "_",names[n],"_anlm.mnc"), mask_path, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]])
            ifelse(slices_axes[a] == "x", micro_tmp$anatomy_df$x_toplot <- micro_tmp$anatomy_df$y, ifelse(slices_axes[a] == "y", micro_tmp$anatomy_df$x_toplot <- micro_tmp$anatomy_df$x, ifelse(slices_axes[a] == "z", micro_tmp$anatomy_df$x_toplot <- micro_tmp$anatomy_df$x, NA)))
            ifelse(slices_axes[a] == "x", micro_tmp$anatomy_df$y_toplot <- micro_tmp$anatomy_df$z, ifelse(slices_axes[a] == "y", micro_tmp$anatomy_df$y_toplot <- micro_tmp$anatomy_df$z, ifelse(slices_axes[a] == "z", micro_tmp$anatomy_df$y_toplot <- micro_tmp$anatomy_df$y, NA)))
            micro_tmp$anatomy_df$slice_index = micro_tmp$anatomy_df$slice_index + n_slices

            micro_tmp$anatomy_df$Micro = names[n]
            micro$anatomy_df = rbind(micro$anatomy_df, micro_tmp$anatomy_df)
        }
        micro$anatomy_df$Micro = as.factor(micro$anatomy_df$Micro)
        micro_final = rbind(micro_final, micro$anatomy_df)

        label_ukb = prepare_label_contours(paste0("../../../UKB/temporary_template/bison_ukbb/RF0_ukbb_bison_Label_2mm.mnc"), mask_path, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]], labels = seq(3,9))
        ifelse(slices_axes[a] == "x", label_ukb$contours_df$x_toplot <- label_ukb$contours_df$y, ifelse(slices_axes[a] == "y", label_ukb$contours_df$x_toplot <- label_ukb$contours_df$x, ifelse(slices_axes[a] == "z", label_ukb$contours_df$x_toplot <- label_ukb$contours_df$x, NA)))
        ifelse(slices_axes[a] == "x", label_ukb$contours_df$y_toplot <- label_ukb$contours_df$z, ifelse(slices_axes[a] == "y", label_ukb$contours_df$y_toplot <- label_ukb$contours_df$z, ifelse(slices_axes[a] == "z", label_ukb$contours_df$y_toplot <- label_ukb$contours_df$y, NA)))
        label_ukb$contours_df$slice_index = label_ukb$contours_df$slice_index + n_slices
        label_ukb_final = rbind(label_ukb_final, label_ukb$contours_df)

        # Keep track of total number of slices
        n_slices = n_slices + length(slices[[a]])
    }

    micro_final$Micro = factor(micro_final$Micro, levels=names)

    # Make invidiual plots by image types

    plot_zscore = function(df, cont_df) {
        zscores = ggplot(data=df, aes(x = x_toplot, y = y_toplot)) +
            geom_raster(aes(fill = intensity), alpha = 1) +
            scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, limits=c(-up_thresh,up_thresh), breaks=c(-up_thresh,0,up_thresh), oob=squish, name="Z-score") +
            geom_path(aes(group=interaction(label, obj)), data=cont_df, size=0.1, color="black") +
            coord_fixed(ratio = 1) +
            facet_grid(slice_index~Micro) +
            labs(title = "Zscores") +
            theme_void() +
            theme(plot.background = element_rect(fill = "white"),
                plot.title = element_text(hjust = 0.5),
                strip.text.y=element_text(size=0),
                strip.text.x=element_text(size=12),
                panel.spacing = unit(0, "lines"),
                legend.position="none")
        # ggsave("test_zscore.png")

        return(zscores)
    }

    zscores = plot_zscore(micro_final %>% filter(mask_value > 0.5), label_ukb_final)

    png(paste0(output, "_max",up_thresh,".png"), width=2500, height=1500, res=300)
    print(zscores)
    dev.off()
}

for (d in 1:ncol(dx)) {
    if (sum(dx[,d]>0, na.rm=TRUE) >= 30) {
        print(colnames(dx)[d])
        input = paste0("./results/",colnames(dx)[d],"_avg")
        output = paste0("./visualization/",colnames(dx)[d],"_avg")

        zscore_png(input, output, 1)
        zscore_png(input, output, 2)
        zscore_png(input, output, 3)
    }
}

