
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
# install.packages("magick")
# library(magick)
# devtools::install_github("thomasp85/transformr", lib="/home/m/mchakrav/parent41/R/x86_64-pc-linux-gnu-library/4.1")
# library(transformr)
# devtools::install_github("thomasp85/gganimate", lib="/home/m/mchakrav/parent41/R/x86_64-pc-linux-gnu-library/4.1")
# writeLines(output, "install_gganimate_output.txt")

names = c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")
names_long = c("dti_MD", "NODDI_ISOVF", "dti_FA", "NODDI_ICVF", "NODDI_OD", "T2star", "QSM")

# Load data

dx = as.data.frame(fread("../../../UKB/QC/exclusion_lists/dx_single.csv"))
colnames(dx) = c("Stroke", "TIA", "Subdural_H", "Subarachnoid_H", "Head_trauma", "Psychiatric", "Infection_nerv", "Abscess", 
                "Encephalitis", "Meningitis", "Chronic_neurol", "Motor_neuron_disease", "MS", "PD", "Cog_imp", "Epilepsy",
                "Head_injury", "Alcoholism", "Opioid_dep", "Other_dep", "Other_neurol", "Haemorrhage", "Ischaemic_stroke", "Fracture")

dx <- dx[, colSums(!is.na(dx)) > 0]

anatVol = mincArray(mincGetVolume("../../../UKB/temporary_template/avg.020_2mm.mnc"))
mask = mincGetVolume("../../../UKB/temporary_template/Mask_2mm.mnc")
mask_path = "../../../UKB/temporary_template/Mask_2mm.mnc"

# Source yohan's ggplot for MRI extension
source("../../yohan_ggplot/plotting_functions.R")
source("../../yohan_ggplot/plotting_functions_labels.R")

# Make png files for each subject grouped by dx

zscore_png = function(input, output, up_thresh) {
    input_id = sub(".*/(\\d{7}).*", "\\1", input)
    # print(input_id)

    # Load data for slices across different axes
    slices = list()
    slices[[1]] = seq(-70, 90, by=2) # Slices for z-axis
    slices[[2]] = seq(-80, 80, by=2) # Slices for x-axis (sagital)
    slices[[3]] = seq(-90, 70, by=2) # Slices for y-axis (coronal)

    slices_axes = c("z", "x", "y")

    flair_final = tibble(matrix(nrow = 0, ncol = 11))
    contour_final = tibble(matrix(nrow = 0, ncol = 11))
    micro_final = tibble(matrix(nrow = 0, ncol = 12))
    label_final = tibble(matrix(nrow = 0, ncol = 11))
    label_ukb_final = tibble(matrix(nrow = 0, ncol = 11))
    
    n_slices = 0
    for (a in 1:length(slices)) {
        print(slices_axes[a])
        slices_num = seq(n_slices + 1, n_slices + length(slices[[a]]))

        # Load image data
        flair = prepare_masked_anatomy(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-",input_id,"_ses-2_flair_UKB.mnc"), mask_path, slice_axis = slices_axes[a], slices[[a]])
        ifelse(slices_axes[a] == "x", flair$anatomy_df$x_toplot <- flair$anatomy_df$y, ifelse(slices_axes[a] == "y", flair$anatomy_df$x_toplot <- flair$anatomy_df$x, ifelse(slices_axes[a] == "z", flair$anatomy_df$x_toplot <- flair$anatomy_df$x, NA)))
        ifelse(slices_axes[a] == "x", flair$anatomy_df$y_toplot <- flair$anatomy_df$z, ifelse(slices_axes[a] == "y", flair$anatomy_df$y_toplot <- flair$anatomy_df$z, ifelse(slices_axes[a] == "z", flair$anatomy_df$y_toplot <- flair$anatomy_df$y, NA)))
        # flair$anatomy_df$slice_index = flair$anatomy_df$slice_index + n_slices
        flair_final = rbind(flair_final, flair$anatomy_df)

        contour = prepare_label_contours(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-", input_id,"_ses-2_Label_UKB.mnc"), mask_path, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]], labels = seq(3,9))
        ifelse(slices_axes[a] == "x", contour$contours_df$x_toplot <- contour$contours_df$y, ifelse(slices_axes[a] == "y", contour$contours_df$x_toplot <- contour$contours_df$x, ifelse(slices_axes[a] == "z", contour$contours_df$x_toplot <- contour$contours_df$x, NA)))
        ifelse(slices_axes[a] == "x", contour$contours_df$y_toplot <- contour$contours_df$z, ifelse(slices_axes[a] == "y", contour$contours_df$y_toplot <- contour$contours_df$z, ifelse(slices_axes[a] == "z", contour$contours_df$y_toplot <- contour$contours_df$y, NA)))
        # contour$contours_df$slice_index = contour$contours_df$slice_index + n_slices
        contour_final = rbind(contour_final, contour$contours_df)

        micro = prepare_masked_anatomy(paste0(input, "_",names[1],".mnc"), mask_path, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]])
        ifelse(slices_axes[a] == "x", micro$anatomy_df$x_toplot <- micro$anatomy_df$y, ifelse(slices_axes[a] == "y", micro$anatomy_df$x_toplot <- micro$anatomy_df$x, ifelse(slices_axes[a] == "z", micro$anatomy_df$x_toplot <- micro$anatomy_df$x, NA)))
        ifelse(slices_axes[a] == "x", micro$anatomy_df$y_toplot <- micro$anatomy_df$z, ifelse(slices_axes[a] == "y", micro$anatomy_df$y_toplot <- micro$anatomy_df$z, ifelse(slices_axes[a] == "z", micro$anatomy_df$y_toplot <- micro$anatomy_df$y, NA)))
        # micro$anatomy_df$slice_index = micro$anatomy_df$slice_index + n_slices

        micro$anatomy_df$Micro = names[1]
        for (n in 2:length(names)) {
            micro_tmp = prepare_masked_anatomy(paste0(input, "_",names[n],".mnc"), mask_path, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]])
            ifelse(slices_axes[a] == "x", micro_tmp$anatomy_df$x_toplot <- micro_tmp$anatomy_df$y, ifelse(slices_axes[a] == "y", micro_tmp$anatomy_df$x_toplot <- micro_tmp$anatomy_df$x, ifelse(slices_axes[a] == "z", micro_tmp$anatomy_df$x_toplot <- micro_tmp$anatomy_df$x, NA)))
            ifelse(slices_axes[a] == "x", micro_tmp$anatomy_df$y_toplot <- micro_tmp$anatomy_df$z, ifelse(slices_axes[a] == "y", micro_tmp$anatomy_df$y_toplot <- micro_tmp$anatomy_df$z, ifelse(slices_axes[a] == "z", micro_tmp$anatomy_df$y_toplot <- micro_tmp$anatomy_df$y, NA)))
            # micro_tmp$anatomy_df$slice_index = micro_tmp$anatomy_df$slice_index + n_slices

            micro_tmp$anatomy_df$Micro = names[n]
            micro$anatomy_df = rbind(micro$anatomy_df, micro_tmp$anatomy_df)
        }
        micro$anatomy_df$Micro = as.factor(micro$anatomy_df$Micro)
        micro_final = rbind(micro_final, micro$anatomy_df)

        label = prepare_masked_anatomy(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-", input_id,"_ses-2_Label_UKB.mnc"), mask_path, slice_axis = slices_axes[a], slices[[a]])
        ifelse(slices_axes[a] == "x", label$anatomy_df$x_toplot <- label$anatomy_df$y, ifelse(slices_axes[a] == "y", label$anatomy_df$x_toplot <- label$anatomy_df$x, ifelse(slices_axes[a] == "z", label$anatomy_df$x_toplot <- label$anatomy_df$x, NA)))
        ifelse(slices_axes[a] == "x", label$anatomy_df$y_toplot <- label$anatomy_df$z, ifelse(slices_axes[a] == "y", label$anatomy_df$y_toplot <- label$anatomy_df$z, ifelse(slices_axes[a] == "z", label$anatomy_df$y_toplot <- label$anatomy_df$y, NA)))
        # label$anatomy_df$slice_index = label$anatomy_df$slice_index + n_slices
        label$anatomy_df$intensity = round(label$anatomy_df$intensity,1)
        label_final = rbind(label_final, label$anatomy_df)

        label_ukb = prepare_label_contours(paste0("../../../UKB/temporary_template/bison_ukbb/RF0_ukbb_bison_Label_2mm.mnc"), mask_path, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]], labels = seq(3,9))
        ifelse(slices_axes[a] == "x", label_ukb$contours_df$x_toplot <- label_ukb$contours_df$y, ifelse(slices_axes[a] == "y", label_ukb$contours_df$x_toplot <- label_ukb$contours_df$x, ifelse(slices_axes[a] == "z", label_ukb$contours_df$x_toplot <- label_ukb$contours_df$x, NA)))
        ifelse(slices_axes[a] == "x", label_ukb$contours_df$y_toplot <- label_ukb$contours_df$z, ifelse(slices_axes[a] == "y", label_ukb$contours_df$y_toplot <- label_ukb$contours_df$z, ifelse(slices_axes[a] == "z", label_ukb$contours_df$y_toplot <- label_ukb$contours_df$y, NA)))
        # label_ukb$contours_df$slice_index = label_ukb$contours_df$slice_index + n_slices
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
            facet_grid(slice_axis~Micro) +
            labs(title = "Zscores") +
            theme_void() +
            theme(plot.background = element_rect(fill = "white"),
                plot.title = element_text(hjust = 0.5),
                strip.text.y=element_text(size=0),
                strip.text.x=element_text(size=12),
                panel.spacing = unit(0, "lines"),
                legend.position="none")

        return(zscores)
    }

    plot_flair = function(df, cont_df) {
        flair = ggplot(data = df, aes(x = x_toplot, y = y_toplot)) +
            geom_raster(aes(fill = intensity), alpha = 1) +
            scale_fill_gradient(low='black', high= 'white', limits=c(0,150), oob=squish, guide = "none") +
            geom_path(aes(group=interaction(label, obj)), data=cont_df, size=0.1, color="red") +
            coord_fixed(ratio = 1) +
            facet_wrap(slice_axis~., nrow=3) +
            labs(title = "FLAIR") +
            theme_void() +
            theme(plot.background = element_rect(fill = "white"),
                plot.title = element_text(hjust = 0.5),
                strip.text.y=element_text(size=12),
                strip.text.x=element_text(size=0),
                panel.spacing = unit(0, "lines"))
        # ggsave("test_flair.png")
    }

    plot_bison = function(df) {
        names_bison = c('Ventricules', 'CSF', 'Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH')
        color_scale_bison = tibble(
            labels=seq(from=0, to=9), 
            colors=c("#000000", "#6a009d", "#0035dd", "#00a4bb", "#009b0f", "#00e100", "#ccf900", "#ffb000", "#e50000", "#fffefb"))
        
        bison = ggplot(data = df %>% left_join(color_scale_bison, by=c("intensity"="labels")), aes(x = x_toplot, y = y_toplot)) +
            geom_raster(aes(fill = colors), alpha = 1) +
            scale_fill_identity(guide = "none") +
            coord_fixed(ratio = 1) +
            facet_wrap(slice_axis~., nrow=3) +
            labs(title = "Labels") +
            theme_void() +
            theme(plot.background = element_rect(fill = "white"),
                plot.title = element_text(hjust = 0.5),
                strip.text.y=element_text(size=12),
                strip.text.x=element_text(size=0),
                panel.spacing = unit(0, "lines"))
        # ggsave("test_bison.png")

        return(bison)
    }

    # Make temporary pngs for each slice
    dir.create(paste0("./tmp/",input_id), showWarnings=FALSE)
    for (i in 1:length(unique(flair_final$slice_index))) {
        print(i)

        zscores = plot_zscore(micro_final %>% filter(slice_index == i), contour_final %>% filter(slice_index == i))
        flair = plot_flair(flair_final %>% filter(slice_index == i), label_ukb_final %>% filter(slice_index == i))
        bison = plot_bison(label_final %>% filter(slice_index == i))

        # Combine image types with patchwork
        plt_layout <- "
        ABCCCCCCC
        ABCCCCCCC
        "
        
        final_plot <- flair + bison + zscores + plot_layout(design=plt_layout)

        # png(paste0(output, "_zscore_max",up_thresh,".png"), width=2500, height=1500, res=300)
        png(paste0("./tmp/",input_id,"/slice_",i,".png"), width=2500, height=1000, res=300)
        print(final_plot)
        dev.off()
    }

    # Make GIF with ImageMagick
    command = paste0("ffmpeg -y -framerate 10 -i ./tmp/",input_id,"/slice_%00d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p ",output,".mp4")
    system(command)
}

raw_png = function(input, output, anlm) {

    raw_up_thresh = c(0.003, 0.5, 0.85, 0.85, 1, 75, 150)
    raw_down_thresh = c(0, 0, 0, 0, 0, 0, -100)

    input_id = str_extract(input, "(?<=sub-)\\d+")
    # print(input_id)

    # Function to rescale intensities between 0 and 1
    rescale_vector <- function(vector, up_t, down_t) {
        adjusted_vec = pmax(pmin(vector, up_t), down_t)
        new_min=0
        new_max=1
        min_value <- min(adjusted_vec)
        max_value <- max(adjusted_vec)
        scaled_vector <- (adjusted_vec - min_value) / (max_value - min_value) * (new_max - new_min) + new_min
        return(scaled_vector)
    }

    # Load data for slices across different axes
    slices = list()
    slices[[1]] = seq(-70, 90, by=2) # Slices for z-axis
    slices[[2]] = seq(-80, 80, by=2) # Slices for x-axis (sagital)
    slices[[3]] = seq(-90, 70, by=2) # Slices for y-axis (coronal)

    slices_axes = c("z", "x", "y")

    flair_final = tibble(matrix(nrow = 0, ncol = 11))
    contour_final = tibble(matrix(nrow = 0, ncol = 11))
    micro_final = tibble(matrix(nrow = 0, ncol = 12))
    label_final = tibble(matrix(nrow = 0, ncol = 11))
    label_ukb_final = tibble(matrix(nrow = 0, ncol = 11))
    
    n_slices = 0
    for (a in 1:length(slices)) {
        print(slices_axes[a])
        slices_num = seq(n_slices + 1, n_slices + length(slices[[a]]))

        # Load image data
        flair = prepare_masked_anatomy(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-",input_id,"_ses-2_flair_UKB.mnc"), mask_path, slice_axis = slices_axes[a], slices[[a]])
        ifelse(slices_axes[a] == "x", flair$anatomy_df$x_toplot <- flair$anatomy_df$y, ifelse(slices_axes[a] == "y", flair$anatomy_df$x_toplot <- flair$anatomy_df$x, ifelse(slices_axes[a] == "z", flair$anatomy_df$x_toplot <- flair$anatomy_df$x, NA)))
        ifelse(slices_axes[a] == "x", flair$anatomy_df$y_toplot <- flair$anatomy_df$z, ifelse(slices_axes[a] == "y", flair$anatomy_df$y_toplot <- flair$anatomy_df$z, ifelse(slices_axes[a] == "z", flair$anatomy_df$y_toplot <- flair$anatomy_df$y, NA)))
        # flair$anatomy_df$slice_index = flair$anatomy_df$slice_index + n_slices
        flair_final = rbind(flair_final, flair$anatomy_df)

        contour = prepare_label_contours(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-", input_id,"_ses-2_Label_UKB.mnc"), mask_path, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]], labels = seq(3,9))
        ifelse(slices_axes[a] == "x", contour$contours_df$x_toplot <- contour$contours_df$y, ifelse(slices_axes[a] == "y", contour$contours_df$x_toplot <- contour$contours_df$x, ifelse(slices_axes[a] == "z", contour$contours_df$x_toplot <- contour$contours_df$x, NA)))
        ifelse(slices_axes[a] == "x", contour$contours_df$y_toplot <- contour$contours_df$z, ifelse(slices_axes[a] == "y", contour$contours_df$y_toplot <- contour$contours_df$z, ifelse(slices_axes[a] == "z", contour$contours_df$y_toplot <- contour$contours_df$y, NA)))
        # contour$contours_df$slice_index = contour$contours_df$slice_index + n_slices
        contour_final = rbind(contour_final, contour$contours_df)

        # micro = prepare_masked_anatomy(list.files(input_dir, pattern = paste0("sub-", input_id, "_ses-2_",names_long[1],"*"), full.names = TRUE), mask_path, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]])
        micro = prepare_masked_anatomy(paste0(input, "_",names_long[1],"_UKB",anlm,".mnc"), mask_path, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]])
        ifelse(slices_axes[a] == "x", micro$anatomy_df$x_toplot <- micro$anatomy_df$y, ifelse(slices_axes[a] == "y", micro$anatomy_df$x_toplot <- micro$anatomy_df$x, ifelse(slices_axes[a] == "z", micro$anatomy_df$x_toplot <- micro$anatomy_df$x, NA)))
        ifelse(slices_axes[a] == "x", micro$anatomy_df$y_toplot <- micro$anatomy_df$z, ifelse(slices_axes[a] == "y", micro$anatomy_df$y_toplot <- micro$anatomy_df$z, ifelse(slices_axes[a] == "z", micro$anatomy_df$y_toplot <- micro$anatomy_df$y, NA)))
        # micro$anatomy_df$slice_index = micro$anatomy_df$slice_index + n_slices
        micro$anatomy_df$intensity = rescale_vector(micro$anatomy_df$intensity, raw_up_thresh[1], raw_down_thresh[1])

        micro$anatomy_df$Micro = names[1]
        for (n in 2:length(names)) {
            micro_tmp = prepare_masked_anatomy(paste0(input, "_",names_long[n],"_UKB",anlm,".mnc"), mask_path, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]])
            ifelse(slices_axes[a] == "x", micro_tmp$anatomy_df$x_toplot <- micro_tmp$anatomy_df$y, ifelse(slices_axes[a] == "y", micro_tmp$anatomy_df$x_toplot <- micro_tmp$anatomy_df$x, ifelse(slices_axes[a] == "z", micro_tmp$anatomy_df$x_toplot <- micro_tmp$anatomy_df$x, NA)))
            ifelse(slices_axes[a] == "x", micro_tmp$anatomy_df$y_toplot <- micro_tmp$anatomy_df$z, ifelse(slices_axes[a] == "y", micro_tmp$anatomy_df$y_toplot <- micro_tmp$anatomy_df$z, ifelse(slices_axes[a] == "z", micro_tmp$anatomy_df$y_toplot <- micro_tmp$anatomy_df$y, NA)))
            # micro_tmp$anatomy_df$slice_index = micro_tmp$anatomy_df$slice_index + n_slices
            micro_tmp$anatomy_df$intensity = rescale_vector(micro_tmp$anatomy_df$intensity, raw_up_thresh[n], raw_down_thresh[n])

            micro_tmp$anatomy_df$Micro = names[n]
            micro$anatomy_df = rbind(micro$anatomy_df, micro_tmp$anatomy_df)
        }
        micro$anatomy_df$Micro = as.factor(micro$anatomy_df$Micro)
        micro_final = rbind(micro_final, micro$anatomy_df)

        label = prepare_masked_anatomy(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-", input_id,"_ses-2_Label_UKB.mnc"), mask_path, slice_axis = slices_axes[a], slices[[a]])
        ifelse(slices_axes[a] == "x", label$anatomy_df$x_toplot <- label$anatomy_df$y, ifelse(slices_axes[a] == "y", label$anatomy_df$x_toplot <- label$anatomy_df$x, ifelse(slices_axes[a] == "z", label$anatomy_df$x_toplot <- label$anatomy_df$x, NA)))
        ifelse(slices_axes[a] == "x", label$anatomy_df$y_toplot <- label$anatomy_df$z, ifelse(slices_axes[a] == "y", label$anatomy_df$y_toplot <- label$anatomy_df$z, ifelse(slices_axes[a] == "z", label$anatomy_df$y_toplot <- label$anatomy_df$y, NA)))
        # label$anatomy_df$slice_index = label$anatomy_df$slice_index + n_slices
        label$anatomy_df$intensity = round(label$anatomy_df$intensity,1)
        label_final = rbind(label_final, label$anatomy_df)

        label_ukb = prepare_label_contours(paste0("../../../UKB/temporary_template/bison_ukbb/RF0_ukbb_bison_Label_2mm.mnc"), mask_path, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]], labels = seq(3,9))
        ifelse(slices_axes[a] == "x", label_ukb$contours_df$x_toplot <- label_ukb$contours_df$y, ifelse(slices_axes[a] == "y", label_ukb$contours_df$x_toplot <- label_ukb$contours_df$x, ifelse(slices_axes[a] == "z", label_ukb$contours_df$x_toplot <- label_ukb$contours_df$x, NA)))
        ifelse(slices_axes[a] == "x", label_ukb$contours_df$y_toplot <- label_ukb$contours_df$z, ifelse(slices_axes[a] == "y", label_ukb$contours_df$y_toplot <- label_ukb$contours_df$z, ifelse(slices_axes[a] == "z", label_ukb$contours_df$y_toplot <- label_ukb$contours_df$y, NA)))
        # label_ukb$contours_df$slice_index = label_ukb$contours_df$slice_index + n_slices
        label_ukb_final = rbind(label_ukb_final, label_ukb$contours_df)

        # Keep track of total number of slices
        n_slices = n_slices + length(slices[[a]])
    }

    micro_final$Micro = factor(micro_final$Micro, levels=names)

    # Make invidiual plots by image types

    plot_micro = function(df, cont_df) {
        zscores = ggplot(data=df, aes(x = x_toplot, y = y_toplot)) +
            geom_raster(aes(fill = intensity), alpha = 1) +
            scale_fill_gradient(low="white", high="red", limits=c(0,1), breaks=c(0,1), oob=squish, name="Micro") +
            geom_path(aes(group=interaction(label, obj)), data=cont_df, size=0.1, color="black") +
            coord_fixed(ratio = 1) +
            facet_grid(slice_axis~Micro) +
            labs(title = "Zscores") +
            theme_void() +
            theme(plot.background = element_rect(fill = "white"),
                plot.title = element_text(hjust = 0.5),
                strip.text.y=element_text(size=0),
                strip.text.x=element_text(size=12),
                panel.spacing = unit(0, "lines"),
                legend.position="none")

        return(zscores)
    }

    plot_flair = function(df, cont_df) {
        flair = ggplot(data = df, aes(x = x_toplot, y = y_toplot)) +
            geom_raster(aes(fill = intensity), alpha = 1) +
            scale_fill_gradient(low='black', high= 'white', limits=c(0,150), oob=squish, guide = "none") +
            geom_path(aes(group=interaction(label, obj)), data=cont_df, size=0.1, color="red") +
            coord_fixed(ratio = 1) +
            facet_wrap(slice_axis~., nrow=3) +
            labs(title = "FLAIR") +
            theme_void() +
            theme(plot.background = element_rect(fill = "white"),
                plot.title = element_text(hjust = 0.5),
                strip.text.y=element_text(size=12),
                strip.text.x=element_text(size=0),
                panel.spacing = unit(0, "lines"))
        # ggsave("test_flair.png")
    }

    plot_bison = function(df) {
        names_bison = c('Ventricules', 'CSF', 'Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH')
        color_scale_bison = tibble(
            labels=seq(from=0, to=9), 
            colors=c("#000000", "#6a009d", "#0035dd", "#00a4bb", "#009b0f", "#00e100", "#ccf900", "#ffb000", "#e50000", "#fffefb"))
        
        bison = ggplot(data = df %>% left_join(color_scale_bison, by=c("intensity"="labels")), aes(x = x_toplot, y = y_toplot)) +
            geom_raster(aes(fill = colors), alpha = 1) +
            scale_fill_identity(guide = "none") +
            coord_fixed(ratio = 1) +
            facet_wrap(slice_axis~., nrow=3) +
            labs(title = "Labels") +
            theme_void() +
            theme(plot.background = element_rect(fill = "white"),
                plot.title = element_text(hjust = 0.5),
                strip.text.y=element_text(size=12),
                strip.text.x=element_text(size=0),
                panel.spacing = unit(0, "lines"))
        # ggsave("test_bison.png")

        return(bison)
    }

    # Make temporary pngs for each slice
    dir.create(paste0("./tmp/",input_id), showWarnings=FALSE)
    for (i in 1:length(unique(flair_final$slice_index))) {
        print(i)

        micros = plot_micro(micro_final %>% filter(slice_index == i), contour_final %>% filter(slice_index == i))
        flair = plot_flair(flair_final %>% filter(slice_index == i), label_ukb_final %>% filter(slice_index == i))
        bison = plot_bison(label_final %>% filter(slice_index == i))

        # Combine image types with patchwork
        plt_layout <- "
        ABCCCCCCC
        ABCCCCCCC
        "
        
        final_plot <- flair + bison + micros + plot_layout(design=plt_layout)

        # png(paste0(output, "_zscore_max",up_thresh,".png"), width=2500, height=1500, res=300)
        png(paste0("./tmp/",input_id,"/slice_raw_",i,".png"), width=2500, height=1000, res=300)
        print(final_plot)
        dev.off()
    }

    # Make GIF with ImageMagick
    command = paste0("ffmpeg -y -framerate 10 -i ./tmp/",input_id,"/slice_raw_%00d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p ",output,"_raw.mp4")
    system(command)
}

# Run

cores = detectCores()
parallelCluster <- makeCluster(cores-5, type = "SOCK", methods = FALSE)
setDefaultCluster(parallelCluster)
registerDoParallel(parallelCluster)

up_thresholds <- c(3)

id_list = c("4519919", "5836644", "1421334")
viz_dir = paste0("./visualization/mp4")

# for (i in 1:length(id_list)) {
foreach(i=1:length(id_list), .packages = c('RMINC', 'RColorBrewer', 'tidyverse', 'scales', 'patchwork', 'ggnewscale')) %dopar% {

    print(id_list[i])

    dx_id = na.omit(colnames(dx)[apply(dx == id_list[i], 2, any)])

    res_dir = paste0("./results/",dx_id[1], "/", id_list[i])
    dir.create(viz_dir, showWarnings=FALSE)
    
    # Make PNG (raw)
    input_anlm = paste0(res_dir, "/", id_list[i],"_anlm")
    output_anlm = paste0(viz_dir, "/", id_list[i], "_",dx_id[1],"_anlm")
    zscore_png(input_anlm, output_anlm, 3)

    input_anlm = paste0("../../maps_UKB_space_anlm_all/sub-",id_list[i],"_ses-2")
    output_anlm = paste0(viz_dir, "/", id_list[i], "_",dx_id[1], "_anlm")
    raw_png(input_anlm, output_anlm, "_anlm")

}

