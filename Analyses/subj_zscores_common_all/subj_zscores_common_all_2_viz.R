
library(data.table)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(RMINC)
library(viridis)
library(RColorBrewer)
library(foreach)
library(doParallel)
library(ggnewscale)
library(scales)

args = commandArgs(trailingOnly=TRUE)

dir.create("./visualization/indiv", showWarnings=FALSE)

# args = c()
# args[1] = 32

num_chunk = as.numeric(args[1])

names = c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM", "jacobians_abs", "jacobians_rel")
names_long = c("dti_MD_UKB", "NODDI_ISOVF_UKB", "dti_FA_UKB", "NODDI_ICVF_UKB", "NODDI_OD_UKB", "T2star_UKB", "QSM_UKB", "jacobians_abs_2mm", "jacobians_rel_2mm")
names_plot = c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM", "J(abs)", "J(rel)")

# Load data

inclusions = as.data.frame(fread(paste0("./tmp/ids_micro_c",num_chunk,"_FA_anlm.txt")))
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

rm(demo)

# Source yohan's ggplot for MRI extension
source("../../yohan_ggplot/plotting_functions.R")
source("../../yohan_ggplot/plotting_functions_labels.R")

# Make png files for each subject grouped by dx

zscore_png = function(input, output, up_thresh) {
    input_id = sub(".*/(\\d{7}).*", "\\1", input)
    # print(input_id)

    # Load data for slices across different axes
    slices = list()
    slices[[1]] = c(0,20) # Slices for z-axis
    slices[[2]] = c(-5, -30) # Slices for x-axis
    slices[[3]] = c(-20) # Slices for y-axis

    slices_axes = c("z", "x", "y")

    n_slices = 0

    flair_final = tibble(matrix(nrow = 0, ncol = 11))
    contour_final = tibble(matrix(nrow = 0, ncol = 11))
    micro_final = tibble(matrix(nrow = 0, ncol = 12))
    label_final = tibble(matrix(nrow = 0, ncol = 11))
    label_ukb_final = tibble(matrix(nrow = 0, ncol = 11))

    for (a in 1:length(slices)) {
        # print(slices_axes[a])
        slices_num = seq(n_slices + 1, n_slices + length(slices[[a]]))

        # Load image data
        flair = prepare_masked_anatomy(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-",input_id,"_ses-2_flair_UKB.mnc"), mask_path, slice_axis = slices_axes[a], slices[[a]])
        ifelse(slices_axes[a] == "x", flair$anatomy_df$x_toplot <- flair$anatomy_df$y, ifelse(slices_axes[a] == "y", flair$anatomy_df$x_toplot <- flair$anatomy_df$x, ifelse(slices_axes[a] == "z", flair$anatomy_df$x_toplot <- flair$anatomy_df$x, NA)))
        ifelse(slices_axes[a] == "x", flair$anatomy_df$y_toplot <- flair$anatomy_df$z, ifelse(slices_axes[a] == "y", flair$anatomy_df$y_toplot <- flair$anatomy_df$z, ifelse(slices_axes[a] == "z", flair$anatomy_df$y_toplot <- flair$anatomy_df$y, NA)))
        flair$anatomy_df$slice_index = flair$anatomy_df$slice_index + n_slices
        flair_final = rbind(flair_final, flair$anatomy_df)

        contour = prepare_label_contours(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-", input_id,"_ses-2_Label_UKB.mnc"), mask_path, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]], labels = seq(3,9))
        ifelse(slices_axes[a] == "x", contour$contours_df$x_toplot <- contour$contours_df$y, ifelse(slices_axes[a] == "y", contour$contours_df$x_toplot <- contour$contours_df$x, ifelse(slices_axes[a] == "z", contour$contours_df$x_toplot <- contour$contours_df$x, NA)))
        ifelse(slices_axes[a] == "x", contour$contours_df$y_toplot <- contour$contours_df$z, ifelse(slices_axes[a] == "y", contour$contours_df$y_toplot <- contour$contours_df$z, ifelse(slices_axes[a] == "z", contour$contours_df$y_toplot <- contour$contours_df$y, NA)))
        contour$contours_df$slice_index = contour$contours_df$slice_index + n_slices
        contour_final = rbind(contour_final, contour$contours_df)

        micro = prepare_masked_anatomy(paste0(input, "_",names[1],".mnc"), mask_path, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]])
        ifelse(slices_axes[a] == "x", micro$anatomy_df$x_toplot <- micro$anatomy_df$y, ifelse(slices_axes[a] == "y", micro$anatomy_df$x_toplot <- micro$anatomy_df$x, ifelse(slices_axes[a] == "z", micro$anatomy_df$x_toplot <- micro$anatomy_df$x, NA)))
        ifelse(slices_axes[a] == "x", micro$anatomy_df$y_toplot <- micro$anatomy_df$z, ifelse(slices_axes[a] == "y", micro$anatomy_df$y_toplot <- micro$anatomy_df$z, ifelse(slices_axes[a] == "z", micro$anatomy_df$y_toplot <- micro$anatomy_df$y, NA)))
        micro$anatomy_df$slice_index = micro$anatomy_df$slice_index + n_slices

        micro$anatomy_df$Micro = names[1]
        for (n in 2:length(names)) {
            micro_tmp = prepare_masked_anatomy(paste0(input, "_",names[n],".mnc"), mask_path, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]])
            ifelse(slices_axes[a] == "x", micro_tmp$anatomy_df$x_toplot <- micro_tmp$anatomy_df$y, ifelse(slices_axes[a] == "y", micro_tmp$anatomy_df$x_toplot <- micro_tmp$anatomy_df$x, ifelse(slices_axes[a] == "z", micro_tmp$anatomy_df$x_toplot <- micro_tmp$anatomy_df$x, NA)))
            ifelse(slices_axes[a] == "x", micro_tmp$anatomy_df$y_toplot <- micro_tmp$anatomy_df$z, ifelse(slices_axes[a] == "y", micro_tmp$anatomy_df$y_toplot <- micro_tmp$anatomy_df$z, ifelse(slices_axes[a] == "z", micro_tmp$anatomy_df$y_toplot <- micro_tmp$anatomy_df$y, NA)))
            micro_tmp$anatomy_df$slice_index = micro_tmp$anatomy_df$slice_index + n_slices

            micro_tmp$anatomy_df$Micro = names[n]
            micro$anatomy_df = rbind(micro$anatomy_df, micro_tmp$anatomy_df)
        }
        micro$anatomy_df$Micro = as.factor(micro$anatomy_df$Micro)
        micro_final = rbind(micro_final, micro$anatomy_df)

        label = prepare_masked_anatomy(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-", input_id,"_ses-2_Label_UKB.mnc"), mask_path, slice_axis = slices_axes[a], slices[[a]])
        ifelse(slices_axes[a] == "x", label$anatomy_df$x_toplot <- label$anatomy_df$y, ifelse(slices_axes[a] == "y", label$anatomy_df$x_toplot <- label$anatomy_df$x, ifelse(slices_axes[a] == "z", label$anatomy_df$x_toplot <- label$anatomy_df$x, NA)))
        ifelse(slices_axes[a] == "x", label$anatomy_df$y_toplot <- label$anatomy_df$z, ifelse(slices_axes[a] == "y", label$anatomy_df$y_toplot <- label$anatomy_df$z, ifelse(slices_axes[a] == "z", label$anatomy_df$y_toplot <- label$anatomy_df$y, NA)))
        label$anatomy_df$slice_index = label$anatomy_df$slice_index + n_slices
        label$anatomy_df$intensity = round(label$anatomy_df$intensity,1)
        label_final = rbind(label_final, label$anatomy_df)

        # label_ukb = prepare_label_contours(paste0("../../../UKB/temporary_template/bison_ukbb/RF0_ukbb_bison_Label_2mm.mnc"), mask_path, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]], labels = seq(3,9))
        # ifelse(slices_axes[a] == "x", label_ukb$contours_df$x_toplot <- label_ukb$contours_df$y, ifelse(slices_axes[a] == "y", label_ukb$contours_df$x_toplot <- label_ukb$contours_df$x, ifelse(slices_axes[a] == "z", label_ukb$contours_df$x_toplot <- label_ukb$contours_df$x, NA)))
        # ifelse(slices_axes[a] == "x", label_ukb$contours_df$y_toplot <- label_ukb$contours_df$z, ifelse(slices_axes[a] == "y", label_ukb$contours_df$y_toplot <- label_ukb$contours_df$z, ifelse(slices_axes[a] == "z", label_ukb$contours_df$y_toplot <- label_ukb$contours_df$y, NA)))
        # label_ukb$contours_df$slice_index = label_ukb$contours_df$slice_index + n_slices
        # label_ukb_final = rbind(label_ukb_final, label_ukb$contours_df)

        # Keep track of total number of slices
        n_slices = n_slices + length(slices[[a]])
    }

    micro_final$Micro = factor(micro_final$Micro, levels=names)
    levels(micro_final$Micro) = names_plot

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

    # plot_flair = function(df, cont_df) {
    plot_flair = function(df) {
        flair = ggplot(data = df, aes(x = x_toplot, y = y_toplot)) +
            geom_raster(aes(fill = intensity), alpha = 1) +
            scale_fill_gradient(low='black', high= 'white', limits=c(0,150), oob=squish, guide = "none") +
            # geom_path(aes(group=interaction(label, obj)), data=cont_df, size=0.1, color="red") +
            coord_fixed(ratio = 1) +
            facet_wrap(~slice_index, ncol=1) +
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
            facet_wrap(~slice_index, ncol=1) +
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

    zscores = plot_zscore(micro_final %>% filter(mask_value > 0.5), contour_final)
    # flair = plot_flair(flair_final %>% filter(mask_value > 0.5), label_ukb_final)
    flair = plot_flair(flair_final %>% filter(mask_value > 0.5))
    bison = plot_bison(label_final %>% filter(mask_value > 0.5))

    # Combine image types with patchwork
    plt_layout <- "
    ABCCCCCCCCC
    ABCCCCCCCCC
    "
    
    final_plot <- flair + bison + zscores + plot_layout(design=plt_layout)

    png(paste0(output, "_zscore_max",up_thresh,".png"), width=3000, height=1500, res=300)
    print(final_plot)
    dev.off()
    print(paste0(output, "_zscore_max",up_thresh,".png"))
}

raw_png = function(input, output, anlm) {

    raw_up_thresh = c(0.003, 0.5, 0.85, 0.85, 1, 75, 150, 1, 1)
    raw_down_thresh = c(0, 0, 0, 0, 0, 0, -100, -1, -1)

    input_id = str_extract(input, "(?<=sub-)\\d+")
    mask_id = paste0("../../masks_tissue/sub-",input_id,"_ses-2_mask_tissue.mnc")
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
    slices[[1]] = c(0,20) # Slices for z-axis
    slices[[2]] = c(-5, -30) # Slices for x-axis
    slices[[3]] = c(-20) # Slices for y-axis

    slices_axes = c("z", "x", "y")

    n_slices = 0

    flair_final = tibble(matrix(nrow = 0, ncol = 11))
    contour_final = tibble(matrix(nrow = 0, ncol = 11))
    micro_final = tibble(matrix(nrow = 0, ncol = 12))
    label_final = tibble(matrix(nrow = 0, ncol = 11))
    label_ukb_final = tibble(matrix(nrow = 0, ncol = 11))

    for (a in 1:length(slices)) {
        # print(slices_axes[a])
        slices_num = seq(n_slices + 1, n_slices + length(slices[[a]]))

        # Load image data
        flair = prepare_masked_anatomy(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-",input_id,"_ses-2_flair_UKB.mnc"), mask_path, slice_axis = slices_axes[a], slices[[a]])
        ifelse(slices_axes[a] == "x", flair$anatomy_df$x_toplot <- flair$anatomy_df$y, ifelse(slices_axes[a] == "y", flair$anatomy_df$x_toplot <- flair$anatomy_df$x, ifelse(slices_axes[a] == "z", flair$anatomy_df$x_toplot <- flair$anatomy_df$x, NA)))
        ifelse(slices_axes[a] == "x", flair$anatomy_df$y_toplot <- flair$anatomy_df$z, ifelse(slices_axes[a] == "y", flair$anatomy_df$y_toplot <- flair$anatomy_df$z, ifelse(slices_axes[a] == "z", flair$anatomy_df$y_toplot <- flair$anatomy_df$y, NA)))
        flair$anatomy_df$slice_index = flair$anatomy_df$slice_index + n_slices
        flair_final = rbind(flair_final, flair$anatomy_df)

        contour = prepare_label_contours(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-", input_id,"_ses-2_Label_UKB.mnc"), mask_path, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]], labels = seq(3,9))
        ifelse(slices_axes[a] == "x", contour$contours_df$x_toplot <- contour$contours_df$y, ifelse(slices_axes[a] == "y", contour$contours_df$x_toplot <- contour$contours_df$x, ifelse(slices_axes[a] == "z", contour$contours_df$x_toplot <- contour$contours_df$x, NA)))
        ifelse(slices_axes[a] == "x", contour$contours_df$y_toplot <- contour$contours_df$z, ifelse(slices_axes[a] == "y", contour$contours_df$y_toplot <- contour$contours_df$z, ifelse(slices_axes[a] == "z", contour$contours_df$y_toplot <- contour$contours_df$y, NA)))
        contour$contours_df$slice_index = contour$contours_df$slice_index + n_slices
        contour_final = rbind(contour_final, contour$contours_df)

        # micro = prepare_masked_anatomy(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-",input_id,"_ses-2_",names_long[1],"_UKB.mnc"), mask_path, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]])
        micro = prepare_masked_anatomy(paste0(input, "_",names_long[1],anlm,".mnc"), mask_id, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]])
        ifelse(slices_axes[a] == "x", micro$anatomy_df$x_toplot <- micro$anatomy_df$y, ifelse(slices_axes[a] == "y", micro$anatomy_df$x_toplot <- micro$anatomy_df$x, ifelse(slices_axes[a] == "z", micro$anatomy_df$x_toplot <- micro$anatomy_df$x, NA)))
        ifelse(slices_axes[a] == "x", micro$anatomy_df$y_toplot <- micro$anatomy_df$z, ifelse(slices_axes[a] == "y", micro$anatomy_df$y_toplot <- micro$anatomy_df$z, ifelse(slices_axes[a] == "z", micro$anatomy_df$y_toplot <- micro$anatomy_df$y, NA)))
        micro$anatomy_df$slice_index = micro$anatomy_df$slice_index + n_slices
        micro$anatomy_df$intensity = rescale_vector(micro$anatomy_df$intensity, raw_up_thresh[1], raw_down_thresh[1])

        micro$anatomy_df$Micro = names[1]
        for (n in 2:length(names)) {
            # micro_tmp = prepare_masked_anatomy(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-",input_id,"_ses-2_",names_long[n],"_UKB.mnc"), mask_path, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]])
            micro_tmp = prepare_masked_anatomy(paste0(input, "_",names_long[n],anlm,".mnc"), mask_id, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]])
            ifelse(slices_axes[a] == "x", micro_tmp$anatomy_df$x_toplot <- micro_tmp$anatomy_df$y, ifelse(slices_axes[a] == "y", micro_tmp$anatomy_df$x_toplot <- micro_tmp$anatomy_df$x, ifelse(slices_axes[a] == "z", micro_tmp$anatomy_df$x_toplot <- micro_tmp$anatomy_df$x, NA)))
            ifelse(slices_axes[a] == "x", micro_tmp$anatomy_df$y_toplot <- micro_tmp$anatomy_df$z, ifelse(slices_axes[a] == "y", micro_tmp$anatomy_df$y_toplot <- micro_tmp$anatomy_df$z, ifelse(slices_axes[a] == "z", micro_tmp$anatomy_df$y_toplot <- micro_tmp$anatomy_df$y, NA)))
            micro_tmp$anatomy_df$slice_index = micro_tmp$anatomy_df$slice_index + n_slices
            micro_tmp$anatomy_df$intensity = rescale_vector(micro_tmp$anatomy_df$intensity, raw_up_thresh[n], raw_down_thresh[n])

            micro_tmp$anatomy_df$Micro = names[n]
            micro$anatomy_df = rbind(micro$anatomy_df, micro_tmp$anatomy_df)
        }
        micro$anatomy_df$Micro = as.factor(micro$anatomy_df$Micro)
        micro_final = rbind(micro_final, micro$anatomy_df)

        label = prepare_masked_anatomy(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-", input_id,"_ses-2_Label_UKB.mnc"), mask_path, slice_axis = slices_axes[a], slices[[a]])
        ifelse(slices_axes[a] == "x", label$anatomy_df$x_toplot <- label$anatomy_df$y, ifelse(slices_axes[a] == "y", label$anatomy_df$x_toplot <- label$anatomy_df$x, ifelse(slices_axes[a] == "z", label$anatomy_df$x_toplot <- label$anatomy_df$x, NA)))
        ifelse(slices_axes[a] == "x", label$anatomy_df$y_toplot <- label$anatomy_df$z, ifelse(slices_axes[a] == "y", label$anatomy_df$y_toplot <- label$anatomy_df$z, ifelse(slices_axes[a] == "z", label$anatomy_df$y_toplot <- label$anatomy_df$y, NA)))
        label$anatomy_df$slice_index = label$anatomy_df$slice_index + n_slices
        label$anatomy_df$intensity = round(label$anatomy_df$intensity,1)
        label_final = rbind(label_final, label$anatomy_df)

        # label_ukb = prepare_label_contours(paste0("../../../UKB/temporary_template/bison_ukbb/RF0_ukbb_bison_Label_2mm.mnc"), mask_path, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]], labels = seq(3,9))
        # ifelse(slices_axes[a] == "x", label_ukb$contours_df$x_toplot <- label_ukb$contours_df$y, ifelse(slices_axes[a] == "y", label_ukb$contours_df$x_toplot <- label_ukb$contours_df$x, ifelse(slices_axes[a] == "z", label_ukb$contours_df$x_toplot <- label_ukb$contours_df$x, NA)))
        # ifelse(slices_axes[a] == "x", label_ukb$contours_df$y_toplot <- label_ukb$contours_df$z, ifelse(slices_axes[a] == "y", label_ukb$contours_df$y_toplot <- label_ukb$contours_df$z, ifelse(slices_axes[a] == "z", label_ukb$contours_df$y_toplot <- label_ukb$contours_df$y, NA)))
        # label_ukb$contours_df$slice_index = label_ukb$contours_df$slice_index + n_slices
        # label_ukb_final = rbind(label_ukb_final, label_ukb$contours_df)

        # Keep track of total number of slices
        n_slices = n_slices + length(slices[[a]])
    }

    micro_final$Micro = factor(micro_final$Micro, levels=names)
    levels(micro_final$Micro) = names_plot

    # Make invidiual plots by image types

    plot_micro = function(df, cont_df) {
        zscores = ggplot(data=df, aes(x = x_toplot, y = y_toplot)) +
            geom_raster(aes(fill = intensity), alpha = 1) +
            # scale_fill_gradient(low="white", high="red", limits=c(0,1), breaks=c(0,1), oob=squish, name="Micro") +
            scale_fill_viridis(option="A", limits=c(0,1), breaks=c(0,1), oob=squish, name="Micro") +
            geom_path(aes(group=interaction(label, obj)), data=cont_df, size=0.1, color="black") +
            coord_fixed(ratio = 1) +
            facet_grid(slice_index~Micro) +
            labs(title = "Raw") +
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

    # plot_flair = function(df, cont_df) {
    plot_flair = function(df) {
        flair = ggplot(data = df, aes(x = x_toplot, y = y_toplot)) +
            geom_raster(aes(fill = intensity), alpha = 1) +
            scale_fill_gradient(low='black', high= 'white', limits=c(0,150), oob=squish, guide = "none") +
            # geom_path(aes(group=interaction(label, obj)), data=cont_df, size=0.1, color="red") +
            coord_fixed(ratio = 1) +
            facet_wrap(~slice_index, ncol=1) +
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
            facet_wrap(~slice_index, ncol=1) +
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

    micros = plot_micro(micro_final %>% filter(mask_value > 0.5), contour_final)
    # flair = plot_flair(flair_final %>% filter(mask_value > 0.5), label_ukb_final)
    flair = plot_flair(flair_final %>% filter(mask_value > 0.5))
    bison = plot_bison(label_final %>% filter(mask_value > 0.5))

    # Combine image types with patchwork
    plt_layout <- "
    ABCCCCCCCCC
    ABCCCCCCCCC
    "
    
    final_plot <- flair + bison + micros + plot_layout(design=plt_layout)

    png(paste0(output, "_raw.png"), width=3000, height=1500, res=300)
    print(final_plot)
    dev.off()
    print(paste0(output, "_raw.png"))

    # plot_annotation(
    #     title = ""
    # )
}

annot_demo_dx = function(df, output) {
    id_line = paste0("ID: ", df$ID[1], "     Sex: ",df$Sex[1],"     Age: ",df$Age[1])
    
    plt = ggplot() + 
        geom_text(aes(x=0, y=0.8, label=id_line), size=3) + 
        theme_void()
    
    if (nrow(df)>1) {
        dx_lines <- data.frame(
            x = rep(0, nrow(df)),
            y = seq(from=0.7, to=(0.8-(nrow(df)/10)), by = -0.1),
            label = paste0(df$icd_code, ": ", df$dx_name, " (", round(df$days_mri_dx/365,2), " years)")
            )

        plt = plt + geom_text(data = dx_lines, aes(x = x, y = y, label = label), size = 2)
    }

    if (nrow(df) == 1 & is.na(df$dx_name) == FALSE) {
        dx_lines <- data.frame(
            x = rep(0, nrow(df)),
            y = 0.7,
            label = paste0(df$icd_code, ": ", df$dx_name, " (", round(df$days_mri_dx/365,2), " years)")
            )

        plt = plt + geom_text(data = dx_lines, aes(x = x, y = y, label = label), size = 2)
    }

    ggsave(paste0(output, ".png"), plot = plt, bg="white", width=3000, height=600, units="px")
}

# Run
cores = detectCores()
parallelCluster <- makeCluster(cores-5, type = "SOCK", methods = FALSE)
setDefaultCluster(parallelCluster)
registerDoParallel(parallelCluster)

foreach(i=1:length(ids), .packages = c('RMINC', 'RColorBrewer', 'tidyverse', 'scales', 'patchwork', 'ggnewscale', 'viridis')) %dopar% {
# for (i in 1:length(ids)) {
    cat(paste0("\nID = ", ids[i], "\n"))

    dx_id = dx[which(dx$ID %in% ids[i]),]
    # dx_id <- dx_id[order(substr(dx_id$icd_code, 1, 1), as.numeric(substr(dx_id$icd_code, 2, 3))), ]
    dx_id <- dx_id[order(as.numeric(dx_id$days_mri_dx)), ]

    id_string = paste0(ids[i], "_", dx_id$Age[1], ifelse(dx_id$Sex[1] == "Female", "F", "M"))
    if(is.na(dx_id$icd_code[1]) == FALSE) {
        for (d in 1:length(dx_id$icd_code)) {
            id_string = paste0(id_string, "_", dx_id$icd_code[d])
        }
    }

    vis_dir=paste0("./visualization/indiv/",id_string)
    dir.create(vis_dir, showWarnings=FALSE)

    # Make PNG
    thresh = 3
    input_z = paste0("./results/mnc/",ids[i],"_zscore_anlm")
    output_z = paste0(vis_dir,"/",id_string)
    zscore_png(input_z, output_z, thresh)

    input_r = paste0("../../maps_UKB_space_anlm_all/sub-",ids[i],"_ses-2")
    output_r = paste0(vis_dir,"/",id_string)
    raw_png(input_r, output_r, "_anlm")

    output_a = paste0(vis_dir,"/",id_string,"_annot")
    annot_demo_dx(dx_id, output_a)

    # Paste pngs together
    command = paste0("convert -gravity Center ",output_a, ".png ", output_z,"_zscore_max",thresh,".png ",output_r,"_raw.png -append ",vis_dir,"/",id_string,"_all.png")
    system(command)
    print(paste0(vis_dir,"/",id_string,"_all.png"))
}

stopCluster(parallelCluster)
