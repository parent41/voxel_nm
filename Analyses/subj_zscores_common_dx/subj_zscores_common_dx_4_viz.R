
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

names = c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")
names_long = c("dti_MD", "NODDI_ISOVF", "dti_FA", "NODDI_ICVF", "NODDI_OD", "T2star", "QSM")

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
mask_path = "../../../UKB/temporary_template/Mask_2mm_dil2.mnc"

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

    plot_flair = function(df) {
        flair = ggplot(data = df, aes(x = x_toplot, y = y_toplot)) +
            geom_raster(aes(fill = intensity), alpha = 1) +
            scale_fill_gradient(low='black', high= 'white', limits=c(0,150), oob=squish, guide = "none") +
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
    flair = plot_flair(flair_final %>% filter(mask_value > 0.5))
    bison = plot_bison(label_final %>% filter(mask_value > 0.5))

    # Combine image types with patchwork
    plt_layout <- "
    ABCCCCCCC
    ABCCCCCCC
    "
    
    final_plot <- flair + bison + zscores + plot_layout(design=plt_layout)

    png(paste0(output, "_zscore_max",up_thresh,".png"), width=2500, height=1500, res=300)
    print(final_plot)
    dev.off()
}

raw_png = function(input, output) {

    raw_up_thresh = c(0.003, 0.5, 0.85, 0.85, 1, 75, 150)
    raw_down_thresh = c(0, 0, 0, 0, 0, 0, -100)

    input_id = sub(".*/(\\d{7}).*", "\\1", input)
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

        micro = prepare_masked_anatomy(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-",input_id,"_ses-2_",names_long[1],"_UKB.mnc"), mask_path, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]])
        ifelse(slices_axes[a] == "x", micro$anatomy_df$x_toplot <- micro$anatomy_df$y, ifelse(slices_axes[a] == "y", micro$anatomy_df$x_toplot <- micro$anatomy_df$x, ifelse(slices_axes[a] == "z", micro$anatomy_df$x_toplot <- micro$anatomy_df$x, NA)))
        ifelse(slices_axes[a] == "x", micro$anatomy_df$y_toplot <- micro$anatomy_df$z, ifelse(slices_axes[a] == "y", micro$anatomy_df$y_toplot <- micro$anatomy_df$z, ifelse(slices_axes[a] == "z", micro$anatomy_df$y_toplot <- micro$anatomy_df$y, NA)))
        micro$anatomy_df$slice_index = micro$anatomy_df$slice_index + n_slices
        micro$anatomy_df$intensity = rescale_vector(micro$anatomy_df$intensity, raw_up_thresh[1], raw_down_thresh[1])

        micro$anatomy_df$Micro = names[1]
        for (n in 2:length(names)) {
            micro_tmp = prepare_masked_anatomy(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-",input_id,"_ses-2_",names_long[n],"_UKB.mnc"), mask_path, slice_axis = slices_axes[a], slice_axis_coordinates = slices[[a]])
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

    plot_flair = function(df) {
        flair = ggplot(data = df, aes(x = x_toplot, y = y_toplot)) +
            geom_raster(aes(fill = intensity), alpha = 1) +
            scale_fill_gradient(low='black', high= 'white', limits=c(0,150), oob=squish, guide = "none") +
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
    flair = plot_flair(flair_final %>% filter(mask_value > 0.5))
    bison = plot_bison(label_final %>% filter(mask_value > 0.5))

    # Combine image types with patchwork
    plt_layout <- "
    ABCCCCCCC
    ABCCCCCCC
    "
    
    final_plot <- flair + bison + micros + plot_layout(design=plt_layout)

    png(paste0(output, "_raw.png"), width=2500, height=1500, res=300)
    # png(paste0("test_raw.png"), width=2500, height=1500, res=300)
    print(final_plot)
    dev.off()
}

# Run
cores = detectCores()
parallelCluster <- makeCluster(cores-5, type = "SOCK", methods = FALSE)
setDefaultCluster(parallelCluster)
registerDoParallel(parallelCluster)

up_thresholds <- c(3)

for (d in 1:ncol(dx)) {
    cat(paste0("\nDiagnosis = ", colnames(dx)[d]))
    dx_ids = dx[,d]
    dx_ids = na.omit(dx_ids)

    dir.create(paste0("./visualization/",colnames(dx)[d]), showWarnings=FALSE)

    foreach(i=1:length(dx_ids), .packages = c('RMINC', 'RColorBrewer', 'tidyverse', 'scales', 'patchwork', 'ggnewscale')) %dopar% {
    # for (i in 1:length(dx_ids)) {
        cat(paste0("\n\tID = ", dx_ids[i], "\t"))
        id_row = which(inclusions$ID == dx_ids[i])

        # If dx ID is in inclusions
        if(length(id_row)>0) {

            res_dir = paste0("./results/",colnames(dx)[d], "/", dx_ids[i])
            viz_dir = paste0("./visualization/",colnames(dx)[d], "/", dx_ids[i])
            dir.create(viz_dir, showWarnings=FALSE)
            
            # Make PNG (raw)
            input = paste0(res_dir, "/", dx_ids[i])
            output = paste0(viz_dir, "/", dx_ids[i])
            input_anlm = paste0(res_dir, "/", dx_ids[i],"_anlm")
            output_anlm = paste0(viz_dir, "/", dx_ids[i], "_anlm")
            # For every threshold
            # for (t in 1:length(down_thresholds)) {
            for (t in 1:length(up_thresholds)) {
                # Not denoised
                zscore_png(input, output, up_thresholds[t])
                cat(paste0("\n\tZ-score PNG = ", output, "_zscore_max",up_thresholds[t],".png"))
                # Denoised
                zscore_png(input_anlm, output_anlm, up_thresholds[t])
                cat(paste0("\n\tZ-score PNG = ", output_anlm, "_zscore_max",up_thresholds[t],".png"))
            }
            raw_png(input, output)
            cat(paste0("\n\tRaw PNG = ", output, "_raw.png"))
            raw_png(input_anlm, output_anlm)
            cat(paste0("\n\tRaw PNG = ", output_anlm, "_raw.png"))
        }
    }
}

stopCluster(parallelCluster)
















# # Old function with MRIcrotome

# zscore_png = function(input, output, up_thresh, down_thresh) {
#     input_id = sub(".*/(\\d{7}).*", "\\1", input)
#     print(input_id)

#     color_scale_div_2 = colorRampPalette(brewer.pal(9,"Blues"))(255)
#     color_scale_div_1 = colorRampPalette(brewer.pal(9,"Reds"))(255)

#     anatVol = mincArray(mincGetVolume(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-",input_id,"_ses-2_flair_UKB.mnc")))
#     bisonVol = mincArray(mincGetVolume(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-", input_id,"_ses-2_Label_UKB.mnc")))
#     micro_subj = list()
#     for (n in 1:length(names)) {
#         # print(names[n])
#         micro_subj[[n]] = mincArray(mincGetVolume(paste0(input, "_",names[n],".mnc")))
#     }

#     # png(file=paste0(output, "_all_zscore_min",down_thresh,"_max",up_thresh,".png"), width=2000, height=4000, pointsize = 150)
#     # sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
#     #     addtitle("FLAIR") %>%
#     #     anatomy(anatVol, low=10, high=200) %>%
#     # sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
#     #     addtitle(names[1]) %>%
#     #     anatomy(anatVol, low=10, high=200) %>%
#     #     overlay(micro_subj[[1]], low=down_thresh, high=up_thresh,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
#     #     legend("Z-values") %>%
#     #     contours(bisonVol, levels=c(7,8), col=c("green"), lwd=2) %>%
#     #     # contours(bisonVol, levels=c(2), col=c("green"), lwd=10) %>%
#     # draw()
#     # dev.off()
#     # print(paste0(output, "_all_zscore_min",down_thresh,"_max",up_thresh,".png"))

#     png(file=paste0(output, "_all_zscore_min",down_thresh,"_max",up_thresh,".png"), width=8500, height=4000, pointsize = 150)
#     # png(file=paste0("test_mricrotome.png"), width=8500, height=4000, pointsize = 150)
#     sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
#         addtitle(names[1]) %>%
#         anatomy(anatVol, low=10, high=200) %>%
#         overlay(micro_subj[[1]], low=down_thresh, high=up_thresh,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
#     sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
#         addtitle(names[2]) %>%
#         anatomy(anatVol, low=10, high=200) %>%
#         overlay(micro_subj[[2]], low=down_thresh, high=up_thresh,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
#     sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
#         addtitle(names[3]) %>%
#         anatomy(anatVol, low=10, high=200) %>%
#         overlay(micro_subj[[3]], low=down_thresh, high=up_thresh,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
#     sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
#         addtitle(names[4]) %>%
#         anatomy(anatVol, low=10, high=200) %>%
#         overlay(micro_subj[[4]], low=down_thresh, high=up_thresh,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
#     sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
#         addtitle(names[5]) %>%
#         anatomy(anatVol, low=10, high=200) %>%
#         overlay(micro_subj[[5]], low=down_thresh, high=up_thresh,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
#     sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
#         addtitle(names[6]) %>%
#         anatomy(anatVol, low=10, high=200) %>%
#         overlay(micro_subj[[6]], low=down_thresh, high=up_thresh,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
#     sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
#         addtitle(names[7]) %>%
#         anatomy(anatVol, low=10, high=200) %>%
#         overlay(micro_subj[[7]], low=down_thresh, high=up_thresh,col=color_scale_div_1, rCol=color_scale_div_2, symmetric=TRUE) %>%
#     legend("Z-values") %>%
#     draw()
#     dev.off()
# }



# raw_png = function(input, output) {
#     input_id = sub(".*/", "", input)

#     color_scale_div_2 = colorRampPalette(brewer.pal(9,"Blues"))(255)
#     color_scale_div_1 = colorRampPalette(brewer.pal(9,"Reds"))(255)

#     down_lims = c(0.05,0.05,0.05,0.05,0.05,0.7,0.001)
#     up_lims = c(0.95,0.95,0.95,0.95,0.95,0.99,0.999)

#     anatVol = mincArray(mincGetVolume(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-",input_id,"_ses-2_flair_UKB.mnc")))
#     micro_subj = list()
#     for (n in 1:length(names)) {
#         # print(names[n])
#         micro_subj[[n]] = mincArray(mincGetVolume(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-",input_id,"_ses-2_",names_long[n],"_UKB.mnc")))
#     }

#     png(file=paste0(output, "_all_raw.png"), width=8500, height=4000, pointsize = 150)
#     sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
#         addtitle("FLAIR") %>%
#         anatomy(anatVol, low=10, high=200) %>%
#     sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
#         addtitle(names[1]) %>%
#         anatomy(micro_subj[[1]], low=quantile(micro_subj[[1]][][micro_subj[[1]][]!=0],down_lims[1]), high=quantile(micro_subj[[1]][][micro_subj[[1]][]!=0], up_lims[1])) %>%
#     sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
#         addtitle(names[2]) %>%
#         anatomy(micro_subj[[2]], low=quantile(micro_subj[[2]][][micro_subj[[2]][]!=0],down_lims[2]), high=quantile(micro_subj[[2]][][micro_subj[[2]][]!=0], up_lims[2])) %>%
#     sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
#         addtitle(names[3]) %>%
#         anatomy(micro_subj[[3]], low=quantile(micro_subj[[3]][][micro_subj[[3]][]!=0],down_lims[3]), high=quantile(micro_subj[[3]][][micro_subj[[3]][]!=0], up_lims[3])) %>%
#     sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
#         addtitle(names[4]) %>%
#         anatomy(micro_subj[[4]], low=quantile(micro_subj[[4]][][micro_subj[[4]][]!=0],down_lims[4]), high=quantile(micro_subj[[4]][][micro_subj[[4]][]!=0], up_lims[4])) %>%
#     sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
#         addtitle(names[5]) %>%
#         anatomy(micro_subj[[5]], low=quantile(micro_subj[[5]][][micro_subj[[5]][]!=0],down_lims[5]), high=quantile(micro_subj[[5]][][micro_subj[[5]][]!=0], up_lims[5])) %>%
#     sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
#         addtitle(names[6]) %>%
#         anatomy(micro_subj[[6]], low=quantile(micro_subj[[6]][][micro_subj[[6]][]!=0],down_lims[6]), high=quantile(micro_subj[[6]][][micro_subj[[6]][]!=0], up_lims[6])) %>%
#     sliceSeries(nrow=4, ncol=1, begin=42, end=52, dimension=3) %>%
#         addtitle(names[7]) %>%
#         anatomy(micro_subj[[7]], low=quantile(micro_subj[[7]][][micro_subj[[7]][]!=0],down_lims[7]), high=quantile(micro_subj[[7]][][micro_subj[[7]][]!=0], up_lims[7])) %>%
#     legend("Micro") %>%
#     draw()
#     dev.off()
# }

