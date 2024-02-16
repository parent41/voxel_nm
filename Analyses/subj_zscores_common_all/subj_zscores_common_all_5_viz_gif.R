
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

names = c("MD","FA")
names_long = c("dti_MD", "dti_FA")

# Load data
inclusions = as.data.frame(fread("../../../WMH_micro_spatial/QC/inclusions_without_excluding_dx_new.txt"))
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

# Source yohan's ggplot for MRI extension
source("../../yohan_ggplot/plotting_functions.R")
source("../../yohan_ggplot/plotting_functions_labels.R")

# Select slices
slices = list()
slices[["z"]] = seq(-70, 90, by=2) # Slices for z-axis
slices[["x"]] = seq(-80, 80, by=2) # Slices for x-axis (sagital)
slices[["y"]] = seq(-90, 70, by=2) # Slices for y-axis (coronal)

rm(demo)

# Make png files for each subject grouped by dx

zscore_png = function(input, output, up_thresh, axis) {
    input_id = sub(".*/(\\d{7}).*", "\\1", input)
    # print(input_id)

    # Load data for slices across different axes
    contour_final = tibble(matrix(nrow = 0, ncol = 11))
    micro_final = tibble(matrix(nrow = 0, ncol = 12))
    
    n_slices = 0
    for (a in 1:length(axis)) {
        # print(axis)
        slices_num = seq(n_slices + 1, n_slices + length(slices[[axis]]))

        # Load image data
        contour = prepare_label_contours(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-", input_id,"_ses-2_Label_UKB.mnc"), mask_path, slice_axis = axis, slice_axis_coordinates = slices[[axis]], labels = seq(3,9))
        ifelse(axis == "x", contour$contours_df$x_toplot <- contour$contours_df$y, ifelse(axis == "y", contour$contours_df$x_toplot <- contour$contours_df$x, ifelse(axis == "z", contour$contours_df$x_toplot <- contour$contours_df$x, NA)))
        ifelse(axis == "x", contour$contours_df$y_toplot <- contour$contours_df$z, ifelse(axis == "y", contour$contours_df$y_toplot <- contour$contours_df$z, ifelse(axis == "z", contour$contours_df$y_toplot <- contour$contours_df$y, NA)))
        # contour$contours_df$slice_index = contour$contours_df$slice_index + n_slices
        contour_final = rbind(contour_final, contour$contours_df)

        micro = prepare_masked_anatomy(paste0(input, "_",names[1],".mnc"), mask_path, slice_axis = axis, slice_axis_coordinates = slices[[axis]])
        ifelse(axis == "x", micro$anatomy_df$x_toplot <- micro$anatomy_df$y, ifelse(axis == "y", micro$anatomy_df$x_toplot <- micro$anatomy_df$x, ifelse(axis == "z", micro$anatomy_df$x_toplot <- micro$anatomy_df$x, NA)))
        ifelse(axis == "x", micro$anatomy_df$y_toplot <- micro$anatomy_df$z, ifelse(axis == "y", micro$anatomy_df$y_toplot <- micro$anatomy_df$z, ifelse(axis == "z", micro$anatomy_df$y_toplot <- micro$anatomy_df$y, NA)))
        # micro$anatomy_df$slice_index = micro$anatomy_df$slice_index + n_slices

        micro$anatomy_df$Micro = names[1]
        for (n in 2:length(names)) {
            micro_tmp = prepare_masked_anatomy(paste0(input, "_",names[n],".mnc"), mask_path, slice_axis = axis, slice_axis_coordinates = slices[[axis]])
            ifelse(axis == "x", micro_tmp$anatomy_df$x_toplot <- micro_tmp$anatomy_df$y, ifelse(axis == "y", micro_tmp$anatomy_df$x_toplot <- micro_tmp$anatomy_df$x, ifelse(axis == "z", micro_tmp$anatomy_df$x_toplot <- micro_tmp$anatomy_df$x, NA)))
            ifelse(axis == "x", micro_tmp$anatomy_df$y_toplot <- micro_tmp$anatomy_df$z, ifelse(axis == "y", micro_tmp$anatomy_df$y_toplot <- micro_tmp$anatomy_df$z, ifelse(axis == "z", micro_tmp$anatomy_df$y_toplot <- micro_tmp$anatomy_df$y, NA)))
            # micro_tmp$anatomy_df$slice_index = micro_tmp$anatomy_df$slice_index + n_slices

            micro_tmp$anatomy_df$Micro = names[n]
            micro$anatomy_df = rbind(micro$anatomy_df, micro_tmp$anatomy_df)
        }
        micro$anatomy_df$Micro = as.factor(micro$anatomy_df$Micro)
        micro_final = rbind(micro_final, micro$anatomy_df)

        # Keep track of total number of slices
        n_slices = n_slices + length(slices[[axis]])
    }

    micro_final$Micro = factor(micro_final$Micro, levels=names)

    # For z-axis, invert slice index to go from bottom to top of brain
    if (axis == "z") {
        micro_final$slice_index = max(micro_final$slice_index) + 1 - micro_final$slice_index
        contour_final$slice_index = (max(contour_final$slice_index) + 1 - contour_final$slice_index)
    }

    # Make zscore images

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

    # Make temporary pngs for each slice
    for (i in 1:length(unique(micro_final$slice_index))) {
        print(i)

        zscores = plot_zscore(micro_final %>% filter(slice_index == i), contour_final %>% filter(slice_index == i))

        png(paste0(output,"_slice_",i,".png"), width=1500, height=1000, res=300)
        print(zscores)
        dev.off()
        print(paste0(output,"_slice_",i,".png"))
    }
}

raw_png = function(input, output, anlm, axis) {

    raw_up_thresh = c(0.003, 0.85)
    raw_down_thresh = c(0, 0)

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
    slices_axes = c("z", "x", "y")

    contour_final = tibble(matrix(nrow = 0, ncol = 11))
    micro_final = tibble(matrix(nrow = 0, ncol = 12))

    n_slices = 0
    for (a in 1:length(axis)) {
        # print(axis)
        slices_num = seq(n_slices + 1, n_slices + length(slices[[axis]]))

        # Load image data
        contour = prepare_label_contours(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-", input_id,"_ses-2_Label_UKB.mnc"), mask_path, slice_axis = axis, slice_axis_coordinates = slices[[axis]], labels = seq(3,9))
        ifelse(axis == "x", contour$contours_df$x_toplot <- contour$contours_df$y, ifelse(axis == "y", contour$contours_df$x_toplot <- contour$contours_df$x, ifelse(axis == "z", contour$contours_df$x_toplot <- contour$contours_df$x, NA)))
        ifelse(axis == "x", contour$contours_df$y_toplot <- contour$contours_df$z, ifelse(axis == "y", contour$contours_df$y_toplot <- contour$contours_df$z, ifelse(axis == "z", contour$contours_df$y_toplot <- contour$contours_df$y, NA)))
        # contour$contours_df$slice_index = contour$contours_df$slice_index + n_slices
        contour_final = rbind(contour_final, contour$contours_df)

        # micro = prepare_masked_anatomy(list.files(input_dir, pattern = paste0("sub-", input_id, "_ses-2_",names_long[1],"*"), full.names = TRUE), mask_path, slice_axis = axis, slice_axis_coordinates = slices[[axis]])
        micro = prepare_masked_anatomy(paste0(input, "_",names_long[1],"_UKB",anlm,".mnc"), mask=mask_id, slice_axis = axis, slice_axis_coordinates = slices[[axis]])
        ifelse(axis == "x", micro$anatomy_df$x_toplot <- micro$anatomy_df$y, ifelse(axis == "y", micro$anatomy_df$x_toplot <- micro$anatomy_df$x, ifelse(axis == "z", micro$anatomy_df$x_toplot <- micro$anatomy_df$x, NA)))
        ifelse(axis == "x", micro$anatomy_df$y_toplot <- micro$anatomy_df$z, ifelse(axis == "y", micro$anatomy_df$y_toplot <- micro$anatomy_df$z, ifelse(axis == "z", micro$anatomy_df$y_toplot <- micro$anatomy_df$y, NA)))
        # micro$anatomy_df$slice_index = micro$anatomy_df$slice_index + n_slices
        micro$anatomy_df$intensity = rescale_vector(micro$anatomy_df$intensity, raw_up_thresh[1], raw_down_thresh[1])

        micro$anatomy_df$Micro = names[1]
        for (n in 2:length(names)) {
            micro_tmp = prepare_masked_anatomy(paste0(input, "_",names_long[n],"_UKB",anlm,".mnc"), mask=mask_id, slice_axis = axis, slice_axis_coordinates = slices[[axis]])
            ifelse(axis == "x", micro_tmp$anatomy_df$x_toplot <- micro_tmp$anatomy_df$y, ifelse(axis == "y", micro_tmp$anatomy_df$x_toplot <- micro_tmp$anatomy_df$x, ifelse(axis == "z", micro_tmp$anatomy_df$x_toplot <- micro_tmp$anatomy_df$x, NA)))
            ifelse(axis == "x", micro_tmp$anatomy_df$y_toplot <- micro_tmp$anatomy_df$z, ifelse(axis == "y", micro_tmp$anatomy_df$y_toplot <- micro_tmp$anatomy_df$z, ifelse(axis == "z", micro_tmp$anatomy_df$y_toplot <- micro_tmp$anatomy_df$y, NA)))
            # micro_tmp$anatomy_df$slice_index = micro_tmp$anatomy_df$slice_index + n_slices
            micro_tmp$anatomy_df$intensity = rescale_vector(micro_tmp$anatomy_df$intensity, raw_up_thresh[n], raw_down_thresh[n])

            micro_tmp$anatomy_df$Micro = names[n]
            micro$anatomy_df = rbind(micro$anatomy_df, micro_tmp$anatomy_df)
        }
        micro$anatomy_df$Micro = as.factor(micro$anatomy_df$Micro)
        micro_final = rbind(micro_final, micro$anatomy_df)

        # Keep track of total number of slices
        n_slices = n_slices + length(slices[[axis]])
    }

    micro_final$Micro = factor(micro_final$Micro, levels=names)

    # For z-axis, invert slice index to go from bottom to top of brain
    if (axis == "z") {
        micro_final$slice_index = max(micro_final$slice_index) + 1 - micro_final$slice_index
        contour_final$slice_index = (max(contour_final$slice_index) + 1 - contour_final$slice_index)
    }

    # Make invidiual plots by image types

    plot_micro = function(df, cont_df) {

        df$intensity = as.numeric(df$intensity)

        zscores = ggplot(data=df, aes(x = x_toplot, y = y_toplot)) +
            geom_tile(aes(fill = intensity), alpha = 1) +
            # scale_fill_gradient(low="white", high="red", limits=c(0,1), breaks=c(0,1), oob=squish, name="Micro") +
            scale_fill_viridis(option="A", limits=c(0,1), breaks=c(0,1), oob=squish, name="Micro", na.value = "white") +
            geom_path(aes(group=interaction(label, obj)), data=cont_df, size=0.1, color="black") +
            coord_fixed(ratio = 1) +
            facet_grid(slice_axis~Micro) +
            labs(title = "Raw") +
            theme_void() +
            theme(plot.background = element_rect(fill = "white"),
                plot.title = element_text(hjust = 0.5),
                strip.text.y=element_text(size=0),
                strip.text.x=element_text(size=12),
                panel.spacing = unit(0, "lines"),
                legend.position="none")
        print(zscores)

        return(zscores)
    }

    # Make temporary pngs for each slice
    for (i in 1:length(unique(micro_final$slice_index))) {
        print(i)

        micros = plot_micro(micro_final %>% filter(slice_index == i)
                                        %>% mutate(intensity = ifelse(mask_value == 0, NA, intensity)),
                            contour_final %>% filter(slice_index == i))

        png(paste0(output,"_slice_",i,".png"), width=1500, height=1000, res=300)
        print(micros)
        dev.off()
        print(paste0(output,"_slice_",i,".png"))
    }
}

flair_png = function(input, output, axis) {
    input_id = sub(".*/(\\d{7}).*", "\\1", input)
    # print(input_id)

    # Load data for slices across different axes
    slices_axes = c("z", "x", "y")

    n_slices = 0

    flair_final = tibble(matrix(nrow = 0, ncol = 11))

    for (a in 1:length(axis)) {
        # print(axis)
        slices_num = seq(n_slices + 1, n_slices + length(slices[[a]]))

        # Load image data
        flair = prepare_masked_anatomy(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-",input_id,"_ses-2_flair_UKB.mnc"), mask_path, slice_axis = axis, slices[[axis]])
        ifelse(axis == "x", flair$anatomy_df$x_toplot <- flair$anatomy_df$y, ifelse(axis == "y", flair$anatomy_df$x_toplot <- flair$anatomy_df$x, ifelse(axis == "z", flair$anatomy_df$x_toplot <- flair$anatomy_df$x, NA)))
        ifelse(axis == "x", flair$anatomy_df$y_toplot <- flair$anatomy_df$z, ifelse(axis == "y", flair$anatomy_df$y_toplot <- flair$anatomy_df$z, ifelse(axis == "z", flair$anatomy_df$y_toplot <- flair$anatomy_df$y, NA)))
        flair$anatomy_df$slice_index = flair$anatomy_df$slice_index + n_slices
        flair_final = rbind(flair_final, flair$anatomy_df)

        # Keep track of total number of slices
        n_slices = n_slices + length(slices[[a]])
    }

    # For z-axis, invert slice index to go from bottom to top of brain
    if (axis == "z") {
        flair_final$slice_index = max(flair_final$slice_index) + 1 - flair_final$slice_index
    }


    plot_flair = function(df) {
        flair = ggplot(data = df, aes(x = x_toplot, y = y_toplot)) +
            geom_raster(aes(fill = intensity), alpha = 1) +
            scale_fill_gradient(low='black', high= 'white', limits=c(0,150), oob=squish, guide = "none") +
            # geom_path(aes(group=interaction(label, obj)), data=cont_df, size=0.1, color="red") +
            coord_fixed(ratio = 1) +
            # facet_wrap(~slice_index, ncol=1) +
            labs(title = "FLAIR") +
            theme_void() +
            theme(plot.background = element_rect(fill = "white"),
                plot.title = element_text(hjust = 0.5, size=20),
                strip.text.y=element_text(size=12),
                strip.text.x=element_text(size=0),
                panel.spacing = unit(0, "lines"))
        # ggsave("test_flair.png")

        return(flair)
    }

    # Make temporary pngs for each slice
    for (i in 1:length(unique(flair_final$slice_index))) {
        print(i)

        flair = plot_flair(flair_final %>% filter(slice_index == i))

        png(paste0(output,"_slice_",i,".png"), width=1000, height=1500, res=300)
        print(flair)
        dev.off()
        print(paste0(output,"_slice_",i,".png"))
    }
}


# Run

cores = detectCores()
parallelCluster <- makeCluster(cores-5, type = "SOCK", methods = FALSE)
setDefaultCluster(parallelCluster)
registerDoParallel(parallelCluster)

id_list = c("1235478", "1776711", "4519919", "5836644", "1421334", "1607044")
axes_to_viz = c("z", "z", "z", "z", "x", "x")

viz_dir = paste0("./visualization/anim")
dir.create(viz_dir, showWarnings=FALSE)

foreach(i=1:length(id_list), .packages = c('RMINC', 'RColorBrewer', 'tidyverse', 'scales', 'patchwork', 'ggnewscale', 'viridis')) %dopar% {
# for (i in 1:length(id_list)) {

    print(id_list[i])

    dx_id = dx[which(dx$ID %in% id_list[i]),]
    dx_id <- dx_id[order(substr(dx_id$icd_code, 1, 1), as.numeric(substr(dx_id$icd_code, 2, 3))), ]

    id_string = paste0(id_list[i], "_", dx_id$Age[1], ifelse(dx_id$Sex[1] == "Female", "F", "M"))
    if(is.na(dx_id$icd_code[1]) == FALSE) {
        for (d in 1:length(dx_id$icd_code)) {
            id_string = paste0(id_string, "_", dx_id$icd_code[d])
        }
    }

    # Make PNGs
    viz_tmp = paste0("./tmp/",id_string)
    dir.create(viz_tmp, showWarnings=FALSE)

    # Zscores
    print("Zscores PNGs")
    thresh = 3
    input_z = paste0("./results/mnc/", id_list[i],"_zscore_anlm")
    output_z = paste0(viz_tmp, "/",id_string,"_zscore_anlm")
    zscore_png(input_z, output_z, thresh, axes_to_viz[i])

    # Raw
    print("Raw PNGs")
    input_r = paste0("../../maps_UKB_space_anlm_all/sub-",id_list[i],"_ses-2")
    output_r = paste0(viz_tmp, "/",id_string,"_raw_anlm")
    raw_png(input_r, output_r, "_anlm", axes_to_viz[i])

    # FLAIR
    print("FLAIR PNGs")
    input_f = input_z
    output_f = paste0(viz_tmp, "/",id_string,"_flair")
    flair_png(input_f, output_f, axes_to_viz[i])

    # Concatenate images together for each slice
    print("Concatenate images")

    for (s in 1:length(slices[[axes_to_viz[i]]])) {
        print(paste0("Slice: ",s))
        command = paste0("convert -gravity center ",output_r, "_slice_",s,".png ",output_z, "_slice_",s,".png -append ",output_r,"_zscores_slice_",s,".png")
        system(command)
        command = paste0("convert -gravity center ",output_f,"_slice_",s,".png ",output_r,"_zscores_slice_",s,".png +append ",output_r,"_zscores_flair_slice_",s,".png")
        system(command)
    }

    # Make MP4 with ffmpeg
    print("Make MP4s and GIFs")

    command = paste0("ffmpeg -y -framerate 8 -i ",output_r,"_zscores_flair_slice_%00d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p ",viz_dir,"/",id_string,"_allslices_",axes_to_viz[i],".mp4")
    system(command)

    # Convert to GIF
    command = paste0("ffmpeg -y -i ",viz_dir,"/",id_string,"_allslices_",axes_to_viz[i],".mp4 -vf \"fps=8,scale=1920:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse\" -loop 0 ",viz_dir,"/",id_string,"_allslices_",axes_to_viz[i],".gif")
    system(command)

    # Delete tmp directory
    unlink(viz_tmp, recursive = TRUE)
}



