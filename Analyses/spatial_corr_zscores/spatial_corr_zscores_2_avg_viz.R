
library(data.table)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)
library(Hmisc)
library(corrplot)

names = c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM", "jacobians_rel", "jacobians_abs")
regions = c('Whole_brain', 'Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH', 'Cerebral_WM')

# Load data

corr_r = as.data.frame(fread("./results/indiv_concat_raw.tsv"))
corr_z = as.data.frame(fread("./results/indiv_concat_zscores.tsv"))

# corr_r[,c(1,2,3,6)] = as.data.frame(lapply(corr_r[,c(1,2,3,6)], as.factor))
# corr_z[,c(1,2,3,6)] = as.data.frame(lapply(corr_z[,c(1,2,3,6)], as.factor))

mypalette = colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                   "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                   "#4393C3", "#2166AC", "#053061")))


# Correlation matrices by region

corrmat_by_reg = function(df, type) {

    # Correlation matrices for each region
    for (l in 1:length(levels(factor(df$Label)))) {
        print(levels(factor(df$Label))[l])

        df_corr = df %>%
            filter(Label == levels(factor(df$Label))[l]) %>%
            bind_rows(
                df %>%
                    filter(Label == levels(factor(df$Label))[l]) %>%
                    mutate(across(c(Marker_1, Marker_2), ~ ifelse(Marker_1 == ., Marker_2, Marker_1)))
            ) %>%
            mutate(marker_1_2 = factor(paste0(Marker_1,"_with_",Marker_2))) %>%
            group_by(marker_1_2) %>%
            dplyr::summarize(mean = mean(cor, na.rm=TRUE)) %>%
            separate(marker_1_2, into = c("Marker_1", "Marker_2"), sep = "_with_") %>%
            arrange(factor(Marker_1, levels = names), factor(Marker_2, levels = names)) %>%
            pivot_wider(names_from = Marker_2, values_from = mean) %>%
            column_to_rownames("Marker_1") %>%
            select(names) #%>%
            # glimpse()

        png(height=20, width=20, units="in", res=300, pointsize=10, file=paste0("./visualization/corrmat_by_region/corrmat_",type,"_",levels(factor(df$Label))[l],".png"))
        corrplot(as.matrix(df_corr), type="full", order="original", method="circle", addCoef.col = "black", addgrid.col = TRUE,
                tl.srt = 45, tl.col = "black", col=mypalette(200), diag=FALSE, tl.cex=3, cl.length = 5, cl.cex = 3, number.cex=3, cex.main=5,
                title=levels(factor(df$Label))[l], mar=c(0,0,4,0), )
        dev.off()
    }

    # Append correlation matrices
    command = paste0("convert -gravity center ")
    for (l in 1:length(levels(factor(df$Label, levels=regions)))) {
        command = paste0(command, "./visualization/corrmat_by_region/corrmat_",type,"_",levels(factor(df$Label, levels=regions))[l],".png ")
    }
    command = paste0(command, "-append ./visualization/corrmat_by_region/corrmat_",type,"_allregions.png")
    system(command)
}

dir.create("./visualization/corrmat_by_region", showWarnings=FALSE)

# Run for raw and zscores
corrmat_by_reg(corr_r, "raw")
corrmat_by_reg(corr_z, "zscores")

# Append
command = paste0("convert -gravity center ./visualization/corrmat_by_region/corrmat_raw_allregions.png ./visualization/corrmat_by_region/corrmat_zscores_allregions.png +append ./visualization/corrmat_by_region/corrmat_all_allregions.png")
system(command)

# Interactions between regions within each marker

slopes_across_reg = function(df, type) {

    df_tmp = df %>%
        bind_rows(
            df %>% mutate(across(c(Marker_1, Marker_2), ~ ifelse(Marker_1 == ., Marker_2, Marker_1)))
        ) %>%
        mutate(marker_1_2 = factor(paste0(Marker_1,"_with_",Marker_2))) %>%
        group_by(marker_1_2, Label) %>%
        dplyr::summarize(mean = mean(cor, na.rm=TRUE)) %>%
        separate(marker_1_2, into = c("Marker_1", "Marker_2"), sep = "_with_") %>%
        arrange(factor(Marker_1, levels = names), factor(Marker_2, levels = names)) #%>%
        # pivot_wider(names_from = Marker_2, values_from = mean) %>%
        # column_to_rownames("Marker_1") %>%
        # select(names) #%>%
        # glimpse()

    df_tmp[,c(1,2,3)] = as.data.frame(lapply(df_tmp[,c(1,2,3)], as.factor))

    ggplot(df_tmp, aes(x = seq(-1,1, by=0.01), y = seq(-1,1, by=0.01))) +
        geom_abline(aes(slope = mean, intercept=0, color=Label), size=0.3, alpha=0.7) +
        geom_abline(data = df_tmp[df_tmp$Label == "Whole_brain",], aes(slope = mean, intercept=0), color="black", size=0.8, alpha=0.7) +
        facet_grid(factor(Marker_1, levels=names) ~ factor(Marker_2, levels=names)) + 
        scale_y_continuous(limits = c(-1,1), breaks=c(-1,0,1), name="") + 
        scale_x_continuous(limits = c(-1,1), breaks=c(-1,0,1), name="")
    ggsave(paste0("./visualization/slopes_across_reg/slopes_across_reg_",type,".png"), width=12, height=10)
}

dir.create("./visualization/slopes_across_reg", showWarnings=FALSE)

slopes_across_reg(corr_r, "raw")
slopes_across_reg(corr_z, "zscores")




