
library(RMINC)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(data.table)
library(ggpubr)


tissue_nm = c('Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM')
tissue_all=c('Ventricules', 'CSF', 'Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH')

names = c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")

# Load results
results_dx = as.data.frame(fread("./results/final_results_dx.tsv"))

results_dx[,c(2,3,7,8)] = as.data.frame(lapply(results_dx[,c(2,3,7,8)], as.factor))

dx_prev_order = as.data.frame(table(results_dx$dx))
dx_prev_order = dx_prev_order[order(-dx_prev_order$Freq),]
dx_prev_order$Freq = dx_prev_order$Freq / 147
dx_prev_order = subset(dx_prev_order, Freq >= 30, select="Var1")
dx_prev_order = as.character(dx_prev_order$Var1)

# Bar plots: one plot by micro and threshold (across dx and regions)

for (i in 1:length(levels(results_dx$threshold))) {
  print(levels(results_dx$threshold)[i])

  plots = list()
  for (n in 1:length(names)) {
    df_tmp = results_dx %>% filter(threshold == levels(results_dx$threshold)[i] & micro == names[n] & label_value != "9")
    dx_prev_tmp = as.data.frame(table(df_tmp$dx) / nrow(df_tmp %>% filter(ID == df_tmp$ID[1])))
    dx_prev_tmp = dx_prev_tmp[order(-dx_prev_tmp$Freq),]
    legend_labels = paste0(as.character(dx_prev_tmp$Var1), " (n=",dx_prev_tmp$Freq, ")")

    plots[[n]] = ggplot(results_dx %>% filter(threshold == levels(results_dx$threshold)[i] & micro == names[n] & label_value != "9"),
          aes(x=label_value, y=perc_vox_above_thresh, fill=factor(dx, levels=c(dx_prev_order)))) + 
          geom_boxplot(outlier.shape = NA) + 
          scale_fill_discrete(name="Diagnosis", labels=legend_labels) + 
          scale_x_discrete(labels = tissue_all[seq(3,8)], name="") + 
          scale_y_continuous(name = paste0("% abnormal voxels above Z=",levels(results_dx$threshold)[i]), limits = c(0, quantile(df_tmp$perc_vox_above_thresh, 0.90))) +
          ggtitle(paste0(names[n])) +
          theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5))
    ggsave(paste0("./visualization/perc_abn_vox_Z",levels(results_dx$threshold)[i],"_",names[n],".png"), width=20, height=5)
  }

  wrap_plots(plots, ncol=1)
  ggsave(paste0("./visualization/perc_abn_vox_Z",levels(results_dx$threshold)[i],"_allmicro.png"), width=20, height=(5*length(names)))
}

# Radar plots

