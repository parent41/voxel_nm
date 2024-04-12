
library(RMINC)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(data.table)
library(ggpubr)
library(scales)
# install.packages(c("fmsb"))
# install.packages("devtools")
# devtools::install_github("ricardo-bion/ggradar", dependencies = TRUE, lib="/gpfs/fs1/home/m/mchakrav/parent41/R/x86_64-pc-linux-gnu-library/4.1")
# library(fmsb)
# library(ggradar)

tissue_nm = c('Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM')
tissue_all=c('Ventricules', 'CSF', 'Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH', 'Cerebral_WM')
tissue_abn = c('Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH', 'Cerebral_WM')
tissue_toplot = c('Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_WM')

names = c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM", "jacobians_rel", "jacobians_abs")

# Load data
results = as.data.frame(fread("./results/results.tsv"))
inclusions = as.data.frame(fread("../../../WMH_micro_spatial/QC/inclusions_without_excluding_dx_new.txt"))
colnames(inclusions) = "ID"

demo = as.data.frame(fread("../../../UKB/tabular/df_demo/UKBB_demo_wider.tsv"))
demo = subset(demo, InstanceID == 2, select=c('SubjectID', 'Sex_31_0', 'Age_when_attended_assessment_centre_21003_0'))
colnames(demo) = c("ID", "Sex", "Age")
demo = merge(inclusions, demo, by="ID", all.x=TRUE)

df_dx = as.data.frame(fread("../../../WMH_micro_spatial/Analyses_nm/predict_firstocc_categ/results/firstocc_categ.tsv"))

no_dx_ids = as.data.frame(inclusions$ID[!(inclusions$ID %in% df_dx$ID)])
colnames(no_dx_ids) = "ID"
df_dx = merge(no_dx_ids, df_dx, by="ID", all=TRUE)

df = merge(demo, df_dx, all.y=TRUE, by="ID")
df[,c(2,5,6,9,10)] = as.data.frame(lapply(df[,c(2,5,6,9,10)], as.factor))

results[,c(2,3,7)] = as.data.frame(lapply(results[,c(2,3,7)], as.factor))

icd_codes_list <- readRDS("../../../WMH_micro_spatial/Analyses_nm/predict_firstocc_categ/results/icd_codes_list.rds")

lm_results = as.data.frame(fread("./results/lm_results.tsv"))
lm_results[,c(1,2,3,4)] = as.data.frame(lapply(lm_results[,c(1,2,3,4)], as.factor))

# Plot beta coefficients and pvalues

plot_beta_pval = function(dx_name, icd_list, name) {
    color_scale = c(MD = "#04319E", ISOVF = "#8298CF",
                    FA = "#2B520B" , ICVF= "#448312", OD="#A2C189",
                    T2star="#D36108", QSM="#E9B084",
                    jacobians_rel="#5900cb", jacobians_abs="#893cef")
    n_dx = length(unique(df$ID[which(df$icd_code %in% icd_list)]))

    df_toplot = subset(lm_results, dx == dx_name & Label %in% c(3,4,5,6,7,10))

    df_toplot$fdr_group = p.adjust(df_toplot$pval_group, method="fdr")

    plot_pval = ggplot(df_toplot, aes(x=Label, y=fdr_group, color=factor(Micro, levels=names), group = factor(Micro, levels=names))) + 
          # geom_point(position = position_jitter(w = 0.4, h = 0)) +
          geom_point(position = position_dodge(0.8)) +
          # geom_line(position = position_dodge(0.8), alpha=0.2) +
          scale_color_manual(name="Micro", values=color_scale) + 
          scale_x_discrete(labels = tissue_toplot, name="") +
          scale_y_continuous(limits = c(0, 0.2), breaks = c(0.01, 0.05, 0.1, 0.2), oob=squish, name="\n\nP-value (FDR-corrected)") +
          geom_hline(yintercept = 0.01, alpha=1) + 
          geom_hline(yintercept = 0.05, alpha=0.7) + 
          geom_hline(yintercept = 0.1, alpha=0.4) + 
          geom_vline(xintercept=seq(0,length(levels(df_toplot$Label)))+0.5 ,color="black", alpha=0.2) +
          facet_wrap(~Threshold, ncol = length(unique(df_toplot$Threshold))) +
          # scale_y_continuous(name = paste0("% abnormal voxels above Z=",levels(results$threshold)[i]), limits = c(0, quantile(df_tmp$perc_vox_above_thresh, 0.90))) +
          theme_classic() + 
          theme(text = element_text(size=15),
                plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(angle = 30, vjust = 1, hjust=1),
                panel.spacing.x = unit(1, "cm"))
    
    plot_betas = ggplot(df_toplot, aes(x=Label, y=Estimate_group, color=factor(Micro, levels=names), group = factor(Micro, levels=names))) + 
          # geom_point(position = position_jitter(w = 0.4, h = 0)) +
          geom_point(position = position_dodge(0.8)) +
          # geom_line(position = position_dodge(0.8), alpha=0.2) +
          scale_color_manual(name="Micro", values=color_scale) + 
          scale_x_discrete(labels = tissue_toplot, name="") +
          scale_y_continuous(limits = c(-1,1), breaks = seq(-1, 1, by=0.2), oob=squish, name="\n\nStandardized Beta") +
          geom_hline(yintercept = 0, alpha=1) + 
          geom_vline(xintercept=seq(0,length(levels(df_toplot$Label)))+0.5 ,color="black", alpha=0.2) +
          facet_wrap(~Threshold, ncol = length(unique(df_toplot$Threshold))) +
          # scale_y_continuous(name = paste0("% abnormal voxels above Z=",levels(results$threshold)[i]), limits = c(0, quantile(df_tmp$perc_vox_above_thresh, 0.90))) +
          ggtitle(paste0(dx_name, " (n = ",n_dx,")")) +
          theme_classic() + 
          theme(text = element_text(size=15),
                plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(angle = 30, vjust = 1, hjust=1),
                panel.spacing.x = unit(1, "cm"))
      
    wrap_plots(plot_betas, plot_pval, ncol=1)
    ggsave(name, width=6*(length(unique(df_toplot$Threshold))), height=10)
    print(name)
}


plot_beta_pval_nodbm = function(dx_name, icd_list, name) {
    color_scale = c(MD = "#04319E", ISOVF = "#8298CF",
                    FA = "#2B520B" , ICVF= "#448312", OD="#A2C189",
                    T2star="#D36108", QSM="#E9B084",
                    jacobians_rel="#5900cb", jacobians_abs="#893cef")
    n_dx = length(unique(df$ID[which(df$icd_code %in% icd_list)]))

    micro_list = c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")

    df_toplot = subset(lm_results, dx == dx_name & Label %in% c(3,4,5,6,7,10) & Micro %in% micro_list)

    df_toplot$fdr_group = p.adjust(df_toplot$pval_group, method="fdr")

    plot_pval = ggplot(df_toplot, aes(x=Label, y=fdr_group, color=factor(Micro, levels=micro_list), group = factor(Micro, levels=micro_list))) + 
          # geom_point(position = position_jitter(w = 0.4, h = 0)) +
          geom_point(position = position_dodge(0.8)) +
          # geom_line(position = position_dodge(0.8), alpha=0.2) +
          scale_color_manual(name="Micro", values=color_scale) + 
          scale_x_discrete(labels = tissue_toplot, name="") +
          scale_y_continuous(limits = c(0, 0.2), breaks = c(0.01, 0.05, 0.1, 0.2), oob=squish, name="\n\nP-value (FDR-corrected)") +
          geom_hline(yintercept = 0.01, alpha=1) + 
          geom_hline(yintercept = 0.05, alpha=0.7) + 
          geom_hline(yintercept = 0.1, alpha=0.4) + 
          geom_vline(xintercept=seq(0,length(levels(df_toplot$Label)))+0.5 ,color="black", alpha=0.2) +
          facet_wrap(~Threshold, ncol = length(unique(df_toplot$Threshold))) +
          # scale_y_continuous(name = paste0("% abnormal voxels above Z=",levels(results$threshold)[i]), limits = c(0, quantile(df_tmp$perc_vox_above_thresh, 0.90))) +
          theme_classic() + 
          theme(text = element_text(size=15),
                plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(angle = 30, vjust = 1, hjust=1),
                panel.spacing.x = unit(1, "cm"))
    
    plot_betas = ggplot(df_toplot, aes(x=Label, y=Estimate_group, color=factor(Micro, levels=micro_list), group = factor(Micro, levels=micro_list))) + 
          # geom_point(position = position_jitter(w = 0.4, h = 0)) +
          geom_point(position = position_dodge(0.8)) +
          # geom_line(position = position_dodge(0.8), alpha=0.2) +
          scale_color_manual(name="Micro", values=color_scale) + 
          scale_x_discrete(labels = tissue_toplot, name="") +
          scale_y_continuous(limits = c(-1,1), breaks = seq(-1, 1, by=0.2), oob=squish, name="\n\nStandardized Beta") +
          geom_hline(yintercept = 0, alpha=1) + 
          geom_vline(xintercept=seq(0,length(levels(df_toplot$Label)))+0.5 ,color="black", alpha=0.2) +
          facet_wrap(~Threshold, ncol = length(unique(df_toplot$Threshold))) +
          # scale_y_continuous(name = paste0("% abnormal voxels above Z=",levels(results$threshold)[i]), limits = c(0, quantile(df_tmp$perc_vox_above_thresh, 0.90))) +
          ggtitle(paste0(dx_name, " (n = ",n_dx,")")) +
          theme_classic() + 
          theme(text = element_text(size=15),
                plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(angle = 30, vjust = 1, hjust=1),
                panel.spacing.x = unit(1, "cm"))
      
    wrap_plots(plot_betas, plot_pval, ncol=1)
    ggsave(name, width=5*(length(unique(df_toplot$Threshold))), height=10)
    print(name)
}


for (d in 2:length(icd_codes_list)) {
  print(names(icd_codes_list)[d])

  try(plot_beta_pval(names(icd_codes_list)[d], icd_codes_list[[d]], paste0("./visualization/beta_pvals_",names(icd_codes_list)[d],".png")),silent=TRUE)
  try(plot_beta_pval_nodbm(names(icd_codes_list)[d], icd_codes_list[[d]], paste0("./visualization/beta_pvals_",names(icd_codes_list)[d],"_nodbm.png")), silent=TRUE)
}













# Radar plots for p-values: One plot by dx

dir.create("./visualization/radar_within_dx", showWarnings=FALSE)

df_radar_pval = function(df, name) {
  colnames(df)[ncol(df)] = "pval"

  color_scale = c(MD = "#04319E", ISOVF = "#8298CF", FA = "#2B520B" , ICVF= "#448312", OD="#A2C189", T2star="#D36108", QSM="#E9B084")
  n_dx = length(unique(results[which(results$dx == as.character(unique(df$dx))),1]))

  df_wide = df[,-which(colnames(df) %in% "dx")]
  df_wide = spread(df_wide, key=label_value, value=pval)
  colnames(df_wide) = c("micro", tissue_abn)
  df_wide = df_wide[,-which(colnames(df_wide) %in% "WMH")]
  row.names(df_wide) = df_wide$micro
  df_wide = df_wide[,-which(colnames(df_wide) %in% "micro")]
  df_wide = rbind(rep(1, ncol(df_wide)), rep(0, ncol(df_wide)), rep(0.5, ncol(df_wide)), df_wide)

  png(file=name, width=10, height=10, pointsize = 20, units = "in",res=300)
  par(xpd=TRUE)
  plot = radarchart(df_wide, axistype=1, maxmin=TRUE, seg=length(seq(0,1,0.05))-1,
                    pcol=c("black",color_scale),
                    pty=32, plwd=2, plty=c(0,rep(1, length(names))),
                    cglcol="grey", cglty=1, axislabcol="grey", caxislabels = c("0", rep("", length(seq(0,1,0.05))-3),"1"), cglwd=0.8,
                    vlcex = 1, palcex=0.8)
  dev.off()

  graphics::legend(x="bottom", y=NULL, legend=paste0("Cluster ",seq(1,4)), horiz=TRUE,
                    y.intersp = -0.5, text.width = 0.5,
                    bty="n", pch=20, col=color_scale_clust, text.col="black", cex=2, pt.cex=3)
  dev.off()

  print(name)

  return(plot)
}

for (d in levels(ttest_results$dx)) {
  print(d)
  for (t in levels(ttest_results$threshold)) {
    print(t)
    df_toplot = ttest_results %>% filter(dx == d, threshold == t)
    plot_pval = dx_point_pval(df_toplot %>% select(c("label_value", "micro", "dx", "pval")), paste0("./visualization/radar_within_dx/pval_",d,"_thresh_",t,".png"))
    plot_fdr = dx_point_pval(df_toplot %>% select(c("label_value", "micro","dx", "fdr")), paste0("./visualization/radar_within_dx/fdr_",d,"_thresh_",t,".png"))
  }
}

# Bar plots: one plot by micro and threshold (across dx and regions)

for (i in 1:length(levels(results$threshold))) {
  print(levels(results$threshold)[i])

  plots = list()
  for (n in 1:length(names)) {
    df_tmp = results %>% filter(threshold == levels(results$threshold)[i] & micro == names[n] & label_value != "9")
    dx_prev_tmp = as.data.frame(table(df_tmp$dx) / nrow(df_tmp %>% filter(ID == df_tmp$ID[1])))
    dx_prev_tmp = dx_prev_tmp[order(-dx_prev_tmp$Freq),]
    legend_labels = paste0(as.character(dx_prev_tmp$Var1), " (n=",dx_prev_tmp$Freq, ")")

    plots[[n]] = ggplot(results %>% filter(threshold == levels(results$threshold)[i] & micro == names[n] & label_value != "9"),
          aes(x=label_value, y=perc_vox_above_thresh, fill=factor(dx, levels=c(dx_prev_order)))) + 
          geom_boxplot(outlier.shape = NA) + 
          scale_fill_discrete(name="Diagnosis", labels=legend_labels) + 
          scale_x_discrete(labels = tissue_all[c(seq(3,8),10)], name="") + 
          scale_y_continuous(name = paste0("% abnormal voxels above Z=",levels(results$threshold)[i]), limits = c(0, quantile(df_tmp$perc_vox_above_thresh, 0.90))) +
          ggtitle(paste0(names[n])) +
          theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5))
    ggsave(paste0("./visualization/perc_abn_vox_Z",levels(results$threshold)[i],"_",names[n],".png"), width=20, height=5)
  }

  wrap_plots(plots, ncol=1)
  ggsave(paste0("./visualization/perc_abn_vox_Z",levels(results$threshold)[i],"_allmicro.png"), width=20, height=(5*length(names)))
}
