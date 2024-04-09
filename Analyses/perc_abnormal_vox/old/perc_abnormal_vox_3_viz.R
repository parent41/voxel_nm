
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
library(fmsb)
library(ggradar)

tissue_nm = c('Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM')
tissue_all=c('Ventricules', 'CSF', 'Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH', 'Cerebral_WM')
tissue_abn = c('Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH', 'Cerebral_WM')

names = c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")

# Load results
results_dx = as.data.frame(fread("./results/final_results_dx.tsv"))

results_dx[,c(2,3,7,8)] = as.data.frame(lapply(results_dx[,c(2,3,7,8)], as.factor))

dx_prev_order = as.data.frame(table(results_dx$dx))
dx_prev_order = dx_prev_order[order(-dx_prev_order$Freq),]
dx_prev_order$Freq = dx_prev_order$Freq / 147
dx_prev_order = subset(dx_prev_order, Freq >= 30, select="Var1")
dx_prev_order = as.character(dx_prev_order$Var1)

# One-sided t-test of dx > controls

ttest_results = as.data.frame(matrix(ncol=6, nrow=0))
colnames(ttest_results) = c("label_value", "threshold", "micro", "dx", "tval", "pval")

comparison_levels <- setdiff(levels(results_dx$dx), "HC")

c=1
for (l in levels(results_dx$label_value)) {
  print(l)
  for (t in levels(results_dx$threshold)) {
    print(t)
    for (m in levels(results_dx$micro)) {
      print(m)
      for (d in comparison_levels) {
        print(d)
        ttest_df = subset(results_dx, label_value == l & threshold == t & micro == m & dx %in% c(d, "HC"))
        ttest = t.test(perc_vox_above_thresh ~ dx, data = ttest_df, alternative="less")
        ttest_results[c,] = c(l, t, m, d, ttest$statistic, ttest$p.value)

        c = c+1    
      }
    }
  }
}


ttest_results[,c(1,2,3,4)] = as.data.frame(lapply(ttest_results[,c(1,2,3,4)], as.factor))
ttest_results[,c(5,6)] = as.data.frame(lapply(ttest_results[,c(5,6)], as.numeric))

ttest_results$fdr = p.adjust(ttest_results$pval, method="fdr")

fwrite(ttest_results, "./results/ttest_results.tsv", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

ttest_results = as.data.frame(fread("./results/ttest_results.tsv"))
ttest_results[,c(1,2,3,4)] = as.data.frame(lapply(ttest_results[,c(1,2,3,4)], as.factor))
ttest_results[,c(5,6,7)] = as.data.frame(lapply(ttest_results[,c(5,6,7)], as.numeric))

# Point plots for p-values: One plot by dx

dir.create("./visualization/point_within_dx", showWarnings=FALSE)

dx_point_pval = function(df, name) {
  colnames(df)[ncol(df)] = "pval"

  color_scale = c(MD = "#04319E", ISOVF = "#8298CF", FA = "#2B520B" , ICVF= "#448312", OD="#A2C189", T2star="#D36108", QSM="#E9B084")
  n_dx = length(unique(results_dx[which(results_dx$dx == as.character(unique(df$dx))),1]))

  plot = ggplot(df %>% filter(label_value != "9"),
          aes(x=label_value, y=pval, color=factor(micro, levels=names))) + 
          geom_point(position = position_jitter(w = 0.2, h = 0)) +
          scale_color_manual(name="Micro", values=color_scale) + 
          scale_x_discrete(labels = tissue_all[c(seq(3,8),10)], name="") +
          scale_y_continuous(limits = c(0, 0.2), breaks = seq(0, 0.2, by=0.01), oob=squish) +
          geom_hline(yintercept = 0.01, alpha=1) + 
          geom_hline(yintercept = 0.05, alpha=0.7) + 
          geom_hline(yintercept = 0.1, alpha=0.4) + 
          # scale_y_continuous(name = paste0("% abnormal voxels above Z=",levels(results_dx$threshold)[i]), limits = c(0, quantile(df_tmp$perc_vox_above_thresh, 0.90))) +
          ggtitle(paste0(as.character(unique(df$dx)), " (n = ",n_dx,")")) +
          theme(text = element_text(size=15),
                plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
    ggsave(name, width=5, height=5)
    print(name)

    return(plot)
}

for (d in levels(ttest_results$dx)) {
  print(d)
  for (t in levels(ttest_results$threshold)) {
    print(t)
    df_toplot = ttest_results %>% filter(dx == d, threshold == t)
    plot_pval = dx_point_pval(df_toplot %>% select(c("label_value", "micro", "dx", "pval")), paste0("./visualization/point_within_dx/pval_",d,"_thresh_",t,".png"))
    plot_fdr = dx_point_pval(df_toplot %>% select(c("label_value", "micro","dx", "fdr")), paste0("./visualization/point_within_dx/fdr_",d,"_thresh_",t,".png"))
  }
}

# Radar plots for p-values: One plot by dx

dir.create("./visualization/radar_within_dx", showWarnings=FALSE)

df_radar_pval = function(df, name) {
  colnames(df)[ncol(df)] = "pval"

  color_scale = c(MD = "#04319E", ISOVF = "#8298CF", FA = "#2B520B" , ICVF= "#448312", OD="#A2C189", T2star="#D36108", QSM="#E9B084")
  n_dx = length(unique(results_dx[which(results_dx$dx == as.character(unique(df$dx))),1]))

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
          scale_x_discrete(labels = tissue_all[c(seq(3,8),10)], name="") + 
          scale_y_continuous(name = paste0("% abnormal voxels above Z=",levels(results_dx$threshold)[i]), limits = c(0, quantile(df_tmp$perc_vox_above_thresh, 0.90))) +
          ggtitle(paste0(names[n])) +
          theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5))
    ggsave(paste0("./visualization/perc_abn_vox_Z",levels(results_dx$threshold)[i],"_",names[n],".png"), width=20, height=5)
  }

  wrap_plots(plots, ncol=1)
  ggsave(paste0("./visualization/perc_abn_vox_Z",levels(results_dx$threshold)[i],"_allmicro.png"), width=20, height=(5*length(names)))
}
