
library(data.table)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reticulate)
library(RMINC)
library(MRIcrotome)
library(viridis)
library(RColorBrewer)
library(foreach)
library(doParallel)

# Load data

names = c("FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM", "jacobians_abs", "jacobians_rel")
# names = c("FA")
# names = c("jacobians_abs", "jacobians_rel")

tissues = c('Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM')
tissues_all = c('Ventricules', 'CSF', 'Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH')

mask = mincGetVolume("../../../UKB/temporary_template/Mask_2mm_dil2.mnc")

nsubj = as.data.frame(fread("../../../WMH_micro_spatial/QC/inclusions_excluding_dx_new.txt"))
nsubj = nrow(nsubj)

# PNGs of voxel-wise NM outputs

raw_up_thresh = c(0.85, 0.003, 0.85, 0.5, 1, 0.5, 150, 1, 1)
raw_down_thresh = c(0, 0, 0, 0, 0, 0, -100, -1, -1)
# raw_up_thresh = c(1, 1) # For DBM
# raw_down_thresh = c(-1, -1) # For DBM


nm = as.data.frame(fread(paste0("./results/nm_",names[1],"_",tissues[1],".tsv")))
nm = colnames(nm)

cores = detectCores()
parallelCluster <- makeCluster(cores - 5, type = "SOCK", methods = FALSE)
setDefaultCluster(parallelCluster)
registerDoParallel(parallelCluster)

dir.create("./visualization/mnc", showWarnings=FALSE)
for (n in 1:length(names)) {
# foreach(n=1:length(names), .packages = c('RMINC', 'viridis', 'MRIcrotome', 'tidyverse')) %dopar% {
    print(names[n])
    dir.create(paste0("./visualization/mnc/",names[n]), showWarnings=FALSE)
    
    for (t in 1:length(tissues)) {
        print(tissues[t])
        dir.create(paste0("./visualization/mnc/",names[n],"/",tissues[t]), showWarnings=FALSE)

        # for (c in 1:length(nm)) {
        foreach(c=1:length(nm), .packages = c('RMINC', 'viridis', 'MRIcrotome', 'tidyverse')) %dopar% {
            print(nm[c])
            anatVol = mincArray(mincGetVolume("../../../UKB/temporary_template/avg.020_2mm.mnc"))
            overlayVol = mincArray(mincGetVolume(paste0("./results/mnc/",names[n],"/",tissues[t],"/",nm[c],".mnc")))

            if (grepl("male", nm[c], ignore.case = TRUE)) {

                # Clamp values
                overlayVol[] = pmin(pmax(overlayVol[], raw_down_thresh[n]), raw_up_thresh[n])

                up_thresh = raw_up_thresh[n]
                down_thresh = raw_down_thresh[n]

                # Positive only
                if (names[n] %in% c("QSM", "jacobians_abs", "jacobians_rel")) {
                    overlayVol[] = overlayVol[] + abs(raw_down_thresh[n])
                    down_thresh = down_thresh + abs(raw_down_thresh[n])
                    up_thresh = up_thresh + abs(raw_down_thresh[n])
                }

                png(file=paste0("./visualization/mnc/",names[n],"/",tissues[t], "/",names[n],"_",tissues[t], "_",nm[c],".png"), width=8500, height=4000, pointsize = 150)
                sliceSeries(nrow=1, ncol=1, begin = 51, end = 51, dimension=1) %>%
                    anatomy(anatVol, low=10, high=200) %>%
                    overlay(overlayVol, low=down_thresh, high=up_thresh,col=magma(255)) %>%
                sliceSeries(nrow=1, ncol=1, begin = 59, end = 59, dimension=2) %>%
                    anatomy(anatVol, low=10, high=200) %>%
                    overlay(overlayVol, low=down_thresh, high=up_thresh,col=magma(255)) %>%
                sliceSeries(nrow=1, ncol=1, begin = 40, end = 40, dimension=3) %>%
                    anatomy(anatVol, low=10, high=200) %>%
                    overlay(overlayVol, low=down_thresh, high=up_thresh,col=magma(255)) %>%
                legend(nm[c]) %>%
                draw()
                dev.off()
            } else if (nm[c] == "N") {
                png(file=paste0("./visualization/mnc/",names[n],"/",tissues[t], "/",names[n],"_",tissues[t], "_",nm[c],".png"), width=8500, height=4000, pointsize = 150)
                sliceSeries(nrow=1, ncol=1, begin = 51, end = 51, dimension=1) %>%
                    anatomy(anatVol, low=10, high=200) %>%
                    overlay(overlayVol, low=100, high=nsubj, col=magma(255)) %>%
                sliceSeries(nrow=1, ncol=1, begin = 59, end = 59, dimension=2) %>%
                    anatomy(anatVol, low=10, high=200) %>%
                    overlay(overlayVol, low=100, high=nsubj, col=magma(255)) %>%
                sliceSeries(nrow=1, ncol=1, begin = 40, end = 40, dimension=3) %>%
                    anatomy(anatVol, low=10, high=200) %>%
                    overlay(overlayVol, low=100, high=nsubj, col=magma(255)) %>%
                legend(nm[c]) %>%
                draw()
                dev.off()
            } else {
                png(file=paste0("./visualization/mnc/",names[n],"/",tissues[t], "/",names[n],"_",tissues[t], "_",nm[c],".png"), width=8500, height=4000, pointsize = 150)
                sliceSeries(nrow=1, ncol=1, begin = 51, end = 51, dimension=1) %>%
                    anatomy(anatVol, low=10, high=200) %>%
                    overlay(overlayVol, low=min(overlayVol[][overlayVol[]!=0]), high=max(overlayVol[][overlayVol[]!=0]),col=magma(255)) %>%
                sliceSeries(nrow=1, ncol=1, begin = 59, end = 59, dimension=2) %>%
                    anatomy(anatVol, low=10, high=200) %>%
                    overlay(overlayVol, low=min(overlayVol[][overlayVol[]!=0]), high=max(overlayVol[][overlayVol[]!=0]),col=magma(255)) %>%
                sliceSeries(nrow=1, ncol=1, begin = 40, end = 40, dimension=3) %>%
                    anatomy(anatVol, low=10, high=200) %>%
                    overlay(overlayVol, low=min(overlayVol[][overlayVol[]!=0]), high=max(overlayVol[][overlayVol[]!=0]),col=magma(255)) %>%
                legend(nm[c]) %>%
                draw()
                dev.off()
            }
        }
    }
}

stopCluster(parallelCluster)


