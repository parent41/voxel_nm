
library(data.table)
library(RMINC)

# Load data

names = c("FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM", "jacobians_abs", "jacobians_rel")
# names = c("FA")
# names = c("jacobians_abs", "jacobians_rel")

tissues = c('Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM')
tissues_all = c('Ventricules', 'CSF', 'Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH')

mask = mincGetVolume("../../../UKB/temporary_template/Mask_2mm_dil2.mnc")

# Export NM results to mnc
dir.create("./results/mnc", showWarnings=FALSE)
for (n in 1:length(names)) {
    print(names[n])
    dir.create(paste0("./results/mnc/",names[n]), showWarnings=FALSE)

    for (t in 1:length(tissues)) {
        print(tissues[t])
        dir.create(paste0("./results/mnc/",names[n],"/",tissues[t]), showWarnings=FALSE)

        nm = as.data.frame(fread(paste0("./results/nm_",names[n],"_",tissues[t],".tsv")))

        for (c in 1:ncol(nm)) {
            print(colnames(nm)[c])
            outvol <- mincGetVolume("../../../UKB/temporary_template/avg.020_2mm.mnc")
            outvol[] <- 0
            col = nm[,c]
            col[is.na(col)] = 0
            outvol[mask > 0.5] = col
            mincWriteVolume(outvol, paste0("./results/mnc/",names[n],"/",tissues[t],"/",colnames(nm)[c],".mnc"), clobber=TRUE)
        }
    }
}


