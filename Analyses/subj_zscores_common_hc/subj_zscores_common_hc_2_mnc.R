
library(data.table)
library(RMINC)

# Arg 1: Micro name
args = commandArgs(trailingOnly=TRUE)

# args = c()
# args[1] = "FA"
# args[2] = 32

# names = c("FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM")

# Load data

inclusions = as.data.frame(fread("../../../WMH_micro_spatial/QC/inclusions_excluding_dx_new.txt"))
colnames(inclusions) = "ID"

demo = as.data.frame(fread("../../../UKB/tabular/df_demo/UKBB_demo_wider.tsv"))
demo = subset(demo, InstanceID == 2, select=c('SubjectID', 'Sex_31_0', 'Age_when_attended_assessment_centre_21003_0'))
colnames(demo) = c("ID", "Sex", "Age")
demo = merge(inclusions, demo, by="ID", all=FALSE)

anatVol = mincArray(mincGetVolume("../../../UKB/temporary_template/avg.020_2mm.mnc"))
mask = mincGetVolume("../../../UKB/temporary_template/Mask_2mm_dil2.mnc")

# For raw maps and denoised maps

ids_raw = as.data.frame(fread(paste0("./tmp/ids_micro_c",args[2],"_",args[1],".txt")))
ids_anlm = as.data.frame(fread(paste0("./tmp/ids_micro_c",args[2],"_",args[1],"_anlm.txt")))

print(paste0("Are ID lists identical between raw and anlm? ",all(ids_raw$V1 == ids_anlm$V1)))

zscores_raw = as.data.frame(fread(paste0("./results/zscores_c",args[2],"_",args[1],".tsv"), header = TRUE))
zscores_anlm = as.data.frame(fread(paste0("./results/zscores_c",args[2],"_",args[1],"_anlm.tsv"), header = TRUE))

# Make mnc files for each subject grouped by dx

ids = as.numeric(ids_raw$V1)

dir.create("./results/mnc", showWarnings=FALSE)

for (i in 1:length(ids)) {
    id_row = which(inclusions$ID == ids[i])

    # If dx ID is in inclusions
    if(length(id_row)>0) {

        res_dir = paste0("./results/",colnames(dx)[d], "/", dx_ids[i])
        dir.create(res_dir, showWarnings=FALSE)

        # Write to mnc (raw)
        vol = as.numeric(zscores_raw[i,])
        vol[is.na(vol)] <- 0
        outvol <- mincGetVolume("../../../UKB/temporary_template/avg.020_2mm.mnc")
        outvol[] <- 0
        outvol[mask > 0.5] <- vol
        mincWriteVolume(outvol, paste0("./results/mnc/sub-",ids[i],"_ses-2_",args[1],"_zscore.mnc"), clobber=TRUE, like="../../../UKB/temporary_template/avg.020_2mm.mnc", verbose=FALSE)

        # cat(paste0("\n\t\tMicro = ", res_dir,"/",demo$ID[id_row],"_",args[1],".mnc \t"))

        # Write to mnc (anlm)
        vol = as.numeric(zscores_anlm[i,])
        vol[is.na(vol)] <- 0
        outvol <- mincGetVolume("../../../UKB/temporary_template/avg.020_2mm.mnc")
        outvol[] <- 0
        outvol[mask > 0.5] <- vol
        mincWriteVolume(outvol, paste0("./results/mnc/sub-",ids[i],"_ses-2_",args[1],"_zscore_anlm.mnc"), clobber=TRUE, like="../../../UKB/temporary_template/avg.020_2mm.mnc", verbose=FALSE)

        # cat(paste0("\n\t\tMicro = ", res_dir,"/",demo$ID[id_row],"_anlm_",args[1],".mnc \t"))
    }
}



