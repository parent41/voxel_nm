
library(data.table)
library(RMINC)

# Arg 1: Micro name
args = commandArgs(trailingOnly=TRUE)

# args = c()
# args[1] = "1"

dx_num = as.numeric(args[1])
names = c("FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM")

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

print(paste0("Diagnosis = ", colnames(dx)[dx_num]))
dx_ids = dx[,dx_num]
dx_ids = na.omit(dx_ids)
    
# For raw maps and denoised maps
for (n in 1:length(names)) {
    print(names[n])

    # zscores_raw = as.data.frame(fread(paste0("../subj_zscores_common_dx/results/zscores_",names[n],".tsv"), header = TRUE))
    zscores_anlm = as.data.frame(fread(paste0("../subj_zscores_common_dx/results/zscores_anlm_",names[n],".tsv"), header = TRUE))

    # zscores_raw_dx = as.data.table(zscores_raw[which(inclusions$ID %in% dx_ids),])
    zscores_anlm_dx = as.data.table(zscores_anlm[which(inclusions$ID %in% dx_ids),])

    # zscores_raw_dx_avg = zscores_raw_dx[, lapply(.SD, mean, na.rm = TRUE)]
    zscores_anlm_dx_avg = zscores_anlm_dx[, lapply(.SD, mean, na.rm = TRUE)]

    # fwrite(zscores_raw_dx_avg, paste0("./results/",colnames(dx)[dx_num],"_avg_",names[n],"_raw.tsv"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
    fwrite(zscores_anlm_dx_avg, paste0("./results/",colnames(dx)[dx_num],"_avg_",names[n],"_anlm.tsv"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

    # vol = as.numeric(zscores_raw_dx_avg)
    # vol[is.na(vol)] <- 0
    # outvol <- mincGetVolume("../../../UKB/temporary_template/avg.020_2mm.mnc")
    # outvol[] <- 0
    # outvol[mask > 0.5] <- vol
    # mincWriteVolume(outvol, paste0("./results/",colnames(dx)[dx_num],"_avg_",names[n],"_raw.mnc"), clobber=TRUE, like="../../../UKB/temporary_template/avg.020_2mm.mnc", verbose=FALSE)

    vol = as.numeric(zscores_anlm_dx_avg)
    vol[is.na(vol)] <- 0
    outvol <- mincGetVolume("../../../UKB/temporary_template/avg.020_2mm.mnc")
    outvol[] <- 0
    outvol[mask > 0.5] <- vol
    mincWriteVolume(outvol, paste0("./results/",colnames(dx)[dx_num],"_avg_",names[n],"_anlm.mnc"), clobber=TRUE, like="../../../UKB/temporary_template/avg.020_2mm.mnc", verbose=FALSE)
}


