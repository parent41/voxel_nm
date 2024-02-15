
library(data.table)
library(RMINC)

# Arg 1: Micro name
args = commandArgs(trailingOnly=TRUE)

# args = c()
# args[1] = "FA"

# names = c("FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM")

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

# For raw maps and denoised maps

zscores_raw = as.data.frame(fread(paste0("./results/zscores_",args[1],".tsv"), header = TRUE))
zscores_anlm = as.data.frame(fread(paste0("./results/zscores_anlm_",args[1],".tsv"), header = TRUE))

# Make mnc files for each subject grouped by dx

for (d in 1:ncol(dx)) {
    cat(paste0("\nDiagnosis = ", colnames(dx)[d]))
    dx_ids = dx[,d]
    dx_ids = na.omit(dx_ids)

    dir.create(paste0("./results/",colnames(dx)[d]), showWarnings=FALSE)

    for (i in 1:length(dx_ids)) {
        cat(paste0("\n\tID = ", dx_ids[i], "\t"))
        id_row = which(inclusions$ID == dx_ids[i])

        # If dx ID is in inclusions
        if(length(id_row)>0) {

            res_dir = paste0("./results/",colnames(dx)[d], "/", dx_ids[i])
            dir.create(res_dir, showWarnings=FALSE)

            # Write to mnc (raw)
            vol = as.numeric(zscores_raw[id_row,])
            vol[is.na(vol)] <- 0
            outvol <- mincGetVolume("../../../UKB/temporary_template/avg.020_2mm.mnc")
            outvol[] <- 0
            outvol[mask > 0.5] <- vol
            mincWriteVolume(outvol, paste0(res_dir,"/",demo$ID[id_row],"_",args[1],".mnc"), clobber=TRUE, like="../../../UKB/temporary_template/avg.020_2mm.mnc", verbose=FALSE)

            cat(paste0("\n\t\tMicro = ", res_dir,"/",demo$ID[id_row],"_",args[1],".mnc \t"))

            Write to mnc (anlm)
            vol = as.numeric(zscores_anlm[id_row,])
            vol[is.na(vol)] <- 0
            outvol <- mincGetVolume("../../../UKB/temporary_template/avg.020_2mm.mnc")
            outvol[] <- 0
            outvol[mask > 0.5] <- vol
            mincWriteVolume(outvol, paste0(res_dir,"/",demo$ID[id_row],"_anlm_",args[1],".mnc"), clobber=TRUE, like="../../../UKB/temporary_template/avg.020_2mm.mnc", verbose=FALSE)

            cat(paste0("\n\t\tMicro = ", res_dir,"/",demo$ID[id_row],"_anlm_",args[1],".mnc \t"))
        }
    }
}



