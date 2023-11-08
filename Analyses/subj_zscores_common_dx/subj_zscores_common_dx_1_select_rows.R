

library(RMINC)
library(data.table)
library(matrixStats)

# Arg1: Name of micro
args = commandArgs(trailingOnly=TRUE)

# args=c("OD")

mask = mincGetVolume("../../../UKB/temporary_template/Mask_2mm_dil2.mnc")

tissue_all=c('Ventricules', 'CSF', 'Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH')
tissue_nm=c('Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM')
tissue_brain = c('Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH')
# names = c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM")

print(args[1])

# Load demographics
inclusions = as.data.frame(fread("../../../WMH_micro_spatial/QC/inclusions_only_dx_new.txt"))
colnames(inclusions) = c("ID")

demo = as.data.frame(fread("../../../UKB/tabular/df_demo/UKBB_demo_wider.tsv"))
demo = subset(demo, InstanceID == 2, select=c('SubjectID', 'Sex_31_0', 'Age_when_attended_assessment_centre_21003_0'))
colnames(demo) = c("ID", "Sex", "Age")
demo = merge(inclusions, demo, by="ID", all.x=TRUE)

dx = as.data.frame(fread("../../../UKB/QC/exclusion_lists/dx_single.csv"))

# Load NM
nm=list()
for (t in 1:length(tissue_nm)) {
    print(tissue_nm[t])
    nm[[t]] = as.data.frame(fread(paste0("../norm_models_BLR/results/nm_",args[1],"_",tissue_nm[t],".tsv")))
}

# Load label

load_only = "NO"

ids_label = as.data.frame(fread("../../bison_matrices/ids_ses2_Label_whole_brain.txt"))
colnames(ids_label) = c("ID")

indices_label = which(ids_label$ID %in% inclusions$ID)

if (load_only == "NO") {
    fwrite(as.data.frame(indices_label), "./tmp/idx_label.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

    # Run script to select certain rows from the huge label tsv
    command="./select_rows.sh ../../bison_matrices/ses2_Label_whole_brain.tsv ./results/ses2_Label_whole_brain_dx.tsv ./tmp/idx_label.txt"
    system(command)
}

# label_ids = as.data.frame(fread("./results/ses2_Label_whole_brain_dx.tsv"))
# print(dim(label_ids))

# Load micro

ids_micro = as.data.frame(fread(paste0("../../micro_matrices/ids_ses2_",args[1],".txt")))
colnames(ids_micro) = c("ID")

indices_micro = which(ids_micro$ID %in% inclusions$ID)
print(length(indices_micro))

if (load_only == "NO") {
    fwrite(as.data.frame(indices_micro), paste0("./tmp/idx_",args[1],".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)

    command=paste0("./select_rows.sh ../../micro_matrices/ses2_",args[1],".tsv ./results/ses2_",args[1],"_dx.tsv ./tmp/idx_",args[1],".txt")
    system(command)
}

# micro_ids = as.data.frame(fread(paste0("./results/ses2_",args[1],"_dx.tsv")))
# print(dim(micro_ids))

