
library(data.table)
library(RMINC)
library(splitstackshape)

# Arg 1: Micro name
args = commandArgs(trailingOnly=TRUE)

# args = c()
# args[1] = "FA"

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

# Select random subset to visualize

set.seed(123)
sampled_df = stratified(demo, c("Sex", "Age"), size = 10)
sampled_df = sampled_df[order(sampled_df$ID), ]

# Select rows in huge voxel matrix and load data

# Label
ids_label = as.data.frame(fread("../../bison_matrices/ids_ses2_Label_whole_brain.txt"))
indices_label = which(ids_label$V1 %in% sampled_df$ID)
fwrite(as.data.frame(indices_label), paste0("./tmp/idx_sampled_label_",args[1],".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
# Test with ID file (make sure same IDs as in chunk)
command=paste0("./select_rows.sh ../../bison_matrices/ids_ses2_Label_whole_brain.txt ./tmp/ids_sampled_label_",args[1],".txt ./tmp/idx_sampled_label_",args[1],".txt")
system(command)
# Make sub-matrices of subjects x voxels
command=paste0("./select_rows.sh ../../bison_matrices/ses2_Label_whole_brain.tsv ./tmp/sampled_label_",args[1],".tsv ./tmp/idx_sampled_label_",args[1],".txt")
system(command)
# Load matrix
label=as.data.frame(fread(paste0("./tmp/sampled_label_",args[1],".tsv")))
ids_label = as.data.frame(fread(paste0("./tmp/ids_sampled_label_",args[1],".txt")))

# Z-scores
zscores_raw_list = list()
zscores_anlm_list = list()
ids_micro_list = list()

nchunks = 32

for(i in 0:nchunks) {
    print(i)
    ids_micro_chunk = as.data.frame(fread(paste0("./tmp/ids_micro_c",i,"_",args[1],".txt")))
    indices_micro = which(ids_micro_chunk$V1 %in% sampled_df$ID)
    fwrite(as.data.frame(indices_micro), paste0("./tmp/idx_sampled_micro_c",i,"_",args[1],".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
    # Test with ID file (make sure same IDs as in chunk)
    command=paste0("./select_rows.sh ./tmp/ids_micro_c",i,"_",args[1],".txt ./tmp/ids_sampled_micro_c",i,"_",args[1],".txt ./tmp/idx_sampled_micro_c",i,"_",args[1],".txt")
    system(command)
    # Make sub-matrices of subjects x voxels
    command=paste0("./select_rows.sh ./results/zscores_c",i,"_",args[1],".tsv ./tmp/sampled_label_c",i,"_",args[1],".tsv ./tmp/idx_sampled_micro_c",i,"_",args[1],".txt")
    system(command)
    command=paste0("./select_rows.sh ./results/zscores_c",i,"_",args[1],"_anlm.tsv ./tmp/sampled_label_c",i,"_",args[1],"_anlm.tsv ./tmp/idx_sampled_micro_c",i,"_",args[1],".txt")
    system(command)
    # Load matrices
    zscores_raw_list[[i+1]] = as.data.frame(fread(paste0("./tmp/sampled_label_c",i,"_",args[1],".tsv")))
    zscores_anlm_list[[i+1]] = as.data.frame(fread(paste0("./tmp/sampled_label_c",i,"_",args[1],"_anlm.tsv")))
    ids_micro_list[[i+1]] = as.data.frame(fread(paste0("./tmp/ids_sampled_micro_c",i,"_",args[1],".txt")))
}

zscores_raw <- do.call(rbind, zscores_raw_list)
zscores_anlm <- do.call(rbind, zscores_anlm_list)
ids_micro = do.call(rbind, ids_micro_list)

print(dim(label))
print(dim(ids_label))
print(dim(zscores_raw))
print(dim(zscores_anlm))
print(dim(ids_micro))

print(paste0("Are ID lists identical between label and micro? ",all(ids_label$V1 == ids_micro$V1)))

fwrite(ids_micro, paste0("./tmp/ids_sampled_micro_",args[1],".txt"), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
fwrite(ids_label, paste0("./tmp/ids_sampled_label_",args[1],".txt"), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

# Make mnc files for each subject

ids = as.numeric(ids_micro$V1)

dir.create("./results/mnc", showWarnings=FALSE)

for (i in 1:length(ids)) {
    print(ids[i])
    demo_id = demo[which(demo$ID == ids[i]),]

    # Write to mnc (raw)
    vol = as.numeric(zscores_raw[i,])
    vol[is.na(vol)] <- 0
    outvol <- mincGetVolume("../../../UKB/temporary_template/avg.020_2mm.mnc")
    outvol[] <- 0
    outvol[mask > 0.5] <- vol
    mincWriteVolume(outvol, paste0("./results/mnc/",ids[i],"_",demo_id$Sex,"_",demo_id$Age,"_zscore_",args[1],".mnc"), clobber=TRUE, like="../../../UKB/temporary_template/avg.020_2mm.mnc", verbose=FALSE)

    # cat(paste0("\n\t\tMicro = ", res_dir,"/",demo$ID[id_row],"_",args[1],".mnc \t"))

    # Write to mnc (anlm)
    vol = as.numeric(zscores_anlm[i,])
    vol[is.na(vol)] <- 0
    outvol <- mincGetVolume("../../../UKB/temporary_template/avg.020_2mm.mnc")
    outvol[] <- 0
    outvol[mask > 0.5] <- vol
    mincWriteVolume(outvol, paste0("./results/mnc/",ids[i],"_",demo_id$Sex,"_",demo_id$Age,"_zscore_anlm_",args[1],".mnc"), clobber=TRUE, like="../../../UKB/temporary_template/avg.020_2mm.mnc", verbose=FALSE)

    # cat(paste0("\n\t\tMicro = ", res_dir,"/",demo$ID[id_row],"_anlm_",args[1],".mnc \t"))
}



