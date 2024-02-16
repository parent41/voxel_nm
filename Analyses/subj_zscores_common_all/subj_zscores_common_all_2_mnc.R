
library(data.table)
library(RMINC)
library(splitstackshape)

# Arg 1: Micro name
args = commandArgs(trailingOnly=TRUE)
num_micro = as.numeric(args[1]) + 1
num_chunk = as.numeric(args[2])

dir.create("./results/mnc")

# args = c()
# num_micro = 1
# num_chunk = 32

names = c("FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM", "jacobians_abs", "jacobians_rel")

# Load data

anatVol = mincArray(mincGetVolume("../../../UKB/temporary_template/avg.020_2mm.mnc"))
mask = mincGetVolume("../../../UKB/temporary_template/Mask_2mm_dil2.mnc")

zscores = as.data.frame(fread(paste0("./results/zscores_c",num_chunk,"_",names[num_micro],"_anlm.tsv"), header=TRUE))
zscores_ids = as.data.frame(fread(paste0("./tmp/ids_micro_c",num_chunk,"_",names[num_micro],"_anlm.txt")))
zscores_ids = as.numeric(zscores_ids$V1)

outvol <- mincGetVolume("../../../UKB/temporary_template/avg.020_2mm.mnc")
for (i in 1:nrow(zscores)) {
    print(zscores_ids[i])
    
    time = system.time({
    vol = as.numeric(zscores[i,])
    vol[is.na(vol)] <- 0
    outvol[] <- 0
    outvol[mask > 0.5] <- vol
    mincWriteVolume(outvol, paste0("./results/mnc/",zscores_ids[i],"_zscore_anlm_",names[num_micro],".mnc"), clobber=TRUE, like="../../../UKB/temporary_template/avg.020_2mm.mnc", verbose=FALSE)
    })[3]
    print(time)
}


