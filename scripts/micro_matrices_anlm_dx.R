
library(RMINC)
library(data.table)

# Arg1 = Long name of micro
# Arg2 = Short name of micro
args = commandArgs(trailingOnly=TRUE)

# args=c()
# args[1] = "dti_FA"
# args[2] = "FA"

print(args)

mask=mincGetVolume("../UKB/temporary_template/Mask_2mm_dil2.mnc")

# Load micro maps of chunk of subjects
# fl = list.files("./maps_UKB_space_anlm_dx", pattern=paste0("*",args[1],"_UKB_anlm*"), full.names=TRUE)
fl = list.files("./maps_UKB_space_anlm_dx", pattern=paste0("*",args[1],"_2mm_anlm*"), full.names=TRUE)

ids = sub('.*sub-(\\d+)_ses.*', '\\1', fl)
fwrite(as.data.frame(ids), paste0("./micro_matrices/ids_dx_anlm_",args[2],".txt"), row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE ) # List of subjects included

mat = matrix(nrow=length(fl), ncol=sum(mask[] == 1))

for (i in 1:length(fl)){
    print(fl[i])
    mat[i,]=mincGetVolume(fl[i])[mask>0.5]
}

fwrite(mat, paste0("./micro_matrices/anlm_dx_ses2_",args[2],".tsv"), row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)

