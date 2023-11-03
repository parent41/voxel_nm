
library(RMINC)
library(data.table)

# Arg1 = location of filelist
# Arg2 = mask
# Arg3 = chunk_micro
args = commandArgs(trailingOnly=TRUE)

# args=c()
# args[1] = "./micro_file_lists/sub6_ses2_dti_FA.txt"
# args[2] = "../UKB/temporary_template/Mask_2mm_dil2.mnc"
# args[3] = "sub6_ses2_FA"

print(args)

mask=mincGetVolume(args[2])

# Load micro maps of chunk of subjects
fl=read.table(args[1])
fl=c(fl$V1)

ids = sub('.*sub-(\\d+)_ses.*', '\\1', fl)
fwrite(as.data.frame(ids), paste0("./micro_matrices/ids_",args[3],".txt"), row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE ) # List of subjects included

mat = matrix(nrow=length(fl), ncol=sum(mask[] == 1))

for (i in 1:length(fl)){
    print(fl[i])
    mat[i,]=mincGetVolume(fl[i])[mask>0.5]
}

fwrite(mat, paste0("./micro_matrices/",args[3],".tsv"), row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)

