
library(RMINC)
library(data.table)

# Arg1 = location of filelist
# Arg2 = mask
# Arg3 = outfile
args = commandArgs(trailingOnly=TRUE)
print(args)

# args=c()
# args[1] = "./micro_file_lists/sub6_ses2_Label.txt"
# args[2] = "../UKB/temporary_template/Mask_2mm_dil2.mnc"
# args[3] = "./bison_matrices/sub6_ses2_Label_whole_brain.tsv"

# Load mask
mask=mincGetVolume(args[2])

fl=read.table(args[1])
fl=c(fl$V1)

ids = sub('.*sub-(\\d+)_ses.*', '\\1', fl)
ids_filename = sub('.*/([^/]+)\\..*', '\\1', args[3])
fwrite(as.data.frame(ids), paste0("./bison_matrices/ids_",ids_filename,".txt"), row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE ) # List of subjects included

mat = matrix(nrow=length(fl), ncol=sum(mask[] == 1))

for (i in 1:length(fl)){
    print(fl[i])
    mat[i,]=mincGetVolume(fl[i])[mask>0.5]
}

mat = round(mat,1)
mode(mat) = 'integer'

fwrite(mat, args[3], row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)

