
library(RMINC)
library(data.table)
library(matrixStats)

# Combine sub-matrices into one matrix with all subjects
loc2 = list.files("../WMH_micro_spatial/micro_matrices", pattern="*ses2_Label_whole_brain.tsv", full.names=TRUE)
loc3 = list.files("../WMH_micro_spatial/micro_matrices", pattern="*ses3_Label_whole_brain.tsv", full.names=TRUE)

mask=mincGetVolume("../UKB/temporary_template/Mask_2mm_dil2.mnc")

tissue=c('Ventricules', 'CSF', 'Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH')

prev2=list() # List of tissue types containing matrices of tissue prevalence (chunk * vertices)
prev3=list() # List of tissue types containing matrices of tissue prevalence (chunk * vertices)

# For every tissue type in BISON, make matrix
for(t in 1:length(tissue)) {
    print(tissue[t])
    prev2[[t]]=matrix(nrow=length(loc2), ncol=sum(mask[] == 1))
    prev3[[t]]=matrix(nrow=length(loc3), ncol=sum(mask[] == 1))
}

# Count tissue prevalence for each chunk (ses2)
for (i in 1:length(loc2)){
    print(loc2[i])

    # Read chunk
    df = as.matrix(fread(loc2[i]))
    # df = round(df, 1)

    #Add chunk to final list of matrices
    for(t in 1:length(tissue)) {
        print(tissue[t])
        prev2[[t]][i,]=colSums2(df == t, na.rm=TRUE)
    }
}

# Count tissue prevalence for each chunk (ses3)
for (i in 1:length(loc3)){
    print(loc3[i])

    # Read chunk
    df = as.matrix(fread(loc3[i]))
    # df = round(df, 1)

    #Add chunk to final list of matrices
    for(t in 1:length(tissue)) {
        print(tissue[t])
        prev3[[t]][i,]=colSums2(df == t, na.rm=TRUE)
    }
}

# Final prevalence for each tissue type by summing chunk prevalence

prev_final2 = matrix(nrow=length(tissue), ncol=sum(mask[] == 1))
prev_final3 = matrix(nrow=length(tissue), ncol=sum(mask[] == 1))
prev_final23 = matrix(nrow=length(tissue), ncol=sum(mask[] == 1))

for(t in 1:length(tissue)) {
    print(tissue[t])
    prev_final[t,] = colSums(prev[[t]], na.rm=TRUE)
}

mode(prev_final) = 'integer'
write.table(prev_final, "./results/tissue_prevalence.tsv", row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)

# Write minc files
prev_final = fread("./results/tissue_prevalence.tsv")

for(t in 1:length(tissue)) {
    outvol <- mincGetVolume("../../../UKB/temporary_template/avg.020_2mm.mnc")
    maskvol <- mincGetVolume("../../../UKB/temporary_template/Mask_2mm.mnc")
    outvol[] <- 0
    outvol[maskvol > 0.5] <- as.numeric(prev_final[t,])
    mincWriteVolume(outvol, paste0("./results/",t,"_",tissue[t],"_prevalence.mnc"), clobber=TRUE, like="../../../UKB/temporary_template/Mask_2mm.mnc")
}













