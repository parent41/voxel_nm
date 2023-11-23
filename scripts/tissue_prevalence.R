
library(RMINC)
library(data.table)
library(matrixStats)

# Combine sub-matrices into one matrix with all subjects
loc2 = list.files("./bison_matrices", pattern="*ses2_Label_whole_brain.tsv", full.names=TRUE)
loc3 = list.files("./bison_matrices", pattern="*ses3_Label_whole_brain.tsv", full.names=TRUE)

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
    prev_final2[t,] = colSums(prev2[[t]], na.rm=TRUE)
    prev_final3[t,] = colSums(prev3[[t]], na.rm=TRUE)
    prev_final23[t,] = prev_final2[t,] + prev_final3[t,]
}


mode(prev_final2) = 'integer'
mode(prev_final3) = 'integer'
mode(prev_final23) = 'integer'

fwrite(prev_final2, "./tissue_prevalence/tissue_prevalence_ses2.tsv", row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)
fwrite(prev_final3, "./tissue_prevalence/tissue_prevalence_ses3.tsv", row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)
fwrite(prev_final23, "./tissue_prevalence/tissue_prevalence_ses23.tsv", row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)

# Write minc files
prev_final2 = as.data.frame(fread("./tissue_prevalence/tissue_prevalence_ses2.tsv"))
prev_final3 = as.data.frame(fread("./tissue_prevalence/tissue_prevalence_ses3.tsv"))
prev_final23 = as.data.frame(fread("./tissue_prevalence/tissue_prevalence_ses23.tsv"))

for(t in 1:length(tissue)) {
    outvol <- mincGetVolume("../UKB/temporary_template/avg.020_2mm.mnc")
    maskvol <- mincGetVolume("../UKB/temporary_template/Mask_2mm_dil2.mnc")
    outvol[] <- 0
    outvol[maskvol > 0.5] <- as.numeric(prev_final2[t,])
    mincWriteVolume(outvol, paste0("./tissue_prevalence/",t,"_",tissue[t],"_prevalence_ses2.mnc"), clobber=TRUE, like="../../../UKB/temporary_template/Mask_2mm_dil2.mnc")
}

for(t in 1:length(tissue)) {
    outvol <- mincGetVolume("../UKB/temporary_template/avg.020_2mm.mnc")
    maskvol <- mincGetVolume("../UKB/temporary_template/Mask_2mm_dil2.mnc")
    outvol[] <- 0
    outvol[maskvol > 0.5] <- as.numeric(prev_final3[t,])
    mincWriteVolume(outvol, paste0("./tissue_prevalence/",t,"_",tissue[t],"_prevalence_ses3.mnc"), clobber=TRUE, like="../../../UKB/temporary_template/Mask_2mm_dil2.mnc")
}

for(t in 1:length(tissue)) {
    outvol <- mincGetVolume("../UKB/temporary_template/avg.020_2mm.mnc")
    maskvol <- mincGetVolume("../UKB/temporary_template/Mask_2mm_dil2.mnc")
    outvol[] <- 0
    outvol[maskvol > 0.5] <- as.numeric(prev_final23[t,])
    mincWriteVolume(outvol, paste0("./tissue_prevalence/",t,"_",tissue[t],"_prevalence_ses23.mnc"), clobber=TRUE, like="../../../UKB/temporary_template/Mask_2mm_dil2.mnc")
}

# Write masks of prevalence >0 and >=100

for(t in 1:length(tissue)) {
    print(tissue[t])
    mask_1 = ifelse(prev_final23[t,] > 0, 1, 0)
    mask_100 = ifelse(prev_final23[t,] >= 100, 1, 0)
    fwrite(as.data.frame(mask_1), paste0("./tissue_prevalence/",t,"_",tissue[t],"_prevalence_ses23_min1.tsv"))
    fwrite(as.data.frame(mask_100), paste0("./tissue_prevalence/",t,"_",tissue[t],"_prevalence_ses23_min100.tsv"))
    print(sum(mask_1 == 1))
    print(sum(mask_100 == 1))

    outvol <- mincGetVolume("../UKB/temporary_template/avg.020_2mm.mnc")
    maskvol <- mincGetVolume("../UKB/temporary_template/Mask_2mm_dil2.mnc")
    outvol[] <- 0
    outvol[maskvol > 0.5] <- as.numeric(mask_1)
    mincWriteVolume(outvol, paste0("./tissue_prevalence/",t,"_",tissue[t],"_prevalence_ses23_min1.mnc"), clobber=TRUE, like="../../../UKB/temporary_template/Mask_2mm_dil2.mnc")

    outvol <- mincGetVolume("../UKB/temporary_template/avg.020_2mm.mnc")
    maskvol <- mincGetVolume("../UKB/temporary_template/Mask_2mm_dil2.mnc")
    outvol[] <- 0
    outvol[maskvol > 0.5] <- as.numeric(mask_100)
    mincWriteVolume(outvol, paste0("./tissue_prevalence/",t,"_",tissue[t],"_prevalence_ses23_min100.mnc"), clobber=TRUE, like="../../../UKB/temporary_template/Mask_2mm_dil2.mnc")
}


