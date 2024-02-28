
library(data.table)
library(RMINC)
# install.packages("Hmisc")
library(Hmisc)
library(foreach)
library(doParallel)

args = commandArgs(trailingOnly=TRUE)

# dir.create("./visualization/indiv", showWarnings=FALSE)
dir.create("./results/indiv", showWarnings=FALSE)

# args = c()
# args[1] = 32

num_chunk = as.numeric(args[1])

names = c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM", "jacobians_abs", "jacobians_rel")
names_long = c("dti_MD_UKB", "NODDI_ISOVF_UKB", "dti_FA_UKB", "NODDI_ICVF_UKB", "NODDI_OD_UKB", "T2star_UKB", "QSM_UKB", "jacobians_abs_2mm", "jacobians_rel_2mm")
tissue_all=c('Ventricules', 'CSF', 'Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH')

mask = mincGetVolume("../../../UKB/temporary_template/Mask_2mm_dil2.mnc")

# Load data

inclusions = as.data.frame(fread(paste0("../subj_zscores_common_all/tmp/ids_micro_c",num_chunk,"_FA_anlm.txt")))
colnames(inclusions) = "ID"
ids = inclusions$ID

# Spatial correlations

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    Marker_1 = rownames(cormat)[row(cormat)[ut]],
    Marker_2 = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}

get_spatial_corrs = function(label_df, micro_df) {
    
    # First batch of results with label = 3
    micro_tmp = micro_df[which(label_df[mask[]>0.5] == 3),]
    corrs = rcorr(as.matrix(micro_tmp))
    corrs = flattenCorrMatrix(corrs$r, corrs$P)
    corrs$Label = tissue_all[3]
    result = corrs

    # For the rest of the labels...
    for (l in 4:9) {
        micro_tmp = micro_df[which(label_df[mask[]>0.5] == l),]
        corrs = rcorr(as.matrix(micro_tmp))
        corrs = flattenCorrMatrix(corrs$r, corrs$P)
        corrs$Label = tissue_all[l]
        result = rbind(result, corrs)
    }

    # Across the whole brain
    micro_tmp = micro_df[which(label_df[] %in% seq(3,9)),]
    corrs = rcorr(as.matrix(micro_tmp))
    corrs = flattenCorrMatrix(corrs$r, corrs$P)
    corrs$Label = "Whole_brain"
    result = rbind(result, corrs)

    # Across the whole cerebral white matter
    micro_tmp = micro_df[which(label_df[] %in% seq(8,9)),]
    corrs = rcorr(as.matrix(micro_tmp))
    corrs = flattenCorrMatrix(corrs$r, corrs$P)
    corrs$Label = "Cerebral_WM"
    result = rbind(result, corrs)

    return(result)
}

# Run
cores = detectCores()
parallelCluster <- makeCluster(cores-5, type = "SOCK", methods = FALSE)
setDefaultCluster(parallelCluster)
registerDoParallel(parallelCluster)

foreach(i=1:length(ids), .packages = c('RMINC', 'data.table', 'Hmisc')) %dopar% {
# for (i in 1:length(ids)) {
    cat(paste0("\nID = ", ids[i], "\n"))

    label_id = round(mincGetVolume(paste0("../../../WMH_micro_spatial/maps_UKB_space/sub-",ids[i],"_ses-2_Label_UKB.mnc")),0)

    # For zscores
    zscores_id = as.data.frame(matrix(nrow=length(label_id[mask>0.5]), ncol=length(names)))
    colnames(zscores_id) = names
    for (n in 1:length(names)) {
        zscores_id_tmp = mincGetVolume(paste0("../subj_zscores_common_all/results/mnc/",ids[i],"_zscore_anlm_",names[n],".mnc"))
        zscores_id[,n] = zscores_id_tmp[mask[]>0.5]
    }

    corrs_id = get_spatial_corrs(label_id, zscores_id)
    corrs_id = cbind(ids[i], corrs_id)
    colnames(corrs_id)[1] = "ID"

    fwrite(as.data.frame(corrs_id), paste0("./results/indiv/",ids[i],"_spatial_corrs_zscores.tsv"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

    # For raw micro
    micro_id = as.data.frame(matrix(nrow=length(label_id[mask>0.5]), ncol=length(names)))
    colnames(micro_id) = names
    for (n in 1:length(names)) {
        micro_id_tmp = mincGetVolume(paste0("../../maps_UKB_space_anlm_all/sub-",ids[i],"_ses-2_",names_long[n],"_anlm.mnc"))
        micro_id[,n] = micro_id_tmp[mask[]>0.5]
    }

    corrs_id = get_spatial_corrs(label_id, micro_id)
    corrs_id = cbind(ids[i], corrs_id)
    colnames(corrs_id)[1] = "ID"

    fwrite(as.data.frame(corrs_id), paste0("./results/indiv/",ids[i],"_spatial_corrs_raw.tsv"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

}



