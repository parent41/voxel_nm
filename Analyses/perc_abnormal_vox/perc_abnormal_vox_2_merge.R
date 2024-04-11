
library(RMINC)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(data.table)
library(ggpubr)

tissue_nm = c('Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM')
tissue_all=c('Ventricules', 'CSF', 'Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH')

names = c("MD", "ISOVF", "FA", "ICVF", "OD", "T2star", "QSM", "jacobians_rel", "jacobians_abs")

# Load perc abnormal vox results

results_fl = list.files("./results/raw", pattern="*", full.names = TRUE)

results_list = list()
for (f in 1:length(results_fl)) {
    print(results_fl[f])

    if (grepl("jacobians", results_fl[f])) {
      micro_res = sub(".*_(?:.*_){2}(.*_.*?)_.*", "\\1", results_fl[f])
    } else {
      micro_res = sub(".*_(.*)_anlm.tsv", "\\1", results_fl[f])
    }

    results_list[[f]] = as.data.frame(fread(results_fl[f]))
    results_list[[f]]$micro = micro_res
    # results = rbind(results, results_tmp)
}

results <- do.call(rbind, results_list)

results = results[order(results$ID),]

# Add WMH+NAWM total abnormality

results <- as.data.frame(results %>%
  group_by(ID) %>%
  group_modify(
    ~{
      # Extract the subset for label_value 8 and 9 within each ID
      subset_8_9 <- .x %>%
        filter(label_value %in% c(8, 9)) %>%
        select(-label_value) %>%
        mutate(label_value = 10) %>%
        group_by(threshold, micro) %>%
        summarize(
          # ID = unique(ID),
          count_vox_label = sum(count_vox_label),
          count_vox_above_thresh = sum(count_vox_above_thresh),
          perc_vox_above_thresh = NA
        )
      
      # Combine the original dataset with the subset_8_9 within each ID
      bind_rows(
        .x,
        subset_8_9
      ) %>%
        arrange(threshold, micro, label_value)
    }
  ))

results$label_value[is.na(results$label_value)] <- 10
results$perc_vox_above_thresh[results$label_value == 10] = results$count_vox_above_thresh[results$label_value == 10] / results$count_vox_label[results$label_value == 10]

fwrite(results, "./results/results.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
