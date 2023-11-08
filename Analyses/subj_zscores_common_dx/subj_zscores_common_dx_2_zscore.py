
import datatable as dt
import pandas as pd
import numpy as np
import polars as pl
from joblib import Parallel, delayed
import multiprocessing as mp
import sys

# num_micro = int(sys.argv[1])
num_micro = 0
names = ["FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM"]
name = names[num_micro]

tissue_all=['Ventricules', 'CSF', 'Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH']
tissue_nm=['Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM']

# Load data

inclusions = pd.read_csv("../../../WMH_micro_spatial/QC/inclusions_only_dx_new.txt", header=None)
inclusions.columns = ['ID']

demo = dt.fread("../../../UKB/tabular/df_demo/UKBB_demo_wider.tsv").to_pandas()
demo = demo[(demo['InstanceID'] == 2)]
demo = demo[['SubjectID', 'Sex_31_0', 'Age_when_attended_assessment_centre_21003_0']]
demo.rename(columns={'SubjectID': 'ID', 'Sex_31_0':'Sex', 'Age_when_attended_assessment_centre_21003_0':'Age'}, inplace=True)
# demo['Sex'] = demo['Sex'].replace({'Female': 0, 'Male': 1}) # 0 = Female, 1 = Male
demo = pd.merge(inclusions, demo, on='ID', how='inner')

dx = pd.read_csv("../../../UKB/QC/exclusion_lists/dx_single.csv")

micro = pl.read_csv(f"./results/ses2_{name}_dx.tsv", has_header=False, separator="\t").to_pandas()
label = pl.read_csv(f"./results/ses2_Label_whole_brain_dx.tsv", has_header=False, separator="\t").to_pandas()

print(inclusions.shape)
print(demo.shape)
print(micro.shape)
print(label.shape)

nvox = micro.shape[1]

# Load NM results
nm = []
for t in range(len(tissue_nm)):
    print(tissue_nm[t])
    nm_tmp = pl.read_csv(f"../norm_models_BLR/results/nm_{name}_{tissue_nm[t]}.tsv", has_header=True, separator="\t").to_pandas()
    nm.append(nm_tmp)


def zscore_nm(i):
    print(i)
    demo_id = demo.iloc[i,:]
    micro_id = micro.iloc[i,:]
    label_id = label.iloc[i,:]
    zscores_id = np.full(nvox, np.nan)
    avg_col = f"{demo_id['Sex']}_mean_{int(demo_id['Age'])}"
    sd_col = f"{demo_id['Sex']}_sd_{int(demo_id['Age'])}"
    # Z-score one tissue type at a time
    for t in range(len(tissue_nm)):
        # print(tissue_nm[t])
        # print(t)
        tissue_idx = [i for i, x in enumerate(label_id) if x == (t+3)]
        avg = nm[t][avg_col][tissue_idx].reset_index(drop=True).astype('float64')
        sd = nm[t][sd_col][tissue_idx].reset_index(drop=True).astype('float64')
        micro_t = micro_id.iloc[tissue_idx].reset_index(drop=True).astype('float64')
        zscores_id[tissue_idx] = (micro_t - avg) / sd
    # Z-score WMHs with NAWM averages
    tissue_idx = [i for i, x in enumerate(label_id) if x == 9]
    avg = nm[5][avg_col][tissue_idx].reset_index(drop=True).astype('float64')
    sd = nm[5][sd_col][tissue_idx].reset_index(drop=True).astype('float64')
    micro_t = micro_id.iloc[tissue_idx].reset_index(drop=True).astype('float64')
    zscores_id[tissue_idx] = (micro_t - avg) / sd
    return(zscores_id)

# Run nm for all voxels (in chunk) and all brain tissue types using parallelized processes
num_cpus = mp.cpu_count()

results = Parallel(n_jobs=int(num_cpus/2))(delayed(zscore_nm)(i) for i in range(micro.shape[0]))
results = pd.DataFrame(results)
results.to_csv(f"./results/zscores_{name}.tsv", sep='\t', index=False, na_rep="NA")



