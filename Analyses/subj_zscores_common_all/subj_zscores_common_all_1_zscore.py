
import datatable as dt
import pandas as pd
import numpy as np
import polars as pl
from joblib import Parallel, delayed
import multiprocessing as mp
import sys
import subprocess
import os

# num_micro = 0
# num_chunk = 32

num_micro = int(sys.argv[1])
num_chunk = int(sys.argv[2])

names = ["FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM", "jacobians_abs", "jacobians_rel"]
name = names[num_micro]

tissue_all=['Ventricules', 'CSF', 'Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH']
tissue_nm=['Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM']

# Load data

inclusions = pd.read_csv("../../../WMH_micro_spatial/QC/inclusions_without_excluding_dx_new.txt", header=None)
inclusions.columns = ['ID']

demo = dt.fread("../../../UKB/tabular/df_demo/UKBB_demo_wider.tsv").to_pandas()
demo = demo[(demo['InstanceID'] == 2)]
demo = demo[['SubjectID', 'Sex_31_0', 'Age_when_attended_assessment_centre_21003_0']]
demo.rename(columns={'SubjectID': 'ID', 'Sex_31_0':'Sex', 'Age_when_attended_assessment_centre_21003_0':'Age'}, inplace=True)
# demo['Sex'] = demo['Sex'].replace({'Female': 0, 'Male': 1}) # 0 = Female, 1 = Male
demo = pd.merge(inclusions, demo, on='ID', how='inner')

# Get chunk of IDs
all_idx = [i for i in range(inclusions.shape[0])]
chunks = [all_idx[i:i + 1000] for i in range(0, len(all_idx), 1000)]

idx_chunk = chunks[num_chunk]

# Load NM results
nm = []
for t in range(len(tissue_nm)):
    print(tissue_nm[t])
    nm_tmp = pl.read_csv(f"../norm_models_BLR/results/nm_{name}_{tissue_nm[t]}.tsv", has_header=True, separator="\t").to_pandas()
    nm.append(nm_tmp)


def zscore_nm(i):
    print(i)
    demo_id = demo.iloc[i,:]
    micro_id = micro.iloc[i,:].astype('float64')
    label_id = label.iloc[i,:]
    zscores_id = np.full(nvox, np.nan)
    avg_colname = f"{demo_id['Sex']}_mean_{int(demo_id['Age'])}"
    sd_colname = f"{demo_id['Sex']}_sd_{int(demo_id['Age'])}"
    # Compute average and standard deviation columns for all tissue types
    avg_columns = [nm[t][avg_colname].astype('float64') for t in range(len(tissue_nm))]
    sd_columns = [nm[t][sd_colname].astype('float64') for t in range(len(tissue_nm))]
    # Z-score one tissue type at a time
    for t in range(len(tissue_nm)):
        print(tissue_nm[t])
        # print(t)
        tissue_idx = np.where(label_id == (t + 3))[0]
        avg = avg_columns[t][tissue_idx].reset_index(drop=True)
        sd = sd_columns[t][tissue_idx].reset_index(drop=True)
        micro_t = micro_id.iloc[tissue_idx].reset_index(drop=True)
        zscores_id[tissue_idx] = (micro_t - avg) / sd
    # Z-score WMHs with NAWM averages
    tissue_idx = np.where(label_id == 9)[0]
    avg = avg_columns[5][tissue_idx].reset_index(drop=True)
    sd = sd_columns[5][tissue_idx].reset_index(drop=True)
    micro_t = micro_id.iloc[tissue_idx].reset_index(drop=True)
    zscores_id[tissue_idx] = (micro_t - avg) / sd
    return(zscores_id)

#################### For raw maps ###############################

# Make chunk of micro and label with bash
ids_label = pd.read_csv("../../bison_matrices/ids_ses2_Label_whole_brain.txt", header=None)
# ids_micro = pd.read_csv(f"../../micro_matrices/ids_ses2_{name}.txt", header=None)

indices_label = ids_label[ids_label.iloc[:,0].isin(inclusions['ID'])].index[idx_chunk]
# indices_micro = ids_micro[ids_micro.iloc[:,0].isin(inclusions['ID'])].index[idx_chunk]

indices_label_plus1 = [i + 1 for i in indices_label]
# indices_micro_plus1 = [i + 1 for i in indices_micro]

np.savetxt(f"./tmp/idx_label_c{num_chunk}_{name}.txt", indices_label_plus1, fmt='%d')
# np.savetxt(f"./tmp/idx_micro_c{num_chunk}_{name}.txt", indices_micro_plus1, fmt='%d')

# Test with ID file (make sure same IDs as in chunk)
command = f"./select_rows.sh ../../bison_matrices/ids_ses2_Label_whole_brain.txt ./tmp/ids_label_c{num_chunk}_{name}.txt ./tmp/idx_label_c{num_chunk}_{name}.txt"
subprocess.call(command, shell=True)
# command = f"./select_rows.sh ../../micro_matrices/ids_ses2_{name}.txt ./tmp/ids_micro_c{num_chunk}_{name}.txt ./tmp/idx_micro_c{num_chunk}_{name}.txt"
# subprocess.call(command, shell=True)

# Make sub-matrices of subjects x voxels
command = f"./select_rows.sh ../../bison_matrices/ses2_Label_whole_brain.tsv ./tmp/label_c{num_chunk}_{name}.tsv ./tmp/idx_label_c{num_chunk}_{name}.txt"
subprocess.call(command, shell=True)
# command = f"./select_rows.sh ../../micro_matrices/ses2_{name}.tsv ./tmp/micro_c{num_chunk}_{name}.tsv ./tmp/idx_micro_c{num_chunk}_{name}.txt"
# subprocess.call(command, shell=True)

label = pl.read_csv(f"./tmp/label_c{num_chunk}_{name}.tsv", has_header=False, separator="\t").to_pandas()
# micro = pl.read_csv(f"./tmp/micro_c{num_chunk}_{name}.tsv", has_header=False, separator="\t").to_pandas()

print(inclusions.shape)
print(demo.shape)
# print(micro.shape)
print(label.shape)

nvox = label.shape[1]

# Run nm for all voxels (in chunk) and all brain tissue types using parallelized processes
# num_cpus = mp.cpu_count()

# results = Parallel(n_jobs=int(num_cpus/2))(delayed(zscore_nm)(i) for i in range(micro.shape[0]))
# results = pd.DataFrame(results)
# results.to_csv(f"./results/zscores_c{num_chunk}_{name}.tsv", sep='\t', index=False, na_rep="NA")

#################### For denoised maps ##########################

# Make chunk of micro and label with bash
ids_micro = pd.read_csv(f"../../micro_matrices/ids_ses2_{name}_anlm.txt", header=None)
indices_micro = ids_micro[ids_micro.iloc[:,0].isin(inclusions['ID'])].index[idx_chunk]
indices_micro_plus1 = [i + 1 for i in indices_micro]
np.savetxt(f"./tmp/idx_micro_c{num_chunk}_{name}_anlm.txt", indices_micro_plus1, fmt='%d')

# Test with ID file (make sure same IDs as in chunk)
command = f"./select_rows.sh ../../micro_matrices/ids_ses2_{name}_anlm.txt ./tmp/ids_micro_c{num_chunk}_{name}_anlm.txt ./tmp/idx_micro_c{num_chunk}_{name}_anlm.txt"
subprocess.call(command, shell=True)
# Make sub-matrices of subjects x voxels
command = f"./select_rows.sh ../../micro_matrices/ses2_{name}_anlm.tsv ./tmp/micro_c{num_chunk}_{name}_anlm.tsv ./tmp/idx_micro_c{num_chunk}_{name}_anlm.txt"
subprocess.call(command, shell=True)

# label = pl.read_csv(f"./tmp/label_c{num_chunk}_{name}.tsv", has_header=False, separator="\t").to_pandas()
micro = pl.read_csv(f"./tmp/micro_c{num_chunk}_{name}_anlm.tsv", has_header=False, separator="\t").to_pandas()

print(inclusions.shape)
print(demo.shape)
print(micro.shape)
print(label.shape)

# Calculate
num_cpus = mp.cpu_count()

results = Parallel(n_jobs=int(num_cpus-5))(delayed(zscore_nm)(i) for i in range(micro.shape[0]))
results = pd.DataFrame(results)
results.to_csv(f"./results/zscores_c{num_chunk}_{name}_anlm.tsv", sep='\t', index=False, na_rep="NA")

# Clean tmp directory
if os.path.exists(f"./tmp/label_c{num_chunk}_{name}.tsv") and name != "FA":
    os.remove(f"./tmp/label_c{num_chunk}_{name}.tsv")

if os.path.exists(f"./tmp/micro_c{num_chunk}_{name}.tsv"):
    os.remove(f"./tmp/micro_c{num_chunk}_{name}.tsv")

if os.path.exists(f"./tmp/micro_c{num_chunk}_{name}_anlm.tsv"):
    os.remove(f"./tmp/micro_c{num_chunk}_{name}_anlm.tsv")




