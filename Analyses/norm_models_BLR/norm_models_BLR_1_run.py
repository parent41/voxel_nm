
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# import seaborn as sns
# import joypy
# from sklearn.model_selection import train_test_split
from pcntoolkit.normative import estimate, evaluate
from pcntoolkit.util.utils import create_bspline_basis, compute_MSLL
import datatable as dt
import os
from itertools import product
import random
import math
from joblib import Parallel, delayed
import multiprocessing as mp
import polars as pl
import copy

random.seed(123)

# Paper: https://www.nature.com/articles/s41596-022-00696-5
# Documentation: https://pcntoolkit.readthedocs.io/en/latest/pages/BLR_normativemodel_protocol.html
# Github: https://github.com/predictive-clinical-neuroscience/PCNtoolkit-demo/tree/main
# Location of PCN toolbox: /gpfs/fs1/home/m/mchakrav/parent41/.virtualenvs/PCN_env/lib/python3.9/site-packages/pcntoolkit/

num_micro = int(sys.argv[1])
num_chunk = int(sys.argv[2])
# num_micro = 0
# num_chunk = 0

names = ["FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM"]
tissues = ['Ventricules', 'CSF', 'Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM', 'WMH']
name = names[num_micro]

# Parallelize chunks of voxels
mask = dt.fread(f"../../tissue_prevalence/1_{tissues[0]}_prevalence_ses23_min100.tsv").to_numpy()
all_idx = [i for i in range(mask.shape[1])]
chunks = [all_idx[i:i + 10000] for i in range(0, len(all_idx), 10000)]

print(f'\nMicro: {name}\nChunk: from {chunks[num_chunk][0]} to {chunks[num_chunk][-1]}')

# Load demographic data
ids = dt.fread("../../../WMH_micro_spatial/QC/inclusions_excluding_dx_new.txt").to_pandas()
ids.rename(columns={'C0': 'ID'}, inplace=True)

demo = dt.fread("../../../UKB/tabular/df_demo/UKBB_demo_wider.tsv").to_pandas()
demo = demo[(demo['InstanceID'] == 2)]
demo = demo[['SubjectID', 'Sex_31_0', 'Age_when_attended_assessment_centre_21003_0']]
demo.rename(columns={'SubjectID': 'ID', 'Sex_31_0':'Sex', 'Age_when_attended_assessment_centre_21003_0':'Age'}, inplace=True)
demo['Sex'] = demo['Sex'].replace({'Female': 0, 'Male': 1}) # 0 = Female, 1 = Male
demo = pd.merge(ids, demo, on='ID', how='inner')

# Read all voxel chunk of micro and label matrices

micro = pl.read_csv(f"../../micro_matrices/ses2_{name}.tsv", has_header=False, columns=chunks[num_chunk], separator="\t").to_pandas()
ids_micro = dt.fread(f"../../micro_matrices/ids_ses2_{name}.txt").to_pandas()

label = pl.read_csv(f"../../bison_matrices/ses2_Label_whole_brain.tsv", has_header=False, columns=chunks[num_chunk], separator="\t").to_pandas()
ids_label = dt.fread(f"../../bison_matrices/ids_ses2_Label_whole_brain.txt").to_pandas()

print(micro.shape)
print(label.shape)

# Indices of IDs to remove
common_ids = set(ids['ID']).intersection(set(demo['ID']), set(ids_micro['C0']), set(ids_label['C0']))

exclude_indices_inc = np.where(~ids['ID'].isin(common_ids))[0]
exclude_indices_demo = np.where(~demo['ID'].isin(common_ids))[0]
exclude_indices_micro = np.where(~ids_micro['C0'].isin(common_ids))[0]
exclude_indices_label = np.where(~ids_label['C0'].isin(common_ids))[0]

ids = ids.drop(exclude_indices_inc).reset_index(drop=True)
demo = demo.drop(exclude_indices_demo).reset_index(drop=True)
micro = micro.drop(exclude_indices_micro).reset_index(drop=True)
label = label.drop(exclude_indices_label).reset_index(drop=True)

# Add intercept and B-spline with age to demo
B = create_bspline_basis(demo['Age'].min(), demo['Age'].max(), p=4)
Phi = np.array([B(i) for i in demo.iloc[:,2]])

intercept_pd = pd.DataFrame(np.ones((demo.shape[0],1)), columns = ["Intercept"])
Phi_pd = pd.DataFrame(Phi, columns=[f'Bspline_{i+1}' for i in range(Phi.shape[1])])

demo_pd = pd.concat([demo.iloc[:, 1:], intercept_pd, Phi_pd], axis=1)

# All possibilities of age (integers) and sex, with intercept
possibilities_demo = list(product([0,1], list(range(demo['Age'].min().astype(int), demo['Age'].max().astype(int)+1))))
possibilities_demo = np.hstack((possibilities_demo, np.ones((len(possibilities_demo), 1), dtype=int)))
Phi_demo = np.array([B(i) for i in possibilities_demo[:,1]]) # Add 4th order B-spline to age
possibilities_demo = np.concatenate((possibilities_demo, Phi_demo), axis=1) # Add intercept column to iv

# Output file: fit metrics + age- and sex-specific means and SDs
metrics_columns = ['ROI', 'N', 'MSLL', 'EXPV', 'SMSE', 'RMSE', 'Rho', 'R2']
male_means_columns = [f'Male_mean_{i+1}' for i in range(int(demo['Age'].min())-1, int(demo['Age'].max()))]
male_sd_columns = [f'Male_sd_{i+1}' for i in range(int(demo['Age'].min())-1, int(demo['Age'].max()))]
female_means_columns = [f'Female_mean_{i+1}' for i in range(int(demo['Age'].min())-1, int(demo['Age'].max()))]
female_sd_columns = [f'Female_sd_{i+1}' for i in range(int(demo['Age'].min())-1, int(demo['Age'].max()))]
all_columns = metrics_columns + male_means_columns + male_sd_columns + female_means_columns + female_sd_columns

# results = pd.DataFrame(np.zeros((micro.shape[1], len(all_columns))), columns = all_columns)
# results = results.reset_index(drop=True)

# Iterate over voxels
# for i in range(micro.shape[1]):
def run_nm(i):
    print(i)
    # i = 23045
    vox = micro_tissue.iloc[:,i] # select vox
    count_tissue_voxels = vox.shape[0] - (vox.isna().sum().sum())
    print(f'Number of voxels for specific tissue type in voxel = ',count_tissue_voxels)
    res = pd.DataFrame(np.zeros((1, len(all_columns))), columns=all_columns)
    # If less than 100 voxels of the right tissue type, return NaNs
    if count_tissue_voxels < 100:
        res.iloc[0,:] = np.full(len(all_columns), np.nan)
        res['ROI'] = i
        res['N'] = count_tissue_voxels
    # If more or equal to 100 voxels of the right tissue type, return NaNs
    if count_tissue_voxels >= 100:
        # Concatenate with demo
        vox = pd.concat([demo_pd, vox], axis=1)
        # Drow rows where not tissue type (NaN)
        vox = vox.dropna(subset=[vox.columns[-1]]).reset_index(drop=True)
        # # train_test split
        # X_train, X_test, y_train, y_test = train_test_split(vox.iloc[:, :-1], vox.iloc[:, -1], test_size=0.2, random_state=123)
        # np.savetxt(os.path.join('tmp/cov_bspline_train.txt'), X_train.values)
        # np.savetxt(os.path.join('tmp/cov_bspline_test.txt'), X_test.values)
        # np.savetxt(os.path.join('tmp/resp_train.txt'), y_train.values)
        # np.savetxt(os.path.join('tmp/resp_test.txt'), y_test.values)
        # Save txt files of covariates and response variables
        np.savetxt(os.path.join(f'tmp/{name}_{tissues[t]}_c{num_chunk}_cov_bspline_{i}.txt'), vox.iloc[:,:-1].values)
        np.savetxt(os.path.join(f'tmp/{name}_{tissues[t]}_c{num_chunk}_resp_{i}.txt'), vox.iloc[:,-1].values)
        # Estimate normative model
        # cov_file_train = os.path.join('tmp/cov_bspline_train.txt')
        # resp_file_train = os.path.join('tmp/resp_train.txt')
        # cov_file_test = os.path.join('tmp/cov_bspline_test.txt')
        # resp_file_test = os.path.join('tmp/resp_test.txt')
        cov_file = os.path.join(f'tmp/{name}_{tissues[t]}_c{num_chunk}_cov_bspline_{i}.txt')
        resp_file = os.path.join(f'tmp/{name}_{tissues[t]}_c{num_chunk}_resp_{i}.txt')
        yhat, s2, nm, Z, metrics_te = estimate(cov_file, resp_file, testresp=resp_file, testcov=cov_file, alg = 'blr', optimizer = 'powell', cvfolds = None, savemodel = False, saveoutput = False, standardize = False)
        os.remove(f'tmp/{name}_{tissues[t]}_c{num_chunk}_cov_bspline_{i}.txt')
        os.remove(f'tmp/{name}_{tissues[t]}_c{num_chunk}_resp_{i}.txt')
        # metrics = [i, metrics_te['MSLL'][0], metrics_te['EXPV'][0], metrics_te['SMSE'][0], metrics_te['RMSE'][0], metrics_te['Rho'][0]]
        # vox['yhat'] = yhat
        # subset_test = vox[(vox.iloc[:, 0] == 0) & (vox.iloc[:, 1] == 64)]
        # Test predict for all possibilities of demographics
        pred = np.transpose(nm.predict(possibilities_demo))
        demo_pred = np.concatenate((possibilities_demo[:,[0,1]], pred), axis=1) # Add intercept column to iv
        demo_pred = pd.DataFrame(data = demo_pred, columns = ['Sex','Age','Mean', 'Var'])
        demo_pred['SD'] = np.sqrt(demo_pred['Var'])
        # Save results
        for idx,c in enumerate(res.columns, start=0):
            # print(c)
            if c == 'ROI':
                res[c] = i
            if c == 'N':
                res[c] = count_tissue_voxels
            if c == 'MSLL':
                res[c] = metrics_te['MSLL'][0]
            if c == 'EXPV':
                res[c] = metrics_te['EXPV'][0]
            if c == 'SMSE':
                res[c] = metrics_te['SMSE'][0]
            if c == 'RMSE':
                res[c] = metrics_te['RMSE'][0]
            if c == 'Rho':
                res[c] = metrics_te['Rho'][0]
            if c == 'R2':
                res[c] = math.pow(metrics_te['Rho'][0],2)            
            if c.startswith("Male_"):
                age = int(res.columns[idx][-2:])
                sex = 1
                if 'mean' in c:
                    res[c] = demo_pred[(demo_pred['Sex'] == sex) & (demo_pred['Age'] == age)]['Mean'].values
                if 'sd' in c:
                    res[c] = demo_pred[(demo_pred['Sex'] == sex) & (demo_pred['Age'] == age)]['SD'].values
            elif c.startswith("Female"):
                age = int(res.columns[idx][-2:])
                sex = 0
                if 'mean' in c:
                    res[c] = demo_pred[(demo_pred['Sex'] == sex) & (demo_pred['Age'] == age)]['Mean'].values
                if 'sd' in c:
                    res[c] = demo_pred[(demo_pred['Sex'] == sex) & (demo_pred['Age'] == age)]['SD'].values
    return(res)

# Run nm for all voxels (in chunk) and all brain tissue types using parallelized processes
num_cpus = mp.cpu_count()

for t in range(2,8):
    print(f"\n--------------------------\nTissue type = {tissues[t]}\n--------------------------\n")
    print(t)
    print((label==(t+1)).sum().sum())
    micro_tissue = copy.copy(micro)
    micro_tissue[label != (t+1)] = np.nan
    print((micro_tissue > 0).sum().sum())
    results = Parallel(n_jobs=int(num_cpus/2))(delayed(run_nm)(i) for i in range(micro_tissue.shape[1]))
    results = pd.concat(results, ignore_index=True)
    results.to_csv(f"./results/nm_{name}_{tissues[t]}_vox_{chunks[num_chunk][0]}_{chunks[num_chunk][-1]}.tsv", sep='\t', index=False, na_rep="NA")









# # Plot

# import plotnine as pn

# # fig = (pn.ggplot(vox, pn.aes(x="Age", y="C23045", color="factor(Sex)")) + 
# fig = (pn.ggplot() +
# pn.geom_point(test, pn.aes(x="Age", y="column_184073"), color="blue", alpha=0.01) + 
# pn.geom_line(demo_pred, pn.aes(x="Age", y="Mean", color="factor(Sex)")) +
# pn.geom_line(test_pred, pn.aes(x="Age", y=(test_pred['Mean'] + (test_pred['SD'])), color="factor(Sex)"), alpha=0.8) +
# pn.geom_line(test_pred, pn.aes(x="Age", y=(test_pred['Mean'] - (test_pred['SD'])), color="factor(Sex)"), alpha=0.8) +
# pn.geom_line(test_pred, pn.aes(x="Age", y=(test_pred['Mean'] + (test_pred['SD']*2)), color="factor(Sex)"), alpha=0.5) +
# pn.geom_line(test_pred, pn.aes(x="Age", y=(test_pred['Mean'] - (test_pred['SD']*2)), color="factor(Sex)"), alpha=0.5) +
# pn.geom_line(test_pred, pn.aes(x="Age", y=(test_pred['Mean'] + (test_pred['SD']*3)), color="factor(Sex)"), alpha=0.2) +
# pn.geom_line(test_pred, pn.aes(x="Age", y=(test_pred['Mean'] - (test_pred['SD']*3)), color="factor(Sex)"), alpha=0.2)
# )
# fig.save("test.png", dpi=300)

# Trying different ways to load huge csv, polars was most fastest
# import pyarrow.csv as pc
# import pyarrow as pa
# import pyarrow.parquet as pq
# from dask import dataframe as dd
# from dask.diagnostics import ProgressBar

# micro = dt.fread(f"../../micro_matrices/ses2_{name}.tsv", sep="\t", columns=mask_bool.tolist())
# micro = pd.read_csv(f"../../micro_matrices/ses2_{name}.tsv", sep="\t", header=None, usecols=mask_idx)
# micro = micro.read()
# micro = pq.read_table(f"../../micro_matrices/ses2_{name}.tsv", use_threads=True, columns=mask[0].tolist())
# micro = pc.read_csv(f"../../micro_matrices/ses2_{name}.tsv", parse_options=pc.ParseOptions(delimiter="\t"), convert_options=pc.ConvertOptions(include_columns=mask_idx_40))
# micro = dd.read_csv(f"../../micro_matrices/ses2_{name}.tsv", sep="\t", sample = 25600000, index_col=False, usecols=list(range(40)))
# micro_pd = micro.compute()
# micro = dd.read_csv(f"../../micro_matrices/ses2_{name}.tsv", sep="\t", sample = 25600000)


