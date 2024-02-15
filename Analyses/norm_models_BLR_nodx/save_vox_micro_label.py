
import pandas as pd
import polars as pl

names = ["FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM", "jacobians_abs", "jacobians_rel"]
# names = ["FA"]
# names = ["jacobians_abs", "jacobians_rel"]

samples = pd.read_csv("./visualization/random_samples.txt", header=None)
samples = samples-1

label = pl.read_csv(f"../../bison_matrices/ses2_Label_whole_brain.tsv", has_header=False, columns=samples[0].tolist(), separator="\t").to_pandas()
label.to_csv(f"./visualization/vox_to_viz_label.tsv", index=False, na_rep="NA")

for n in names:
    print(n)
    micro = pl.read_csv(f"../../micro_matrices/ses2_{n}.tsv", has_header=False, columns=samples[0].tolist(), separator="\t").to_pandas()
    micro.to_csv(f"./visualization/{n}_vox_to_viz_micro.tsv", index=False, na_rep="NA")



