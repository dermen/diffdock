
import glob
import numpy as np
import pandas
import os

fnames = glob.glob("zinc_smi/*/*pkl")
Nsamp = 1e5
out_csv = "3i28_100k.csv"
pdb = os.path.abspath("3i28.pdb" )

Ns = {}
for i_f, f in enumerate(fnames):
    df = pandas.read_pickle(f)
    Ns[f] = len(df)
    print(i_f, f, len(df))
Ntot = sum(Ns.values())
Nsamp_per = {f: int(np.ceil(N/Ntot*Nsamp)) for f,N in Ns.items()}
all_samps = []
for f in Nsamp_per:
    df = pandas.read_pickle(f)
    perm = np.random.permutation(len(df))
    n = Nsamp_per[f]
    df_samp = df.iloc[perm[:n]]
    all_samps.append(df_samp)
    print(f)
df_all = pandas.concat(all_samps)
print("Total complexes: %d" % len(df_all))
df_all.reset_index(inplace=True, drop=True)
df_all['protein_path'] = pdb
df_all['ligand_description'] = df_all.smiles
df_all['complex_name'] = ["complex_%d"%i for  i in df_all.index]
df_all['protein_sequence'] = ""
df_all = df_all[['complex_name','protein_path','ligand_description','protein_sequence']]
df_all.to_csv(out_csv, index=False)
print(f"Wrote {out_csv}.")

