
"""
Converts a bunch of smi files into pandas dataframes with smiles strings 
"""

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("smi_glob", type=str, help="glob for .smi files")
args = parser.parse_args()

import os
from mpi4py import MPI
COMM = MPI.COMM_WORLD
from rdkit import Chem
import glob
import pandas

smis = glob.glob(args.smi_glob)
assert smis
Nsmi = len(smis)
for i_smi, smi in enumerate(smis):
    outfile = smi.replace(".smi", ".pkl") 
    if os.path.exists(outfile):
        print(f"File {outfile} exists! Continuing...")
        continue
    supp = Chem.rdmolfiles.SmilesMolSupplier(smi)
    records = []
    Nrec = len(supp)
    seen = 0
    for i in range(Nrec):
        if i % COMM.size != COMM.rank:
            continue
        a,b = supp.GetItemText(i).strip().split()
        records.append((a,b))
        seen += 1
        if COMM.rank==0 and seen%100==0:
            print( f"{smi} ({i_smi+1}/{Nsmi}) record {i+1}/{Nrec}:",a, end="\r", flush=True)
    records = COMM.reduce(records)
    if COMM.rank==0:
        a,b = zip(*records)
        df = pandas.DataFrame({"smiles":a,"zinc":b})
        df.to_pickle(outfile)
        print(f"Wrote {outfile}.")

#pdb = os.path.abspath("3i28.pdb" )
#out_csv = "3i28_100k.csv"
#df['protein_path'] = pdb
#df['ligand_description'] = df.smiles
#df['complex_name'] = ["complex_%d"%i for  i in df.index]
#df['protein_sequence'] = ""
#df2 = df[['complex_name','protein_path','ligand_description','protein_sequence']]
#df2.to_csv(out_csv, index=False)

