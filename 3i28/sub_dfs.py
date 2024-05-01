import os
import pandas

# on 4 NERSC nodes, each with 4 GPUs, 4000 complexes can be run in ~2-3 hours

master_csv = "3i28_100k.csv"  # 100,000 complexes 
df = pandas.read_csv(master_csv)
n_csv = 25  # so ~ 4000 complexes per CSV
sub_rootdir="dfs"  # this folder will be created to store the output and intermediate files
submit_jobs = True  # set to False to simply print the sbatch commands

if not os.path.exists(sub_rootdir):
	os.makedirs(sub_rootdir)

dfs = np.array_split(df,n_csv)
for i,d in enumerate(dfs):
    print(d.iloc[0])
    name = "%s/d%d.csv" % (sub_rootdir,i)
    d.to_csv(name, index=False)

for i in range(n_csv):
    cmd = 'sbatch -C gpu -J d{d} -N4 -A m4326_g --cpus-per-gpu=8 --ntasks-per-node=32 --gpus-per-node=4 -q regular -t 240 -o {dirname}/d{d}.out -e {dirname}/d{d}.err --wrap="time srun -N4 -n16 --ntasks-per-node=4 --cpus-per-gpu=1 --gpus-per-node=4 -c8 python inference.py --config default_inference_args.yaml --protein_ligand_csv {dirname}/d{d}.csv --out_dir {dirname}/d{d} --ndev 4"'
    cmd = cmd.format(d=i, dirname=sub_rootdir)
    print(cmd)
    print("")
    if submit_jobs:
		os.system(cmd)
