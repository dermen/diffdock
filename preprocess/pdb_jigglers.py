import argparse
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("pdbin", type=str, help="Name of the input PDB file")
parser.add_argument("outdir", type=str, help="output folder (created if missing)")
parser.add_argument("--njiggles", type=int, default=50, help="number of jiggle pdbs written to output folder")
parser.add_argument("--seed", type=int,default=None, help="default random seed to use for every jiggly")
parser.add_argument("--shiftRange", type=float, nargs=2, default=[0,5], help="An interval where , for each jiggle, a uniform random number will be drawn, and where that number will then specify the strength of the jiggle (0 means no jiggle)")
parser.add_argument("--pathToJiggleScript", help="path to the jigglepbd awk script", default="./jigglepdb.awk", type=str)
args = parser.parse_args()
import os
import random
os.makedirs(args.outdir, exist_ok=True)

jigglepdb=os.path.abspath(args.pathToJiggleScript)
assert os.path.exists(jigglepdb)
shift_min, shift_max = args.shiftRange
assert shift_max > shift_min, "shiftRange should be --shiftRange minval maxval"
print(f"PDB Jigglers will do {args.njiggles} jiggles")
print("Using shift range:", args.shiftRange)
pdb_code = os.path.basename(args.pdbin).split(".")[0]

cmd=f"""
{jigglepdb} \
 -v shift=byB -v seed=%d \
 -v shift=%f \
 -v frac_thrubond=0.9 -v ncyc_thrubond=500 \
 -v frac_magnforce=1.1 -v ncyc_magnforce=500 \
 -v shift_scale=1 \
 -v dry_shift_scale=1 {args.pdbin} > %s
"""

for i_jiggle in range(args.njiggles):
    jiggled_pdb=os.path.join(args.outdir, "%s_jiggled_%d.pdb" % (pdb_code, i_jiggle))
    seed = random.getrandbits(16) if args.seed is None else args.seed
    shift = random.uniform(*args.shiftRange)
    
    logfile = jiggled_pdb.replace(".pdb", ".sh")
    cmd_i =cmd % (seed, shift, jiggled_pdb) 
    print(cmd_i)

    os.system(f"echo '{cmd_i}' > {logfile};{cmd_i}")
    print(f"jiggles done so far: {i_jiggle+1}/{args.njiggles}")

print(f"Copying reference PDB {args.pdbin} to {args.outdir}")
os.system(f"cp {args.pdbin} {args.outdir}")

print("Done.")
