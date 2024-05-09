# coding: utf-8
import glob
import os
from itertools import groupby
import pandas
import numpy as np
import sys
outdir = sys.argv[1]

print("Gathering all fnames, will take time...")
fnames = glob.glob("%s/d[0-9]*/proc*/complex*/*confidence*" % outdir)
#fnames += glob.glob("2k/proc*/com*/*confidence*")
print("Found %d fnames" % len(fnames))
# the inner most directory is the complex name... 
# info is a dict with key=complex name, value= list of rank*sdf files
key = lambda x: os.path.basename(os.path.dirname(x))
gb = groupby(sorted(fnames, key=key), key=key)
info = {k:list(v) for k,v in gb}
print("loading dataframe")
# combine the input tables:
df = pandas.concat([pandas.read_csv(f) for f in glob.glob("%s/*.csv"%outdir)])
    
print("Computing scores")
#Scores = {}
top_scores = []
sum_scores3 = []
sum_scores5 = []
sum_scores10 = []
names = []
rank1_files = []
rank2_files = []
rank3_files = []
for name in info:
    ranks = info[name]
    # sort files by increasing rank
    ranks = sorted(ranks, key=lambda x: int(os.path.basename(x).split("_")[0].split("rank")[1]))
    # extract the scores
    scores = [float(x.split("confidence")[1].split(".sdf")[0]) for x in ranks]
    #Scores[name] = scores
    top_scores.append(scores[0])
    sum_scores10.append(sum(scores))
    sum_scores5.append(sum(scores[:5]))
    sum_scores3.append(sum(scores[:3]))
    names.append(name)
    rank1_files.append(ranks[0])
    rank2_files.append(ranks[1])
    rank3_files.append(ranks[2])

df_score = pandas.DataFrame({"complex_name":names, "rank1_score":top_scores, 
    "top10_sum_scores": sum_scores10, "top3_sum_scores":sum_scores3, "top5_sum_scores": sum_scores5,
    "rank1_sdf": rank1_files, 'rank2_sdf':rank2_files, 'rank3_sdf':rank3_files})
df_m = pandas.merge(df, df_score, on="complex_name", how="inner")
df_m.to_pickle("%s/scores.pkl"%outdir)
df_m.to_csv("%s/scores.tsv"%outdir, sep="\t")

good = df_m.loc[df_m.top3_sum_scores > 0].reset_index(drop=True)
good['sdfs'] = [ open(f, 'r').read() for f in good.rank1_sdf]
# TODO: average the positive confidence poses ? 
# good['ave_sdfs'] = .. .. 
good.to_pickle("%s/confident_results.pkl" % outdir)
print(f"saved scores.tsv / scores.pkl and confident_results.pkl in {outdir}. ")

