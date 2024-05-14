import pandas as pd
import sys
import os

args = sys.argv
scoref = args[1]
bimf = args[2]
outf=args[3]

print("reading..")
score = pd.read_csv(scoref, sep="\t", dtype=str, header=None)
score.columns = ["CHR", "SNP", "BP", "A1", "A2","score"]

bim_df = pd.read_csv(bimf, sep="\t", header=None)
bim_df.columns = ["chr", "rsID","GP",  "position", "allele1", "allele2"]


# Create a dictionary to map rsID to position from bim_df
rsid_to_position = dict(zip(bim_df["rsID"], bim_df["position"]))
# Map rsID from score_df to position from bim_df
score["BP"] = score["SNP"].map(rsid_to_position).astype(str)

print("converting..")
score["SNP"] = score["CHR"] + ":" + score["BP"] + ":" + score[["A1", "A2"]].apply(lambda x: ":".join(sorted(x)), axis=1)

print("writing..")
score.to_csv(outf, sep="\t", index=None, header=None, encoding="ascii")

