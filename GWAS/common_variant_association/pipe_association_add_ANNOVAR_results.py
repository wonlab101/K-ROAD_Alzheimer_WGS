### add rsID from annovar output

import pandas as pd
import numpy as np
from datetime import datetime
import sys


args = sys.argv
annovar1=args[1]
annovar2=args[2]
f=args[3]
outf=args[4]

markercol = "ID"     # column name of marker ID
A1col = "ALLELE0"    # A1/A2 col in assoc file which used to Alt/Ref on making ANNOVAR input
A2col = "ALLELE1"

print("%s    read data" % (datetime.now()))
assoc = pd.read_csv(f, sep="\t", dtype=str)

av1 = pd.read_csv(annovar1, sep="\t", dtype=str)
av2 = pd.read_csv(annovar2, sep="\t", dtype=str)

# setting
print("%s    data setting" % (datetime.now()))
av1 = av1.loc[av1["avsnp150"]!=".", :]
av2 = av2.loc[av2["avsnp150"]!=".", :]
n_av1 = len(av1)
print("av1: %s rows" % n_av1)

av1[A1col] = av1["Alt"]
av1[A2col] = av1["Ref"]
av1[markercol] = av1["Chr"].map(str) + ":" + av1["End"].map(str) + ":" +av1[A1col].map(str) + ":" + av1[A2col].map(str)

if len(av2)!=0:
    av2[A1col] = av2["Ref"]
    av2[A2col] = av2["Alt"]
    av2[markercol] = av2["Chr"].map(str) + ":" + av2["End"].map(str) + ":" +av2[A1col].map(str) + ":" + av2[A2col].map(str)
    av2 = av2.loc[~av2[markercol].isin(av1[markercol]), :]
    n_av2 = len(av2)
    print("av2: %s rows" % len(av2))
    if len(av2)!=0:
        av1["flipped"] = "."
        av2["flipped"] = "1"
        av = pd.concat([av1, av2])
    else:
        av = av1

else:
    av = av1
    n_av2 = 0

av = av.drop(av.columns.tolist()[:5], axis=1)


# get rsID
print("%s    get rsID" % (datetime.now()))
assoc = pd.merge(assoc, av, on=[markercol, A1col, A2col], how="left")


# write
print("%s    write data" % (datetime.now()))
assoc.to_csv(outf, index=None, header=True, sep="\t", encoding="ascii")
print("%s    Annotated:  %s rows  (%s missing info)" % (datetime.now(), len(av), len(assoc)-len(av)))

