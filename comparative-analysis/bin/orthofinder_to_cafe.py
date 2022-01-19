#!/usr/bin/env python

import argparse
import glob
import sys
import os
import pandas as pd
#from Bio import SeqIO

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="orthofinder_to_cafe.py", description = """Converts the orthology file from Orthofinder to a CAFE input file""", epilog = """written by Philipp Resl 2019""")
pars.add_argument('-i', dest="orth", required=True, help="Path to Orthofinder OrthoGroups.GeneCount.tsv file")
pars.add_argument('-o', dest="outfile", required=True, help="Output file path and name")
args=pars.parse_args()


print("Parsing Orthology file...")
data = pd.read_csv(args.orth, sep="\t")
colnames = data.columns.tolist()
# this line for old cafe versions which did not allow underscores:
#colnames = [name.replace(".proteins", "").replace("_", "") if "Orthogroup" not in name else "Family ID" for name in colnames]
colnames = [name.replace(".proteins", "") if "Orthogroup" not in name else "Family ID" for name in colnames]
data.columns = colnames
data.insert(loc=0, column="Desc", value="(null)")
colnames = data.columns.tolist()
data = data.drop(columns=["Total"])
print("Writing output to", args.outfile)
data.to_csv(args.outfile, sep="\t", index=False)
print("done")
