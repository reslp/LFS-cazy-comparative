#!/usr/bin/env python

import argparse
import glob
import sys
import os
from Bio import SeqIO

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="select_complete_buscos.py", description = """Converts the orthology file from funannotate to a CAFE input file""", epilog = """written by Philipp Resl 2019""")
pars.add_argument('-orth', dest="orth", required=True, help="Path to funannotate orthology_groups.txt file")
pars.add_argument('-ids', dest="ids", required=True, help="File with id mappings")
pars.add_argument('-o', dest="outfile", required=False, help="Output file path and name")
args=pars.parse_args()


print("Opening IDs file", args.ids)
idfile = open(args.ids, "r")

ids = {}
for line in idfile:
	names = line.strip().split("\t")
	ids[names[0]] = names[1]
	
print("Found", len(ids), "ids in file.")
print(sorted(ids))
print("Parsing othology file", args.orth)
count = 0
gene_number = 1
header = "Desc\tFamily ID"

for key in sorted(ids):
	header += "\t"
	header += ids[key].replace("_", "")

outfile = open(args.outfile, "w")
outfile.write(header +"\n")
	
for line in open(args.orth, "r"):
	items = line.strip().split("\t")
	outstring = items[0] + "\t" + str(gene_number)
	seqnames = items[-1].split(", ")
	print(len(seqnames))
	for key in sorted(ids):
		for seq in seqnames:
			if seq.startswith(key):
				count+=1
		outstring += "\t"
		outstring += str(count)
		count = 0
	#outstring = outstring.strip("\t")
	outstring += "\n"
	outfile.write(outstring)
	outstring = ""
	gene_number += 1
		