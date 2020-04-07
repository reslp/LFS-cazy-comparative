#!/usr/bin/env python

import argparse
import glob
import sys
import os
from Bio import SeqIO

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="select_complete_buscos.py", description = """This script alignments for orthologous groups clustered with MCL""", epilog = """written by Philipp Resl""")
pars.add_argument('-dump', dest="dump_path", required=True, help="Path to the dump file created with MCL")
pars.add_argument('-fasta', dest="fasta_path", required=True, help="Path to FASTA files containing sequences (used for all vs. all BLAST)")
pars.add_argument('-ids', dest="ids", required=True, help="Space separated ids used to identify species: eg: -ids \"xylpar xylsor homsap aratha\"")
pars.add_argument('-o', dest="outdir", required=True, help="Output directory")
args=pars.parse_args()

if not os.path.isdir(args.outdir):
	print("Output directory does not exist, will create it.")
	os.mkdir(args.outdir)
	
def get_sequences(orthog):
	seq = ""
	for ortho in orthog:
		seq += ">"
		seq += ortho + "\n"
		seq += str(sequences[ortho].seq) +"\n"
	return seq

	
# convert argparse arguments and open files
dumpfile = open(args.dump_path,"r")
ids = args.ids.split(" ")
print("Found", len(ids), "taxon IDs")

# import and read sequence file
sequences =  SeqIO.to_dict(SeqIO.parse(args.fasta_path, "fasta"))

# traverse through dumpfile to select single copy orthologs
i = 1
for orthogroup in dumpfile:
	ortholist = orthogroup.strip().split("\t")
	is_there = []
	for id in ids:
		is_there.append(sum(id in s for s in ortholist))
	if is_there.count(1) == len(ids):
		seqs = get_sequences(ortholist)
		file = open(args.outdir+"/mcl_ortho_"+str(i)+".fas", "w")
		file.write(seqs)
		file.close()
		i += 1
		
print("Wrote ", str(i), "files with single copy genes.")

