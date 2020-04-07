#!/usr/bin/env python
# This script outputs genes for specific annotations:
# USAGE: select_genes_for_annotation.py -gff gff3_files -fasta fasta_files -item CAZy:GH5 > outfile.fa

import argparse
import glob
import sys
from Bio import SeqIO

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="select_genes_for_annotation.py", description = """This script outputs genes for specific annotations""", epilog = """written by Philipp Resl""")
pars.add_argument('-gff', dest="gff3_path", required=True, help="Path to GFF3 files which should be scanned")
pars.add_argument('-fasta', dest="fasta_path", required=True, help="Path to FASTA files containing sequences")
pars.add_argument('-item', dest="item", required=True, help="Annotation which should be retreived")
args=pars.parse_args()


gff3_files = [file for file in glob.glob(args.gff3_path+"/*") if ".gff3" in file]
sys.stderr.write("Found %s GFF3 files.\n" % str(len(gff3_files)))
list_of_ids = []
for handle in gff3_files:
	sys.stderr.write("Reading file: %s\n" % handle)
	file = open(handle, "r")
	i = 0
	for line in file:
		line = line.strip()
		annotation=line.split("\t")[-1]
		annotation = annotation.replace(",", ";")
		name = annotation.split(";")[0].split("=")[-1]
		annotations =  annotation.split(";")[1:]
		for annot in annotations:
			if args.item == annot:
				list_of_ids.append(name)
				i += 1
	sys.stderr.write("%s found %s times.\n" % (args.item, str(i)))

fasta_files = [file for file in glob.glob(args.fasta_path+"/*") if ".fa" in file]
sys.stderr.write("Found %s FASTA files.\n" % str(len(fasta_files)))

for handle in fasta_files:
	sys.stderr.write("Reading file: %s\n" % handle)
	for seq in SeqIO.parse(handle, "fasta"):
		if seq.id in list_of_ids:
			print(">%s" % seq.id)
			print(seq.seq) 
