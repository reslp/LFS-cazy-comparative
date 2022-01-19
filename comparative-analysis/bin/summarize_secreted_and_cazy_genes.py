#!/usr/bin/env python
# This script selects secreted genes (predicted with signal p) and cazymes
# USAGE: select_secreted_and_cazy_genes.py -gff gff3_files -fasta fasta_files

import argparse
import glob
import sys

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="select_genes_for_annotation.py", description = """This script outputs crosslinked information on secreted enzymes and cazymes.""", epilog = """written by Philipp Resl""")
pars.add_argument('-gff', dest="gff3_path", required=True, help="Path to GFF3 files which should be scanned")
args=pars.parse_args()


gff3_files = [file for file in glob.glob(args.gff3_path+"/*") if ".gff3" in file]
print("Found %s GFF3 files.\n" % str(len(gff3_files)), file=sys.stderr)


n_secreted = 0
n_cazy = 0
n_transmembrane = 0
n_secreted_cazy = 0
n_transmem_cazy = 0
n_cazy_internal = 0

print("species\tn_cazy\tn_secretd\tn_transmembrane\tn_secreted_cazy\tn_transmem_cazy\tn_cazy_internal")
for handle in gff3_files:
	n_secreted = 0
	n_cazy = 0
	n_transmembrane = 0
	n_secreted_cazy = 0
	n_transmem_cazy = 0
	n_cazy_internal = 0
	sp = handle.split("/")[-1].split(".")[0]
	file = open(handle, "r")
	i = 0
	for line in file:
		line = line.strip()
		if "SECRETED" in line:
			n_secreted += 1
		if "CAZy" in line:
			n_cazy += 1
			if "TransMembrane" not in line and "SECRETED" not in line:
				n_cazy_internal += 1
		if "TransMembrane" in line:
			n_transmembrane += 1
		if "SECRETED" in line and "CAZy" in line:
			n_secreted_cazy += 1
		if "TransMembrane" and "CAZy" in line:
			n_transmem_cazy += 1
	print(sp, n_cazy, n_secreted, n_transmembrane, n_secreted_cazy, n_transmem_cazy, n_cazy_internal, sep="\t")

#fasta_files = [file for file in glob.glob(args.fasta_path+"/*") if ".fa" in file]
#sys.stderr.write("Found %s FASTA files.\n" % str(len(fasta_files)))

#for handle in fasta_files:
#	sys.stderr.write("Reading file: %s\n" % handle)
#	for seq in SeqIO.parse(handle, "fasta"):
#		if seq.id in list_of_ids:
#			print(">%s" % seq.id)
#			print(seq.seq)
