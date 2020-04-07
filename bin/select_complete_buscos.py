#!/usr/bin/env python
# This script outputs complete sets fo busco sequences
# USAGE: select_complete_buscos.py

import argparse
import glob
import sys
from Bio import SeqIO

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="select_complete_buscos.py", description = """This script selects and renames sets of complete BUSCOs according to the orthology table created by funannotate""", epilog = """written by Philipp Resl""")
pars.add_argument('-tab', dest="tab_path", required=True, help="Path to orthology table created by funannotate compare")
pars.add_argument('-fasta', dest="fasta_path", required=True, help="Path to FASTA files containing sequences")
pars.add_argument('-outpre', dest="outpre", required=True, help="directory for output files")
pars.add_argument('-idfile', dest="idfile", required=True, help="id file for double taxon check")
args=pars.parse_args()

fasta_files = [file for file in glob.glob(args.fasta_path+"/*") if ".fa" in file]
print("Found %s .fa files. (Should be the same as no. of species in dataset)\n" % str(len(fasta_files)))


table_file = open(args.tab_path, "r")

buscos = {}

for line in table_file:
	if line.strip().split("\t")[3] != "None" and len(line.strip().split("\t")[4].split(",")) == len(fasta_files): #filter only complete buscos
		busco_name = line.strip().split("\t")[3].replace(", ", "_") # account for ambiguous naming
		busco_list = line.strip().split("\t")[4]
		buscos[busco_name] = busco_list

sequences = {}
for handle in fasta_files:
	for seq in SeqIO.parse(handle, "fasta"):
		sequences[seq.id] = str(seq.seq)
print("Found %s total sequences." % str(len(sequences)))

idfile = open(args.idfile, "r")

ids = {}
for line in idfile:
	line = line.strip()
	ids[line.split("\t")[0]] = line.split("\t")[1]
#print(ids)


print("Extracting Busco sequences and creating files:")
i = 0
for busco in buscos.keys():
	#print(buscos[busco])
	#check for duplicated sequence ids
	#print (buscos[busco])
	for id in ids.keys():
		flagged = 0
		for busc in buscos[busco].split(", "):
			if busc.startswith(id):
				flagged += 1 # starts with to account for name in names
		if flagged > 1:
			print(busco+" contains duplicated IDs for single species. File will be skipped. Maybe this is worth checking!")
			break
	if flagged == 1:
		# if all is correct go here
		buscos[busco] = buscos[busco].split(", ")
		outfile = open(args.outpre+"/"+busco+".fa", "w")
		i += 1
		for sequence in buscos[busco]:
			outfile.write(">"+sequence+"\n")
			outfile.write(sequences[sequence]+"\n")
		outfile.close()
print ("Created ",i, " files")
print("done")