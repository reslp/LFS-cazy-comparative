#!/usr/bin/env python
# this script extracts sequences from the Fungal virulence database:
# http://sysbio.unl.edu/DFVF/index.php
# this script is crap

file_name = "/Users/sinnafoch/Dropbox/Philipp/xylographa_comparative_genomics/data/database_fungal_virulence.txt"

file = open(file_name, "r")

new_seq = True
for line in file:
	if new_seq == True:
		if line.startswith("UniProtID:"):
			uniprotname = line.split("\t")[-1]
		if line.startswith("Organism:"):
			orgname = line.split("\t")[-1]
			organme = orgname.replace(" ", "_")
		if line.startswith("Protein Sequence:"):
			seq = line.split("\t")[-1].split(" ")[0]
			new_seq = False
	else:
		if line.startswith("\t\t\t\t"):
			seq += line.split("\t")[-1].split(" ")[0]
		if line.startswith("\n"):
			name = ">"+ uniprotname.strip()+ "_"+ orgname.strip()
			name = name.replace(" ", "")
			print(name)
			print(seq)
			new_seq = True
			