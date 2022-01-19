#!/usr/bin/env python

import sys
from Bio import SeqIO

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")
    

name_ids = sys.argv[1]
file_name = sys.argv[2]

ids_file = open(name_ids, "r")
#seq_file = open(file_name, "r")

ids = {}
for line in ids_file:
	line = line.strip()
	ids[line.split("\t")[0]] = line.split("\t")[1] 

for seq in SeqIO.parse(file_name, "fasta"):
	for key in ids.keys():
		if seq.id.startswith(key):
			seq.id = ids[key]
			print(">"+seq.id)
			print(str(seq.seq))
			break
		
	

