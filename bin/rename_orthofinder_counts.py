#!/usr/bin/env python
# written by Philipp Resl

import sys
import pandas as pd
import argparse

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="rename_orthofinder_counts.py", description = """Script will rename Orthofinder counts to more mmeaningful names based on the presence of characterized sequences present in each orthogroup""", epilog= """written by Philipp Resl""")
pars.add_argument('--orthogroups-file', dest="og", required=True, help="Orthogroups.txt file from Orthofinder")
pars.add_argument('--counts-file', dest="counts", required=True, help="Orthogroups.GeneCount.tsv file from Orthofinder")
pars.add_argument('--mapping-file', dest="mapping", required=True,  help="File with mapping IDs (2 column tab separated). Column 1 should be sequence names of file used in Orthfinder. Column 2 should be meaningful name of sequence.")
args=pars.parse_args()


orthogroups_file = args.og
#ids_file = sys.argv[2]
mapping_file = args.mapping
orthogroup_counts_file = args.counts

#ids = {}
#with open(ids_file, "r") as idf:
#	for line in idf:
#		line = line.strip()
#		ids[line.split("\t")[0]] = line.split("\t")[1]
mapping = {}
with open(mapping_file, "r") as mf:
	for line in mf:
		line = line.strip()
		mapping[line.split("\t")[0]] = line.split("\t")[1]

orthogroups = {}
with open(orthogroups_file, "r") as ogf:
	for line in ogf:
		line = line.strip()
		orthogroups[line.split(":")[0]] = [g for g in line.split(":")[1].split(" ") if g]


def get_mapped_name(og):
	oglist = orthogroups[og]
	ogset = set(oglist)
	name = ""
	for gene in oglist:
		if gene in mapping.keys():
			if name:
				if mapping[gene] not in name:
					name += ";" + mapping[gene]
			else:
				name = mapping[gene]
	return name

with open(orthogroup_counts_file, "r") as ogcf:
	print("Orthogroup_old\t", ogcf.readline().strip())
	for line in ogcf:
		line = line.strip()
		og = line.split("\t")[0]
		name = get_mapped_name(og)
		if name:
			line = line.replace(og, name)	
		line = og + "\t" + line
		print(line)
#print(mapping)
#def get_counts(genes):
#	elements = set(genes)
#	counts = {}
#	for element in elements:
#		counts[element] = genes.count(element)
#	return counts
#
##print(ids)
#orthogroups = {}
#with open(orthogroups_file, "r") as ogf:
#	for line in ogf:
#		line = line.strip()
#		# get gene list
#		genes = [g for g in line.split(":")[1].split(" ")]
#		# remove empty elements
#		genes = [g for g in genes if g]
#		# reduce gene names
#		genes = [g.split("pred")[0] for g in genes]
#		orthogroups[line.split(":")[0]] = get_counts(genes)
#		
#		#print(genes)
#
#
#counts_dict ={key:{ids[k] if k in ids.keys() else k:v for (k,v) in value.items()} for (key,value) in orthogroups.items()}
#
#
#counts_renamed = {}
#for key in counts_dict.keys():
#	print(key)
#	for name in counts_dict[key].keys():
#		print(name)
#		if name in mapping.keys():
#			counts_renamed[mapping[name]] = counts_dict[key]
#		else:
#			counts_renamed[key] = counts_dict[key]
#
#
#print(list(ids.values()))
#counts_renamed = {key:{k:v for (k,v) in value.items() if k in list(ids.values())} for (key,value) in counts_renamed.items()}
#
##counts_renamed = {mapping[key] if key in mapping.keys() else key:value for (key,value) in counts_dict.items()}
#
#df = pd.DataFrame.from_dict(counts_renamed, orient="index")
#
#df.to_csv("counts_test.csv")
#	
