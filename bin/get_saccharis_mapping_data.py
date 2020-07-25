#!/usr/bin/env python

import argparse
import sys
import pandas as pd

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="get_saccharis_mapping_data.py", description = """This script uses multiple files to create a file used to plot information on saccharis gene trees""", epilog = """written by Philipp Resl""")
pars.add_argument('-fasta', dest="fasta", required=True, help="FASTA file containing the renamed sequences included in saccharis.")
pars.add_argument('-info', dest="info", required=True, help="Genome information summary csv file")
pars.add_argument('-map', dest="map", required=True, help="Addditional mapping file containing information per species about lifestyle etc.")
args=pars.parse_args()


# extract fasta sequence names
print("Processing FASTA file", file=sys.stderr)
file = open(args.fasta, "r")

sequences = {}
for line in file:
	if line.startswith(">"):
		saccharis_name, transcript_name = line.strip().split(" ")
		sequences[transcript_name] = saccharis_name.split(">")[1]
file.close()


# load genome information file
print("Processing Genome info file", file=sys.stderr)
data = pd.read_csv(args.info, sep=";")

species_names = data["isolate"].to_list()
group = data["class"].to_list()
tags = data["locus_tag"].to_list()

tag_species = {}
for tag, sp in zip(tags, species_names):
	tag_species[tag] = sp

species_group = {}
for sp, gr in zip(species_names, group):
	species_group[sp] = gr


# load additional mapping file
print("Processing additional mapping file", file=sys.stderr)
data = pd.read_csv(args.map, sep=",")

species_names = data.iloc[:,0].to_list()
lichen = data["lichen"].to_list()

species_lichen = {}
for sp, lich in zip(species_names, lichen):
	species_lichen[sp] = lich
	
print("Writing output", file=sys.stderr)
print("saccharis_name\tgene_name\tspecies_name\tclass\tlichen")

for name in sequences.keys():
	outstr = sequences[name] + "\t" + name + "\t"
	for lt in tag_species.keys():
		if lt in name:
			outstr += tag_species[lt]
			outstr += "\t"
			outstr += species_group[tag_species[lt]]
			outstr += "\t"
			outstr += str(species_lichen[tag_species[lt]])
			break
	print(outstr)
	
