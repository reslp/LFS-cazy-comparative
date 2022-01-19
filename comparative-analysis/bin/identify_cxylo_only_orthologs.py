#!/usr/bin/env python3
# This script will identify complete sets of orthologs present in all species.
# It uses the orthology_groups.txt file produced by funannotate
# written by Philipp Resl

ortho_file = open("orthology_groups.txt", "r")

species_list_not = ["agyruf", "clamac", "clamet", "cyaast", "dibbae", "evepru", "grascr", "gyafla", "lamins", "lashis", "psefur", "tracoa", "umbmue", "umbpus"]
species_list = ["xylbjo", "xylope", "xylpal", "xylpar", "xylsor", "xyltru"]

for line in ortho_file:
	
	line = line.strip()
	name, dnds, kog, none, species = line.split("\t")
	pres = 0
	for sp in species_list: # check if all species are present
		if sp in species:
			pres += 1
	all_length = len(species.split(",")) 
	abs = 0
	for sp in species_list_not:
		if sp not in species:
			abs += 1
	if pres == len(species_list) and all_length==len(species_list) and abs == len(species_list_not): #final check if all species are present and none is duplicated	
		print(line)
		pres = 0
		abs = 0
	else:
		pres = 0
		abs = 0
	