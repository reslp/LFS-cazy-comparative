#!/usr/bin/env python

import argparse
import glob
import sys
import os
from Bio import SeqIO
import pandas as pd

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="get_functional_for_cafe_families.py", description = """This script will gather functional annotations for expanded/contracted CAFE families""", epilog = """written by Philipp Resl""")
pars.add_argument('-fam', dest="fam", required=True, help="Path to CAFE family output file. (Created with cafetutorial_report_analysis.py)")
pars.add_argument('-cm', dest="cm", required=True, help="Path to the CAFE tree mapping file. (Created by visualize_cafe.R)")
pars.add_argument('-ortho', dest="ortho", required=True, help="Path to Orthogroups.tsv file (Created by Orthofinder)")
pars.add_argument('-gff', dest="gff", required=True, help="GFF3 type file with annotation information.")
pars.add_argument('-node', dest="type", required=True, help="Select for which nodes information should be extracted. Possible options are [all, internal no. numbers comma seperated]")
pars.add_argument('-pre', dest="pre", required=True, help="Output prefix to all files")
pars.add_argument('-which', dest='which', required=True, help="+ for expanded or - for contracted gene families")
args=pars.parse_args()

whichnodes = ""
if args.type == "all":
	print("Extracting data from all nodes:")
	whichnodes = "all"
elif args.type == "internal":
	print("Extracting data from internal nodes:")
	whichnodes = "internal"
else:
	print("Extracting data for nodes", args.type)
	whichnodes = args.type.split(",")
	
if args.which == "+":
	print("Extracting expanded families")
if args.which == "-":
	print("Extracting contracted families")
	
fam_file = open(args.fam, "r")
#skip first two lines:
fam_file.readline()
fam_file.readline()

expanded_orthogroup_dict = {}
#contracted_orthogroup_dict = {}

for line in fam_file:
	line = line.strip()
	name = line.split("\t")[0]
	name = name.strip(":")
	orthogroups = line.split("\t")[-1]
	orthogroups = orthogroups.split(",")
	if args.which == "+":
		expanded = [ortho.split("[")[0] for ortho in orthogroups if "+" in ortho]
		expanded_orthogroup_dict[name] = expanded
	if args.which == "-":
		expanded = [ortho.split("[")[0] for ortho in orthogroups if "-" in ortho]
		expanded_orthogroup_dict[name] = expanded
fam_file.close()
#print(expanded_orthogroup_dict)
#print(contracted_orthogroup_dict)

cm_file = open(args.cm, "r")
cafe_node_dict = {}
ref_cafe_node_dict = {}
tree_node_dict = {}
for line in cm_file:
	cafe_node_dict[line.strip().split("\t")[0]] = line.strip().split("\t")[1]
	ref_cafe_node_dict[line.strip().split("\t")[1]] = line.strip().split("\t")[0]
	tree_node_dict[line.strip().split("\t")[1]] = line.strip().split("\t")[2].split(",")
cm_file.close()
#print(cafe_node_dict)
#print(tree_node_dict)

orthodata = pd.read_csv(args.ortho, sep="\t")


seq_names_dict = {}
if type(whichnodes) == str:
	for node in expanded_orthogroup_dict.keys():
		if whichnodes == "internal" and not cafe_node_dict[node].isdigit():
			continue
		key = cafe_node_dict[node]
		print("Extracting data for node", key)
		decendent_species = tree_node_dict[key]
		seqnames = []
		for orthogroup in expanded_orthogroup_dict[node]:
			for species in decendent_species:
				li = orthodata.loc[orthodata["Orthogroup"] == orthogroup, [species+".proteins"]].values.tolist()
				a = li[0]
				b = a[0]
				if type(b) is str:
					seqnames.append(b)
		seq_names_dict[key] = ", ".join(seqnames)
			
elif type(whichnodes) == list:
	for node in whichnodes:
		print("Extracting data for node", node)
		key = ref_cafe_node_dict[node]
		decendent_species = tree_node_dict[node]
		seqnames = []
		for orthogroup in expanded_orthogroup_dict[key]:
			for species in decendent_species:
				li = orthodata.loc[orthodata["Orthogroup"] == orthogroup, [species+".proteins"]].values.tolist()
				a = li[0]
				b = a[0]
				if type(b) is str:
					seqnames.append(b)
				#print(a)
			#print(",".join(seqnames))
		seq_names_dict[node] = ", ".join(seqnames)
		
def extract_interpro(ann_str):
	iplist = []
	for annot in ann_str.split(";"):
		if "InterPro" in annot:
			for interpro in annot.split(","):
				if "InterPro" in interpro.split(":")[0]:
					iplist.append(interpro.split(":")[-1])
	return iplist

def extract_pfam(ann_str):
	pflist = []
	for annot in ann_str.split(";"):
		if "PFAM" in annot:
			for pfam in annot.split(","):
				if "PFAM" in pfam.split(":")[0]:
					pflist.append(pfam.split(":")[-1])
	return pflist

def extract_eggnog(ann_str):
	egglist = []
	for annot in ann_str.split(";"):
		if "EggNog" in annot:
			for egg in annot.split(","):
				if "EggNog" in egg.split(":")[0]:
					egglist.append(egg.split(":")[-1])
	return egglist

def extract_product(ann_str):
	prodlist = []
	for annot in ann_str.split(";"):
		if "product" in annot:
			for prod in annot.split(","):
				if "product" in prod.split("=")[0]:
					prodlist.append(prod.split("=")[-1])
	return prodlist

def extract_cazy(ann_str):
	cazylist = []
	for annot in ann_str.split(";"):
		if "CAZy" in annot:
			for cazy in annot.split(","):
				if "CAZy" in cazy.split(":")[0]:
					cazylist.append(cazy.split(":")[-1])
	return cazylist

def extract_goterms(ann_str):
	golist = []
	for annot in ann_str.split(";"):
		if "GO_function" in annot:
			for go in annot.split(","):
				if "GO_function" in go.split(": ")[0]:
					golist.append(go.split(": ")[-1])
	#print(len(golist))
	return golist

class node_functional:
	def __init__(self, node):
		self.nnode = node
		self.interprolist = []
		self.plist = []
		self.egglist = []
		self.prodlist = []
		self.cazylist = []
		self.golist = []	
		self.totalseq = 0
		self.seqcount = 0
		self.seqnames = []
		
	def stats(self):
		stat = str(self.nnode)+"\t"+ str(len(self.seqnames))+ "\t"+ str(self.seqcount)+"\t"+ str(len(self.interprolist))+ "\t"+ str(len(set(self.interprolist)))+"\t"+ str(len(self.plist))+ "\t"+ str(len(set(self.plist)))+"\t"+str(len(self.egglist))+"\t"+str(len(set(self.egglist)))+"\t"+ str(len(self.prodlist))+ "\t"+ str(len(set(self.prodlist)))+ "\t"+ str(len(self.cazylist))+ "\t"+ str(len(set(self.cazylist)))+ "\t"+ str(len(self.golist))+ "\t"+ str(len(set(self.golist)))
		return stat

gff_file = open(args.gff, "r")
total_node_list = []
for key in seq_names_dict.keys():
	print("Parsing annotations for node:", key)
	this_node = node_functional(key)
	seqnames = seq_names_dict[key]
	seqnamelist = seqnames.split(", ")
	this_node.totalseq = len(seqnamelist)
	seqcount = 0
	for seq in seqnamelist:
		for line in gff_file:
			if seq in line:
				line = line.strip("\n")
				annotation = line.split("\t")[-1]
				flag = 0
				if "InterPro" in annotation:
					l = extract_interpro(annotation)
					iprann = "\n".join(l)
					this_node.interprolist += l
					#print(iprann)
					flag = 1
				if "PFAM" in annotation:
					pl =extract_pfam(annotation)
					pfann = "\n".join(pl)
					this_node.plist += pl 
					flag = 1
				if "EggNog" in annotation:
					eggl = extract_eggnog(annotation)
					eggann = "\n".join(eggl)
					this_node.egglist += eggl
					#print(eggann)
					flag = 1
				if "product" in annotation:
					prodl = extract_product(annotation)
					prodann = "\n".join(prodl)
					this_node.prodlist += prodl
					#print(prodann)
					flag = 1
				if "CAZy" in annotation:
					cazyl = extract_cazy(annotation)
					cazyann = "\n".join(cazyl)
					this_node.cazylist += cazyl
					#print(cazyann)
					flag = 1
				if "GO_function" in annotation:
					gol = extract_goterms(annotation)
					goann = "\n".join(gol)
					this_node.golist += gol
					#print(goann)
					flag = 1
				if flag == 1:
					this_node.seqcount += 1
				this_node.seqnames.append(seq)
				break	
		gff_file.seek(0)
	total_node_list.append(this_node)

print("Writing output files:")
if args.which =="+":
	args.pre += "_expanded"
if args.which == "-":
	args.pre += "_contracted"
statsfile = open(args.pre+"_statistics.txt", "w")
statsfile.write("node"+"\t"+"genes"+"\t"+"annot.genes"+"\t"+"InterproIDs"+"\t"+"unique InterproIDs"+"\t"+"PFAMs"+"\t"+"unique PFAMs"+"\t"+"EggNogs"+"\t"+"unique Eggnogs"+"\t"+"productnames"+"\t"+"unique names"+"\t"+"CAZy"+"\t"+"unique CAZy"+"\t"+"GOterms"+"\t"+"unique GOterms\n")
for node in total_node_list:
	statsfile.write(node.stats()+"\n")	
statsfile.close()

all_iprs = []
for node in total_node_list:
	for ip in set(node.interprolist):
		if ip not in all_iprs:
			all_iprs.append(ip)
	#print(len(set(node.interprolist)))

interprofile = open(args.pre+"_interpro.txt", "w")
interprofile.write("node"+"\t"+"\t".join(all_iprs)+"\n")
for node in total_node_list:
	interprofile.write(str(node.nnode)+"\t"+"\t".join([str(node.interprolist.count(ip)) for ip in all_iprs])+"\n")	
interprofile.close()

all_pfam = []
for node in total_node_list:
	for pf in set(node.plist):
		if pf not in all_pfam:
			all_pfam.append(pf)
	#print(len(set(node.plist)))

pfamfile = open(args.pre+"_pfam.txt", "w")
pfamfile.write("node"+"\t"+"\t".join(all_pfam)+"\n")
for node in total_node_list:
	pfamfile.write(str(node.nnode)+"\t"+ "\t".join([str(node.plist.count(x)) for x in all_pfam])+"\n")	
pfamfile.close()

all_egg = []
for node in total_node_list:
	for x in set(node.egglist):
		if x not in all_egg:
			all_egg.append(x)
	#print(len(set(node.egglist)))

eggfile = open(args.pre+"_egg.txt", "w")
eggfile.write("node"+"\t"+"\t".join(all_egg)+"\n")
for node in total_node_list:
	eggfile.write(str(node.nnode)+"\t"+ "\t".join([str(node.egglist.count(x)) for x in all_egg])+"\n")	
eggfile.close()

all_prod = []
for node in total_node_list:
	for x in set(node.prodlist):
		if x not in all_prod:
			all_prod.append(x)
	#print(len(set(node.prodlist)))

prodfile = open(args.pre+"_prod.txt", "w")
prodfile.write("node"+"\t"+"\t".join(all_prod)+"\n")
for node in total_node_list:
	prodfile.write(str(node.nnode)+"\t"+ "\t".join([str(node.prodlist.count(x)) for x in all_prod])+"\n")	
prodfile.close()

all_cazy = []
for node in total_node_list:
	for x in set(node.cazylist):
		if x not in all_cazy:
			all_cazy.append(x)
	#print(len(set(node.cazylist)))

cazyfile = open(args.pre+"_cazy.txt", "w")
cazyfile.write("node"+"\t"+"\t".join(all_cazy)+"\n")
for node in total_node_list:
	cazyfile.write(str(node.nnode)+"\t"+ "\t".join([str(node.cazylist.count(x)) for x in all_cazy])+"\n")	
cazyfile.close()

all_go = []
for node in total_node_list:
	for x in set(node.golist):
		x = x.split(" ")[0]
		if x not in all_go:
			all_go.append(x)
	#print(len(set(node.golist)))

gofile = open(args.pre+"_go.txt", "w")
gofile.write("node"+"\t"+"\t".join(all_go)+"\n")
for node in total_node_list:
	simple_go_list =[go.split(" ")[0] for go in node.golist]
	gofile.write(str(node.nnode)+"\t"+ "\t".join([str(simple_go_list.count(x)) for x in all_go])+"\n")	
gofile.close()


	