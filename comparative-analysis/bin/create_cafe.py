#!/usr/bin/env python

import argparse
import glob
import sys
import os
from Bio import SeqIO

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="create_cafe.py", description = """This script creates a CAFE input file""", epilog = """written by Philipp Resl""")
pars.add_argument('-t', dest="tree", required=True, help="Path to the ultrametric tree file")
pars.add_argument('-ct', dest="cafe_tree", required=True, help="Path to the CAFE type tree file with rate regimes.")
pars.add_argument('-template', dest="template", required=True, help="Path to the CAFE input file template")
pars.add_argument('-data', dest="data", required=True, help="File with CAFE formatted gene family information")
pars.add_argument('-pre', dest="pre", required=True, help="Prefix for ouput files in CAFE analysis")
#pars.add_argument('-o', dest="outdir", required=True, help="Output directory")
args=pars.parse_args()

# uncomment only for testing:
#print("Treefile:", args.tree)
#print("CAFE style tree:", args.cafe_tree)
#print("Template:", args.template)
#print("DATA file:", args.data)

file = open(args.tree, "r")
tree = file.readline().strip()
tree = tree.replace(";","")
tree = tree.replace("_","")
file.close()
file = open(args.cafe_tree, "r")
cafe_tree = file.readline().strip()
file.close()

file = open(args.template, "r")

outstring = ""
for line in file:
	if "%%tree%%" in line:
		line = line.replace("%%tree%%", tree)
	if "%%prefix%%" in line:
		line = line.replace("%%prefix%%", args.pre)
	if "%%cafetree%%" in line:
		line = line.replace("%%cafetree%%", cafe_tree)
	if "%%datafile%%" in line:
		line = line.replace("%%datafile%%", args.data)
	outstring += line
		
print(outstring)
		
		


