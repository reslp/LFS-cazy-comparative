#!/usr/bin/env python
# This script creates an input file for r8s given a tree and an alignment (for the number of sites)
# USAGE: create_r8s.py

import argparse
import sys

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="create_r8s.py", description = """This script creates a r8s input file from a template""", epilog = """written by Philipp Resl""")
pars.add_argument('-t', dest="template", required=True, help="Path to template file")
pars.add_argument('-tr', dest="tree", required=True, help="Path to treefile (NEWICK)")
pars.add_argument('-l', dest="logfile", required=False, help="Path to iqtree log file from the tree run.")

args=pars.parse_args()

file = open(args.template, "r")

template = ""
for line in file:
	template += line

#print(template)

file.close()

file = open(args.tree, "r")
tree = file.readline().strip()

#print(tree)

file.close()

file = open(args.logfile, "r")

longest = 0
for line in file:
	if line.startswith("Alignment has "):
		#print(line)
		longest += int(line.split(" ")[-5])
#print(longest)

template = template.replace("%%tree%%", tree)
template = template.replace("%%nsites%%", str(longest))

print(template)
