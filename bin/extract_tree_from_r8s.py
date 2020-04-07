#!/usr/bin/env python
# USAGE: extract_tree_from_r8s.py

import argparse
import sys

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="extract_tree_from_r8s.py.py", description = """Extracts tree from a r8s output file""", epilog = """written by Philipp Resl""")
pars.add_argument('-r', dest="r8s", required=True, help="Path to r8s output file")

args=pars.parse_args()

file = open(args.r8s, "r")

tree = ""
for line in file:
		tree = line.strip().split("tree tree_1 =")[-1]
		tree = tree.split(";")[0]
print(tree+";")