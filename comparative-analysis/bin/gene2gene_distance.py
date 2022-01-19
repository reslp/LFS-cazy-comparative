#!/usr/bin/env python
# this script will calculate the distribution of gene to gene distances from a GFF3 file.

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from statistics import median
from os.path import basename
import sys
import argparse

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="gene2gene_distance.py", description = """This will calculate gene to gene distances from funannotate annotation files.""", epilog = """written by Philipp Resl""")
pars.add_argument('-i', dest="input_files", required=True,nargs="+",action="store", help="Input file name(s)")
args=pars.parse_args()

#print(args.input_files)


outfile = open("gene_length_species_medians.txt", "w")
for file in args.input_files:
	data = pd.read_csv(file, sep="\t", header=0)
	start_stop_data = data["scaffold:start-end"].tolist()
	contigs = set([cont.split(":")[0] for cont in start_stop_data])
	starts = [int(std.split(":")[-1].split("-")[0]) for std in start_stop_data]
	stops = [int(std.split(":")[-1].split("-")[1]) for std in start_stop_data]
	distances = []
	for contig in contigs:
		for start, stop in zip(starts, stops):
			distances.append(stop - start)
	#df = pd.DataFrame(distances)
	#plt.figure(figsize=(10,7), dpi= 80)
	#title = basename(file)+"\nMedian: %s bp" % str(median(distances))
	#sns_plot = sns.distplot(df , color="dodgerblue", axlabel="gene length").set_title(title)
	#fig = sns_plot.get_figure()
	#fig.savefig(basename(file).split(".")[0]+"_gene_length.pdf")
	#plt.close()
	outfile.write(basename(file).split(".")[0]+"\t"+str(median(distances))+"\n")
outfile.close()

outfile = open("gene2gene_species_medians.txt", "w")
for file in args.input_files:
	data = pd.read_csv(file, sep="\t", header=0)
	start_stop_data = data["scaffold:start-end"].tolist()
	contigs = set([cont.split(":")[0] for cont in start_stop_data])
	starts = [int(std.split(":")[-1].split("-")[0]) for std in start_stop_data]
	stops = [int(std.split(":")[-1].split("-")[1]) for std in start_stop_data]

	distances = []
	for contig in contigs:
		for i in range(0,len(starts)-1):
			distances.append(starts[i+1]-stops[i])
	#df = pd.DataFrame(distances)
	#plt.figure(figsize=(10,7), dpi= 80)
	#title = basename(file)+"\nMedian: %s bp" % str(median(distances))
	#sns_plot = sns.distplot(df , color="dodgerblue", axlabel="gene2gene distance").set_title(title)
	#fig = sns_plot.get_figure()
	#fig.savefig(basename(file).split(".")[0]+"_gene2gene_distance.pdf")
	#plt.close()
	outfile.write(basename(file).split(".")[0]+"\t"+str(median(distances))+"\n")
outfile.close()
