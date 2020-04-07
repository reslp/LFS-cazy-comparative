# This scripts reroots a phylogenetic tree at the MRCA node of a given set of taxa.

#library(ape)
library(phytools)
library(tidyverse)
#library(ggtree)

args = commandArgs(trailingOnly=TRUE)
treefile <- args[1]
roottaxa <- unlist(strsplit(args[2],","))
#treefile <- "/home/reslp/Dropbox/Philipp/xylographa_comparative_genomics/results/75_genomes/phylogeny/concatenated/concat.treefile"
outfile <- args[3]
print(treefile)
print(roottaxa)
print(outfile)
tree <- read.tree(treefile)
#roottaxa <-  c("Botrytis_cinerea","Erysiphe_necator","Epichloe_typhina","Sordaria_macrospora","Aureobasidium_pullulans","Phyllosticta_citricarpa","Arthonia_radiata")
print(tree)
root <- findMRCA(tree, roottaxa)
tree_new <- ggtree::reroot(tree,root)
tree_out <- ladderize(tree_new, right=T)
tree_out$node.label[1] <- "" # rerooting introduces a new nodename for the created root node. this will get rid of it
ape::write.tree(tree_out, file = outfile)


