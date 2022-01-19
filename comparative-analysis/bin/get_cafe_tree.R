#!/usr/bin/env Rscript

# this script will reformat a newick tree file string to a pseudo-newick string
# which CAFE understands to impose different rate regimes

args <- commandArgs(trailingOnly=TRUE)
print(args)

wd <- args[1]
treefile <- args[2]
outfile1 <- args[3]
outfile2 <- args[4]
print(length(args))
if (length(args) == 5)
{
  species <- unlist(strsplit(args[5],","))
  print(species)
} else {species <- ""}

setwd(wd)
## variables for testing
#setwd("/Users/sinnafoch/Dropbox/Philipp/xylographa_comparative_genomics/test")
#treefile <- "../results/52_genomes/phylogeny/52_genomes_ultra.tre"
#species <- c("Xylographa_parallela", "Xylographa_trunciseda")
#outfile <- "test_tree.tre"

#print(wd)
#print(treefile)
#print(species)
#print(length(species))
#print(outfile1)
#print(outfile2)

library(phytools)
tree <- read.tree(file=treefile)
tree$node.label <- NULL # we don't need support values or other node labels for this
tree_names <- tree # save a copy of the tree for later
tree$edge.length = tree$edge.length / 10
outstring1 <- write.tree(tree)
write(outstring1, file=outfile1)

replstr <- function(str) {
  return("")
}

#tree$tip.label <- unlist(lapply(tree$tip.label, replstr))

if (length(species) == 0)  { # for single rates (no region of interest is specified in the tree)
  print("Single rate analysis")
  for (i in 1:length(tree$edge.length)) {
    tree$edge.length[i] = 1
  }
} else { # two rates for different parts of the tree
  print("Two rate analysis")
  print(tree)
  node <- getMRCA(tree_names, as.vector(species))
  decendent_nodes <- getDescendants(tree, node)
  decendent_nodes <- c(node, decendent_nodes) # add the original node to decendents to get the MRCA edge
  indices <- match(decendent_nodes, tree$edge[,2])
  
  for (i in 1:length(tree$edge.length)) {
    if (i %in% indices) {
      tree$edge.length[i] = 2
    }
    else {tree$edge.length[i] = 1}
  }
}



# needs to do some reformatting because the needed tree is not really newick format
outstring <- write.tree(tree)
outstring <- gsub(":","", outstring)
outstring <- gsub(";","", outstring)
write(outstring, file=outfile2)






