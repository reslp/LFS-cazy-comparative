library(phytools)
library(ape)
library(ggtree)
library(tidyverse)


args <- commandArgs(trailingOnly=TRUE)

wd <- args[1]
trees_file <- args[2]
outfile <- args[3]

trees <- read.nexus(trees_file)

parse_tiplabels <- function(tree) {
  labels <- tree$tip.label
  new_labels <- c()
  for (n in 1:length(labels)) {
    label <- labels[n]
    if (grepl("*", label, fixed=TRUE) == TRUE) {
      name = strsplit(label, "<")[[1]][1]
      value = strsplit(strsplit(label, "<")[[1]][2], ">")[[1]][2]
      new_labels <- c(new_labels, paste(name, "_", value, sep=""))
    }
    else {
      name = strsplit(label, "<")[[1]][1]
      new_labels <- c(new_labels, name)
    }
  }
  tree$tip.label <- new_labels
  return(tree)
}

parse_nodelabels <- function(tree) {
  labels <- tree$node.label
  new_labels <- c()
  for (n in 1:length(labels)) {
    label <- labels[n]
    if (grepl("*", label, fixed=TRUE) == TRUE) {
      count = strsplit(label, "_")[[1]][2]
      new_labels <- c(new_labels, paste(count, "*", sep=""))
    } else {
      count = strsplit(label, "_")[[1]][2]
      new_labels <- c(new_labels, count)
    }
  }
  tree$node.label <- new_labels
  return(tree)
}

get_gene_names <- function(trees) {
  gene_names <- c()
  for (n in 1:length(names(trees))) {
    name <- tail(strsplit(names(trees)[n],"_")[[1]],1)
    gene_names <- c(gene_names, name)
  }
  return(unique(gene_names))
}



# rename tips and nodes to only keep labels where significant changes occur
trees <- lapply(trees, parse_tiplabels)
class(trees) <- "multiPhylo"
trees <- lapply(trees, parse_nodelabels)
class(trees) <- "multiPhylo" 

# extract all the names of families from the trees:
gene_names <- get_gene_names(trees)

# create dataframe before looping through the trees to get values
occurence_data <- data.frame(matrix(nrow=trees[[1]]$Nnode, ncol=length(gene_names)))
occurence_data[is.na(occurence_data)]<- 0
colnames(occurence_data) <- gene_names

for (n in 1:length(trees)) {
  gene <- tail(strsplit(names(trees)[n],"_")[[1]],1)
  for (j in 1:length(trees[[n]]$node.label)) {
    if (grepl("*",trees[[n]]$node.label[j], fixed=TRUE) == TRUE) {
      occurence_data[j, gene] <- occurence_data[j, gene] + 1
    }
  }
}

# need to find a way to automatically get this number
nmodels <- 22
occurence_data[2,]
combined_names <- c()
for (i in 1:length(rownames(occurence_data))) {
  node_label <- c()
  for (family in colnames(occurence_data)) {
    if (occurence_data[i,family] > 0) {
      node_label <- c(node_label, paste(family,occurence_data[i,family], sep=":"))
    } 
  }
  #if (length(node_label) == 0) {
  #  node_label <- c("0")
  #}
  if (length(node_label) > 0) {
    combined_names <- c(combined_names, paste(node_label, collapse=";"))
  }
  
  #print(occurence_data[i, occurence_data[i,] > 0])
}

#get nodes which have values:
nodes_to_plot <- c()
for (i in 1:length(rownames(occurence_data))){
  if (all(occurence_data[i,] == 0) == FALSE) {
    nodes_to_plot <- c(nodes_to_plot, i+82)
  }
}
nodes_to_plot  

#select the first tree for plotting, they are all the same   
tree <- trees[[1]]
tree$node.label <- combined_names
pdf(file=outfile, width=10, height=10)
ggtree(tree, size=0.3) + geom_tiplab(align=TRUE, size=2.5)+ geom_text2(label=combined_names, aes(subset=(node %in% nodes_to_plot)), nudge_x=-0.03,nudge_y=0.7, size=2) + xlim(0, 200)
dev.off()
