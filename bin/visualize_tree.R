#!/usr/bin/env Rscript
#setwd("/Users/sinnafoch/Dropbox/Philipp/xylographa_comparative_genomics/test")
args <- commandArgs(trailingOnly=TRUE)
print(args)

wd <- args[1]
treefile <- args[2]
outfile <- args[3]

setwd(wd)

library(ggtree)
library(ape)
library(ggrepel)
library(patchwork)
library(phytools)
#library(castor)
#install.packages("castor")
#tree <- castor::read_tree(file="../results/debary_52_genomes/phylogeny/52_genomes.cf.tree",underscores_as_blanks=T) 
tree <- read.tree(file=treefile)
#find_root(tree)
new_root <- findMRCA(tree,c("Sordaria_macrospora", "Arthonia_radiata"))
tree_rooted <- ggtree::reroot(tree, new_root)
#tree_rooted <- ape::root(tree, node=new_root, edgelabels=T, r=T, resolve.root=T)
#tree_rooted <- phytools::reroot(tree, node=new_root)
#tree_rooted <- ladderize(castor::root_at_node(tree, new_root_node=new_root-length(tree$tip.label)))

#length(tree$tip.label)


labels <- tree_rooted$node.label[tree_rooted$node.label !=""]
label_of_second_node <- tree_rooted$node.label[getDescendants(tree_rooted, length(tree$tip.label)+1)[2]-length(tree$tip.label)]
labels <- append(labels, label_of_second_node, after=1)
labels[1] <- ""

bootstrap <- sub("([0-9]*).*", "\\1", labels)
gene_concordance <- sub("[0-9]*/([0-9]*\\.*[0-9]*)/.*", "\\1", labels)
site_concordance <- sub("[0-9]*/[0-9]*\\.*[0-9]*/(.*)", "\\1", labels)

tree_rooted$node.label <- labels
tr <- ggtree(tree_rooted)
d <- tr$data
d <- d[!d$isTip,]
d$bootstrap <- as.numeric(bootstrap)
d$gene_concordance <- as.numeric(gene_concordance)
d$site_concordance <- as.numeric(site_concordance)

tr_bs <- ggtree(ladderize(tree_rooted), aes(label=label)) + geom_treescale(x=0,y=-5, width=0.2) + xlim(0,0.9) +geom_tiplab()
tr_bs <- tr_bs + geom_label_repel(data=d, aes(label=bootstrap, fill=bootstrap),size=2, nudge_x=-0.03) + 
  theme(plot.title=element_text(face="bold"),legend.position=c(0.55,0.05), legend.direction = "horizontal") + scale_fill_viridis_c(begin=0.4) + ggtitle(label="A")
tr_bs <- tr_bs + labs(subtitle="  Phylogeny with bootstrap support")

tr_gc <-  ggtree(ladderize(tree_rooted), aes(label=label)) + geom_treescale(x=0.5,y=-5, width=0.2) + xlim(0,0.9) +geom_tiplab()
tr_gc <- tr_gc + geom_label_repel(data=d, aes(label=gene_concordance, fill=gene_concordance),size=2, nudge_x=-0.03) + 
  theme(plot.title=element_text(face="bold"),legend.position=c(0.55,0.05), legend.direction = "horizontal") + scale_fill_viridis_c(begin=0.4) + ggtitle(label="B")
tr_gc <- tr_gc +labs(subtitle="  Phylogeny with gene concordance support")

tr_sc <-  ggtree(ladderize(tree_rooted), aes(label=label)) + geom_treescale(x=0.5,y=-5, width=0.2) + xlim(0,0.9) +geom_tiplab()
tr_sc <- tr_sc + geom_label_repel(data=d, aes(label=site_concordance, fill=site_concordance),size=2, nudge_x=-0.03) + 
  theme(plot.title=element_text(face="bold"),legend.position=c(0.55,0.05), legend.direction = "horizontal") + scale_fill_viridis_c(begin=0.4) + ggtitle(label="C")
tr_sc <- tr_sc +labs(subtitle="  Phylogeny with site concordance support")
tr_sc

pdf(outfile, width=11.7, height=8.3)
tr_bs + tr_gc + tr_sc
dev.off()


