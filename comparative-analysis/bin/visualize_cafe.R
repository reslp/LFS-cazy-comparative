#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

wd <- args[1]
treefile <- args[2]
data <- args[3]
plottreefile <- args[4]
prefix <- args[5]

library(ggtree)
library(ggrepel)
library(patchwork)
library(phangorn)
library(phytools)

# variables for testing
#tree <- read.tree(file="../results/52_genomes/cafe/52_genomes_cafe_labeled_tree.tre")
#data <- read.csv(file="../results/52_genomes/cafe/52_genomes_table.txt", header=T, sep="\t")
#plottree <- read.tree(file="../results/52_genomes/phylogeny/52_genomes_ultra.tre")
#setwd("/Users/sinnafoch/Dropbox/Philipp/xylographa_comparative_genomics/test")
#prefix<-"test"

# load all necessary data and set cwd
tree <- read.tree(treefile)
data <- read.csv(file=data, header=T, sep="\t")
plottree <- read.tree(plottreefile)
setwd(wd)
# after reading the CSV we need to create new rownames with node labels
rownames(data) <- data$X
data$X<-NULL

all_nodes <- c(tree$tip.label ,tree$node.label)
data <- data[match(all_nodes, rownames(data)),]
rownames(data)[length(tree$tip.label)+1] <- tree$node.label[1] # root node is not in table so we need to add it manually
print("Extracting significant expansions and contractions...")
sig_expansions <- data$Sig.Expansions[(length(tree$tip.label)+1):length(data$Sig.Expansions)]
sig_contractions <- data$Sig.Contractions[(length(tree$tip.label)+1):length(data$Sig.Contractions)]
print(sig_expansions)
print(sig_contractions)


cafe_tree <- ggtree(ladderize(plottree))+ xlim(NA,3000) +geom_tiplab()
tr <- ggtree(plottree)
d <- tr$data
d <- d[!d$isTip,]
d$expansions <- sig_expansions
d$contractions <- sig_contractions
print(d)
cafe_tree_expansion <- cafe_tree + geom_label_repel(data=d, aes(label=expansions, fill=expansions),size=2, nudge_x=-0.0)+ scale_fill_gradient(low="white", high="#F63E02")+ ggtitle(label="A")+ labs(subtitle="  No. of significantly expanded gene families")
cafe_tree_contraction <- cafe_tree + geom_label_repel(data=d, aes(label=contractions, fill=contractions),size=2, nudge_x=-0.0)+ scale_fill_gradient(low="white", high="#3DA5D9")+ ggtitle(label="B")+ labs(subtitle="  No. of significantly contracted gene families")

pdf(paste(prefix,"_expanded_contracted_families.pdf",sep=""), width=11.7, height=8.3)
cafe_tree_expansion + cafe_tree_contraction
dev.off()


pdf(paste(prefix,"_tree_with_node_numbers.pdf", sep=""), width=8.3, height=11.7)
cafe_tree +geom_label_repel(data=d, aes(label=node), size=2, nudge_x=-1)+ xlim(NA,1500)
dev.off()

indices <- match(tree$node.label, c(tree$tip.label, tree$node.label))
save_df <- data

rownames(save_df) <- c(plottree$tip.label, indices)


write.csv(save_df[1:length(plottree$tip.label),], file=paste(prefix,"_tips_number_table.csv",sep=""))
write.csv(save_df[(length(plottree$tip.label)+1):length(rownames(save_df)),], file=paste(prefix,"_node_number_table.csv",sep=""))

df <- data.frame(cafe=rownames(data),tree=rownames(save_df))
df$tips = ""
for (i in 1:length(rownames(df))) {
  if (is.na(as.numeric(as.character(df$tree[i])))==FALSE) {
    #print(as.numeric(as.character(df$tree[i])))
    df$tips[i] <- as.character(paste(plottree$tip.lab[unlist(Descendants(plottree, as.numeric(as.character(df$tree[i])), "tips"))], collapse=","))
  }
  else {
    df$tips[i] <- as.character(df$tree[i])
  }
}
write.table(df, file=paste(prefix,"_cafe_tree_node_mapping.txt", sep=""), sep="\t", row.names =F, col.names=F, quote=F)


