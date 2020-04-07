#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

wd <- args[1]
treefile <- args[2]
datafile <- args[3]
prefix <- args[4]
iprfile <- args[5]
pffile <- args[6]
#which <- args[7]
#type <- args[8]
outdir <- args[7]

print(wd)
print(treefile)
print(datafile)
print(prefix)
print(iprfile)
print(pffile)

# get the functional term which should be visualized from the filename
which <- tail(unlist(strsplit(datafile, "_")), n=1)
which <- head(unlist(strsplit(which, "\\.")), n=1)
print(which)
# get whether expansions or contractions should be analysed
type <- as.vector(unlist(strsplit(datafile, "_")))
item <- length(type)-1
type <- type[item]
print(type)
print(outdir)

library(reshape2)
library(ggplot2)
library(ggrepel)
library(ggtree)
library(phytools)
library(patchwork)
library(cowplot)
setwd(wd)
iprmapping <- read.csv(iprfile, header=T, sep="\t")
pfammapping <- read.csv(pffile, sep="\t", header=F)
tree <- read.tree(treefile)

# settings for testing
#setwd("/Users/sinnafoch/Dropbox/Philipp/xylographa_comparative_genomics/test")
#iprmapping <- read.csv("../data/52_genomes/entry.list", header=T, sep="\t")
#eggmapping <- read.csv("../data/52_genomes/fuNOG.annotations.tsv", sep="\t", header=F)
#pfammapping <- read.csv("../data/52_genomes/pfam_all_description.txt", sep="\t", header=F)
#tree <- read.tree("../results/52_genomes/phylogeny/52_genomes_ultra.tre")
#prefix <- "test"
#types <- c("pfam","interpro", "egg", "cazy", "go", "prod")

plot_stats <- function(data, type) {
  for (i in 1:nrow(data)) {
    node <- data[i, ] #select node of interest
    name <- node[[1]]
    print(name)
    ipr <- node[2:length(node)]
    df <- melt(ipr)
    df <- df[df$value != 0,] # remove entries with zero counts
    df <- df[sample(1:nrow(df)),] # shuffle rows for better spread and labeling
    val=40
    if (length(df$value[order(df$value, decreasing=TRUE)])<40) {
      val = length(df$value[order(df$value, decreasing=TRUE)])
    } else {val = 40}
    
    cutoff <- df$value[order(df$value, decreasing=TRUE)[val]] # cutoff for 50 highest values which will be labeled
    df$variable <- factor(df$variable, levels = df$variable) # make factors so that ggplot does not reorder
    p <- ggplot(data=df, aes(x=variable, y=value, label=variable)) + geom_point(color = ifelse(df$value >= cutoff, "red", "black"), stat="identity") + theme_minimal() +
      theme(axis.line = element_line(colour = "black"),axis.ticks.x = element_blank(),axis.text.x  = element_blank(),axis.title.x = element_blank(),plot.title= element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank())+  ylab(paste("No. of ",type,sep=""))+# ylim(-4,max(df$value))+
      geom_text_repel(data= subset(df, value >= cutoff),hjust= 0,nudge_y=0.5,segment.size = 0.1, size=2, max.iter=5000, ylim=c(cutoff-3,NA), force=10, direction="y") +
      ggtitle(paste(type," in genes of ", which, " gene families at node ", name, sep=""))
    ptr <- ggtree(tree, size=0.2) +geom_point2(aes(subset=(node==name)), size=2, shape=23, fill="steelblue")
    pdf(file=paste(outdir,prefix,"_",type,"_",which,"_node_",name,".pdf",sep=""), width=11.7, height=6)
    print(p+ {ptr+ plot_spacer()+plot_layout(ncol = 1, heights=c(1,2))}+plot_layout(widths=c(5,1)))
    dev.off() 
  }
}


data <- read.csv(datafile, sep="\t", header=T)
if (type == "interpro") {
    names <- vector()
    for (i in 1:length(colnames(data))) {
      names <- c(names, as.vector(iprmapping$ENTRY_NAME[iprmapping$ENTRY_AC == colnames(data)[i]]))
    }
    names <- c("node", names)
    colnames(data) <- names
} else if (type == "pfam") {
    names <- vector()
    for (i in 1:length(colnames(data))) {
      ssub <- sub("[.].*", "", as.character(pfammapping$V1))
      names <- c(names, as.vector(pfammapping$V2[ssub == colnames(data)[i]]))
    }
    names <- c("node", names)
    colnames(data) <- names
} 
plot_stats(data, type)

