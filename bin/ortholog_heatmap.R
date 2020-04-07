#!/usr/bin/env Rscript
#setwd("/Users/sinnafoch/Dropbox/Philipp/xylographa_comparative_genomics/test")
args <- commandArgs(trailingOnly=TRUE)
print(args)

wd <- args[1]
treefile <- args[2]
datafile <- args[3]
prefix <- args[4]
set <- args[5]

print(wd)
print(treefile)
print(datafile)
print(prefix)
print(set)

#if (!"patchwork" %in% installed.packages())
#{
#  options(unzip = "internal")
#  devtools::install_github("thomasp85/patchwork", upgrade=F)
#}


library(ape)
library(ggplot2)
library(reshape2)
library(patchwork)
library(ggtree)
library(wesanderson)
library(tidyverse)


tree <- read.tree(treefile)
data <- read.table(datafile)
#data <- read.table("../results/52_genomes/orthofinder/orthofinder/Results_Jun20/Comparative_Genomics_Statistics/OrthologuesStats_one-to-many.tsv")
#data <- read.table("../results/52_genomes/orthofinder/orthofinder/Results_Jun20/Comparative_Genomics_Statistics/Orthogroups_SpeciesOverlaps.tsv")

rownames(data) <- gsub("\\.proteins","\\1",rownames(data))
colnames(data) <- gsub("\\.proteins","\\1",rownames(data))

df <- data[,rev(match(tree$tip.label, colnames(data)))]
#df_melted <- melt(round(cor(df), 2)) # create a correlation map
df_melted <- melt(as.matrix(df))
head(df_melted)
colnames(df_melted) <- c("label", "category", "no.of.OG")
midp <- mean(df_melted$no.of.OG)
#function definitions. These functions have been modified from:
# https://thackl.github.io/ggtree-composite-plots

tree_y <-  function(ggtree, data){
  if(!inherits(ggtree, "ggtree"))
    stop("not a ggtree object")
  left_join(select(data, label), select(ggtree$data, label, y)) %>%
    pull(y)
}

# get the range of the ggtree y-axis data
tree_ylim <- function(ggtree){
  if(!inherits(ggtree, "ggtree"))
    stop("not a ggtree object")
  range(ggtree$data$y)
}


# plot data next to a ggtree aligned by shared labels
ggtreeplot <- function(ggtree, data = NULL, mapping = aes(), flip=FALSE,
                       expand_limits=expand_scale(0,.6), ...){
  if(!inherits(ggtree, "ggtree"))
    stop("not a ggtree object")
  
  # match the tree limits
  limits <- tree_ylim(ggtree)
  limits[1] <- limits[1] + (limits[1] * expand_limits[1]) - expand_limits[2]
  limits[2] <- limits[2] + (limits[2] * expand_limits[3]) + expand_limits[4]
  
  if(flip){
    mapping <- modifyList(aes_(x=~x), mapping)
    data <- mutate(data, x=tree_y(ggtree, data))
    gg <- ggplot(data=data, mapping = mapping, ...) +
      scale_x_continuous(limits=limits, expand=c(0,0))
  }else{
    print("here")
    mapping <- modifyList(aes_(y=~y), mapping)
    data <- mutate(data, y=tree_y(ggtree, data))
    gg <- ggplot(data=data, mapping = mapping, ...) +
      scale_y_continuous(limits=limits, expand=c(0,0))
  }
  gg
}

no_legend <- function() theme(legend.position="none")

scale_y_tree <- function(expand=expand_scale(0, 0.6), ...){
  scale_y_continuous(expand=expand, ...)
}

no_y_axis <- function () 
  theme(axis.line.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        )
no_x_axis <- function () 
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=315, hjust = 0, size=8))

ptree <- ggtree(tree)+ geom_tiplab(align=T, size=2.7)+ scale_x_continuous(expand=expand_scale(0.7))+scale_y_tree()


pal <- wes_palette("Zissou1", length(df_melted$no.of.OG), type = "continuous")
pheatmap <- ggtreeplot(ptree, data=df_melted, aes(x=category)) + geom_tile(aes(fill=no.of.OG ), color="white")+
      scale_fill_gradientn(colors = pal) +no_x_axis() + no_y_axis()+
      ggtitle(paste("No. of Orthogroups ", set, sep=""))


pdf(file=paste(prefix,"_",set, ".pdf", sep=""), width=11.7, height=8.3)
ptree + pheatmap
dev.off()


#save.image(file=paste(prefix,"_",set,"_orthology_stats.RData",sep=""))