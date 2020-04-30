#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

wd <- args[1]
cazy_file <- args[2]
treefile <- args[3]
prefix <- args[4]
set <- args[5]
outdir <- args[6]

print(wd)
print(cazy_file)
print(treefile)
print(prefix)
print(set)
print(outdir)

#if (!"patchwork" %in% installed.packages())
#{
#  devtools::install_github("thomasp85/patchwork", upgrade=F)
#}
#installed.packages()
library(phytools)
library(ggplot2)
library(ggpubr)
library(ggtree)
library(cowplot)
library(reshape2)
library(patchwork)
library(tidyverse)




#wd <- "/Users/sinnafoch/Dropbox/Philipp/xylographa_comparative_genomics/test"
#prefix <- "hemi"
#set <- "hemicellulose"



setwd(wd)
#setwd('/Volumes/sinnafoch/sinnafoch/Dropbox/Philipp/Genomes/02_cazy_distribution/anc')

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

no_y_axis <- function () 
  theme(axis.line.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
no_x_axis <- function () 
  theme(axis.line.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

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

# reformat the tree to ages between 0 and 1:
# Loading tree:

print("Loading tree...")
tree <- read.tree(treefile)
mbt <- max(branching.times(tree))
tree$edge.length <- tree$edge.length / mbt

ptree <- ggtree(tree, size=0.3) + geom_tiplab(align=TRUE, size=2.5) + geom_treescale(x=0.1,y=50, width=0.2, linesize=0.3) + scale_x_continuous(expand=expand_scale(0.2)) + scale_y_tree() +xlim(NA, 1.3) #+geom_treescale(x=0,y=45, width=0.5)#+scale_x_continuous(expand=expand_scale(0.2))

d <- fortify(tree)
d <- subset(d, isTip)
tip_order <- with(d, label[order(y, decreasing=T)])


print("Read cazy summary input file...")
all_cazy <- read.csv(file=cazy_file)
rownames(all_cazy) <- all_cazy$X
all_cazy$X <- NULL
all_cazy <- t(all_cazy)

all_cazy <- all_cazy[match(tip_order, rownames(all_cazy)),]
all_cazy <- all_cazy [nrow(all_cazy ):1,]


#change working directory for output
print("Change WD for output...")
print(paste(wd,"/results/",prefix,"/",outdir,"/",sep=""))
setwd(paste(wd,"/results/",prefix,"/",outdir,"/",sep=""))

print("plotting tree...")
tree_out <- paste(prefix, "_tree.pdf", sep="")
pdf(file=tree_out, height=11.3, width=7.8) 
plot(tree, edge.width=2, cex=1.2, adj=0)
nodelabels()
add.scale.bar()
ptree
dev.off()

#Ancestral State reconstruction with Cellulose degrading set of genes according to Floudas 2012
if (set == "cellulose") {
  cellulose_df <- all_cazy[,colnames(all_cazy) %in% c("GH3","GH5","GH6","GH7","GH10","GH11", "GH12", "GH28","GH43","GH45","GH61","GH74", "CBM1", "AA9", "CE1", "CE16", "CE5", "CE8", "CE12", "CE15")]
  print("Will analyze the cellulose set")
} else if (set == "hemicellulose") {
  cellulose_df <- all_cazy[,colnames(all_cazy) %in% c("GH10", "GH11", "CE1", "CE12", "CE3", "CE2", "CE5", "CE7", "GH43_29", "GH43_34", "GH51", "GH45", "GH62")]
  print("Will analyze the hemicellulose set")
} else {
  cellulose_df <- all_cazy
  print("Will analyze all cazymes")
}

print("Reconstructing ancestral states...")
fit_cel <- list()
obj_cel <- list()
for (i in 1:length(colnames(cellulose_df))){
  print(i)
  trait <- as.matrix(cellulose_df)[,i]
  fit_cel[[i]] <- anc.ML(tree, trait, model="OU")
  cat("convergence for ",i,"=",fit_cel[[i]]$convergence,"\n")
  obj_cel[[i]]<-contMap(tree,trait,plot=FALSE)
  obj_cel[[i]]<-setMap(obj_cel[[i]],colors=c("#f7fcf5","#74c476","#00441b"), space="Lab")
}
names <- colnames(cellulose_df)

# create dataframe to create heatmap for ancestral states:
print("Creating heatmap for ancestral states")
rows <- c("root", "ancestral_Lecanoromycete", "ancestral_Eurotiomycete", "ancestral_Ostropomycetidae","ancestral_Lecanoromycetidae", "ancestral_Xylographa")
ancestral_states <- data.frame(matrix(nrow=length(rows), ncol=length(colnames(cellulose_df))))
rownames(ancestral_states) <- rows
colnames(ancestral_states) <- colnames(cellulose_df)
root_node <- findMRCA(tree, tree$tip.label)
print(root_node)
print("Ancestral internal nodes as follows (for manual check):")
lecanoro_node <- findMRCA(tree, c("Peltigera_leucophlebia", "Xylographa_parallela", "Umbilicaria_muehlenbergii"))
print("Lecanoromycetes:")
print(lecanoro_node)
eurotio_node <- findMRCA(tree, c("Pseudophaeomoniella_oleicola", "Penicillium_chrysogenum"))
print("Eurotiomycetes:")
print(eurotio_node)
ostropo_node <- findMRCA(tree, c("Loxospora_cismonica", "Xylographa_parallela"))
print("Ostropomycetidae:")
print(ostropo_node)
lecanidae_node <- findMRCA(tree, c("Cladonia_metacorallifera", "Peltigera_leucophlebia"))
print("Lecanoromycetidae:")
print(lecanidae_node)
xylo_node <- findMRCA(tree, c("Xylographa_parallela", "Xylographa_trunciseda"))
print("Xylographa:")
print(xylo_node)
my_nodes <- c(root_node, lecanoro_node, eurotio_node, ostropo_node, lecanidae_node, xylo_node)

for (i in 1:length(fit_cel)) {
  ancestral_states["ancestral_Xylographa", names[i]] <- round(fit_cel[[i]]$ace[toString(xylo_node)])
  ancestral_states["ancestral_Ostropomycetidae", names[i]] <- round(fit_cel[[i]]$ace[toString(ostropo_node)])
  ancestral_states["ancestral_Lecanoromycetidae", names[i]] <- round(fit_cel[[i]]$ace[toString(lecanidae_node)])
  ancestral_states["ancestral_Lecanoromycete", names[i]] <- round(fit_cel[[i]]$ace[toString(lecanoro_node)])
  ancestral_states["ancestral_Eurotiomycete", names[i]] <- round(fit_cel[[i]]$ace[toString(eurotio_node)])
  ancestral_states["root", names[i]] <- round(fit_cel[[i]]$ace[toString(root_node)])
}


print("Plotting the all heatmap...")
base_size <- 11
#reorder names for plot
ancestral_states <- ancestral_states[c("ancestral_Xylographa", "ancestral_Ostropomycetidae", "ancestral_Lecanoromycetidae", "ancestral_Lecanoromycete", "ancestral_Eurotiomycete", "root"),]
rownames(ancestral_states)
ancestral_states$name <- rownames(ancestral_states)

plot_df1 <- melt(ancestral_states)
psummary_anc <- ggplot(plot_df1, aes(variable, name)) + geom_tile(aes(fill = value)) + scale_fill_gradient2(low = "#EAF4F7",   high = "#F06449")+ theme_classic()+geom_text(aes(label = ifelse(round(plot_df1$value, 1)>0, round(plot_df1$value, 1), "")), size=2)
psummary_anc <- psummary_anc + labs(x = "", y = "") + scale_x_discrete(position="top") +scale_y_discrete(expand = c(0, 0)) + theme(plot.margin = margin(0, 1, 1, 0, "cm"),legend.position = "none",axis.ticks = element_blank(), axis.text.x = element_text(size = base_size *0.8,angle=45, hjust = 0, colour = "grey50"))
psummary_anc
plot_df_sum <- melt(cellulose_df)

colnames(plot_df_sum) <- c("label", "category", "value")
plot_df_sum <- as_tibble(plot_df_sum)
plot_df_sum$label <- as.character(plot_df_sum$label)
plot_df_sum$category <- as.character(plot_df_sum$category)

psummary <- ggtreeplot(ptree, plot_df_sum, aes(x=category)) + geom_tile(aes(fill=value)) +scale_fill_gradient2(low = "#EAF4F7", high = "#F06449") +theme_classic()+no_y_axis() +no_x_axis()+theme(axis.text.x = element_text(size = base_size *0.8, angle=45, hjust = 0, colour = "grey50")) + no_legend() +scale_x_discrete(position="top")
psummary <- psummary + scale_y_continuous(expand=c(0,0)) + geom_text(aes(label = ifelse(plot_df_sum$value >0, plot_df_sum$value, "")), colour="black", size=2)


pdf(file=paste(prefix,"_",set,"_heatmap.pdf",sep=""), width=11.7, height=8.3)
psummary
dev.off()

pdf(file=paste(prefix,"_",set,"_anc_heatmap.pdf",sep=""), width=11.7, height=8.3)
psummary_anc
dev.off()


# Prepare the combined plot:

#Make some additional changes to the anc df before combining with the other plots:
# To create the labels for the anc values while maintaining subplot alignment 
# I had to create another empty plot with just the labels, but the extend into the plotting area:
psummary_legend <- ggplot(plot_df1, aes(variable, name)) + scale_y_discrete(position = "right") + theme(axis.text.y.right = element_text(hjust = 1, size=base_size*0.6))
psummary_legend <- psummary_legend + theme(legend.position = "none", axis.title.y.right=element_blank(), panel.background = element_blank()) +no_y_axis() +no_x_axis()
psummary_legend <- psummary_legend + theme(axis.text.y.right = element_text(margin = margin(0,0.5,0,-4.5, unit = 'cm'), colour="black"))

# make same changes to the plt margin
psummary_anc <- psummary_anc + no_y_axis() + no_x_axis()+theme(plot.margin=margin(theme_grey()$plot.margin))

#remove margins around plots
#psummary_legend <- psummary_legend + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
#psummary <- psummary + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
#ptree <- ptree + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
#psummary_anc <- psummary_anc + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))


#highlight nodes in tree for plotting:
ptree <- ptree + geom_point2(aes(subset=(node %in% my_nodes)),color="black",size=2)

#plot layout:
layout <- "
AAABB
AAABB
AAABB
AAABB
AAABB
AAABB
AAABB
AAABB
AAABB
AAABB
AAABB
AAABB
AAABB
AAABB
AAABB
CCCDD
"
print("combine output")
pfiller <- ggplot() + geom_blank() + theme_minimal()
pdf(file=paste(prefix,"_",set,"_combined.pdf",sep=""), width=11.7, height=8.3)
ptree + psummary + psummary_legend + psummary_anc + plot_layout(design=layout)
dev.off()

#ptree + psummary + psummary_legend + psummary_anc + plot_layout(design=layout)
#ggarrange(ptree,psummary, psummary_legend,psummary_anc)
print("saving environment")
save.image(file=paste(prefix,"_",set,"_anc_cazy_all.RData",sep=""))
