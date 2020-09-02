#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

wd <- args[1]
prefix <- args[2]
outdir <- args[3]
rdata <- args[4]

print(wd)
print(prefix)
print(outdir)
print(rdata)

library(phytools)
library(ggplot2)
library(ggpubr)
library(ggtree)
library(cowplot)
library(reshape2)
library(patchwork)
library(tidyverse)

load(rdata)
setwd(wd)


# function definitions. These functions have been modified from:
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

## Run plotting routine for each interesting gene set:

sets <- c("cellulose", "hemicellulose", "phylosig")

for (foo in 1:length(sets)) {
  set <- sets[foo]
  if (set == "phylosig") {
    cazies <- read.table(paste(wd,"/results/",prefix,"/phylosig/phylosig_sign_cazymes.txt", sep=""), header=T, sep=",")
    cz <- cazies$cazyme
    cellulose_df <- all_cazy[,colnames(all_cazy) %in% cazies$cazyme]
    print("Will analyze phylosig cazy")
    print(colnames(cellulose_df))
  } else if (set == "cellulose") {
    cz <- c("GH3","GH5","GH6","GH7","GH10","GH11", "GH12", "GH28","GH43","GH45","GH61","GH74", "CBM1", "AA9", "CE1", "CE16", "CE5", "CE8", "CE12", "CE15")
    cellulose_df <- all_cazy[,colnames(all_cazy) %in% c("GH3","GH5","GH6","GH7","GH10","GH11", "GH12", "GH28","GH43","GH45","GH61","GH74", "CBM1", "AA9", "CE1", "CE16", "CE5", "CE8", "CE12", "CE15")]
    print("Will analyze the cellulose set")
  } else if (set == "hemicellulose") {
    cz <- c("GH10", "GH11", "CE1", "CE12", "CE3", "CE2", "CE5", "CE7", "GH43_29", "GH43_34", "GH51", "GH45", "GH62")
    cellulose_df <- all_cazy[,colnames(all_cazy) %in% c("GH10", "GH11", "CE1", "CE12", "CE3", "CE2", "CE5", "CE7", "GH43_29", "GH43_34", "GH51", "GH45", "GH62")]
    print("Will analyze the hemicellulose set")
  } else {
    print("No set specified")
    quit(save = "no", status = 1)
  }
  #print(cellulose_df)
  print("First plotting the heatmaps...")
  base_size <- 11
  #reorder names for plote
  ancestral_states_plot <- ancestral_states[rows,colnames(ancestral_states) %in% cz]
  print(rownames(ancestral_states_plot))  
  rows_in_order <- c("root", "Eurotio_Lecanoro_split", "ancestral_Eurotiomycetes","ancestral_Lecanoromycetes_slat", "ancestral_Lecanoromycetes_sstr","Lecanoro_Ostropo_split","ancestral_Lecanoromycetidae", "ancestral_Ostropomycetidae", "ancestral_Xylographa")
  ancestral_states_plot <- ancestral_states_plot[match(rows_in_order, rownames(ancestral_states_plot)),]
  print(rownames(ancestral_states_plot))
  ancestral_states_plot$name <- rownames(ancestral_states_plot)

  plot_df1 <- melt(ancestral_states_plot)
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


  pdf(file=paste("results/", prefix, "/plots_ancestral_states/", prefix,"_",set,"_heatmap.pdf",sep=""), width=11.7, height=8.3)
  print(psummary)
  dev.off()

  pdf(file=paste("results/", prefix, "/plots_ancestral_states/", prefix,"_",set,"_anc_heatmap.pdf",sep=""), width=11.7, height=8.3)
  print(psummary_anc)
  dev.off()

  # Prepare the combined plot:
  print("Now plot combined plot")
  # Make some additional changes to the anc df before combining with the other plots:
  # To create the labels for the anc values while maintaining subplot alignment 
  # I had to create another empty plot with just the labels, but the extend into the plotting area:
  psummary_legend <- ggplot(plot_df1, aes(variable, name)) + scale_y_discrete(position = "right") + theme(axis.text.y.right = element_text(hjust = 1, size=base_size*0.6))
  psummary_legend <- psummary_legend + theme(legend.position = "none", axis.title.y.right=element_blank(), panel.background = element_blank()) +no_y_axis() +no_x_axis()
  psummary_legend <- psummary_legend + theme(axis.text.y.right = element_text(margin = margin(0,0.5,0,-4.5, unit = 'cm'), colour="black"))

  # make same changes to the plt margin
  psummary_anc <- psummary_anc + no_y_axis() + no_x_axis() + theme(plot.margin=margin(theme_grey()$plot.margin))

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
  CCCDD
  CCCDD
"
  
  print("combine output")
  pfiller <- ggplot() + geom_blank() + theme_minimal()
  pdf(file=paste("results/", prefix, "/plots_ancestral_states/", prefix,"_",set,"_combined.pdf",sep=""), width=11.7, height=8.3)
  print(ptree + psummary + psummary_legend + psummary_anc + plot_layout(design=layout))
  dev.off()

}

print("saving environment")
save.image(file=paste("results/",prefix,"/plots_ancestral_states/ancestral_states.rData", sep=""))


