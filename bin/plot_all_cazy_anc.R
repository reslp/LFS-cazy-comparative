#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

wd <- args[1]
rfile <- args[2]
prefix <- args[3]

#wd <- "/Users/sinnafoch/Dropbox/Philipp/xylographa_comparative_genomics"
#rfile <- "tmp/82_genomes_all_anc_cazy_all.RData"

setwd(wd)


library(phytools)
library(ggplot2)
library(ggpubr)
library(ggtree)
library(cowplot)
library(reshape2)
library(patchwork)
library(tidyverse)

# load the Rdata file from the reconstruction
load(rfile)
base_size <- 7

#make sure the dataframe is divisable by 40 so that all pages look the same:
npages <- length(colnames(cellulose_df)) / 40
wanted_columns <- ceiling(npages) * 40
to_add <- seq(length(colnames(cellulose_df))+1, wanted_columns)
to_add <- paste("Z",as.character(to_add), sep="")
dd <- data.frame(cellulose_df)
dd[,to_add]<-0
# create a dummy empty row to create space between real and anc data
df <- data.frame(t(rep(0,length(colnames(dd)))))
colnames(df) <- colnames(dd)
rownames(df) <- "  "
dd <- rbind(dd, df)
cellulose_df <- as.matrix(dd)
ancestral_states$name <- NULL # column not needed for plotting
ancestral_states[,as.character(to_add)] <- 0
# colnames(ancestral_states) %in% colnames(cellulose_df) # must all be true
npages <- length(colnames(cellulose_df)) / 40 #recalculate and use value in for loop


# define order of y-axis lables so that they are the same as in the tree. This should be easier to follow
label_order <- c(tree$tip.label, "  ", rownames(ancestral_states))
cellulose_df <- rbind(cellulose_df, as.matrix(ancestral_states))

start <- 0
end <- 40

omitt_fake_labels <- function(colnms) { # function to get rid of fake column labels
  for (i in 1:length(colnms)) {
    if (startsWith(colnms[i], "Z")){
      colnms[i] <- " "
    }
  }
  return(colnms)
}

start <- 0
end <- 40

filename <- paste("results/",prefix,"/plots_ancestral_states/",prefix,"_all_cazymes.pdf", sep="")

pdf(filename, width=11.7, height=9.3)

for (i in 1:npages) {
  plot_df_sum <- melt(cellulose_df[,(start+1):end])
  colnames(plot_df_sum) <- c("label", "category", "value")
  plot_df_sum <- as_tibble(plot_df_sum)
  plot_df_sum$label <- as.character(plot_df_sum$label)
  plot_df_sum$category <- as.character(plot_df_sum$category)
  
  y_lab_text_size <- 7
  plot_df_sum$label <- factor(plot_df_sum$label, levels = rev(label_order)) # to maintain order as it is in the tree
  psummary <- ggplot(plot_df_sum, aes(x=category, y=label)) + geom_tile(aes(fill=value))  +scale_fill_gradient2(low = "#EAF4F7",   high = "#F06449") + geom_text(aes(label = ifelse(plot_df_sum$value >0, plot_df_sum$value, "")), colour="grey50", size=2)
  psummary <- psummary + theme_grey(base_size = base_size) + labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0),position="top",labels = omitt_fake_labels(unique(plot_df_sum$category))) +scale_y_discrete(expand = c(0, 0)) + theme(plot.margin = margin(0, 1, 1, 0, "cm"),legend.position = "none",axis.ticks = element_blank(), axis.text.y=element_text(size=y_lab_text_size), axis.text.x = element_text(size = base_size *0.6,angle=45, hjust = 0, colour = "grey50"))
  psummary
  
  #print(all_plots)
  print(psummary)
  
  start <- end
  print(start)
  end <- 40 * (i+1)
  print(end)
  print(" ")
}
dev.off()
save.image(file=paste("results/",prefix,"/plots_ancestral_states/ancestral_states_all.rData", sep=""))
