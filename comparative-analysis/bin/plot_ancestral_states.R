#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

wd <- args[1]
prefix <- args[2]
outdir <- args[3]
rdata <- args[4]
envfile <- args[5]

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

sessionInfo()
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
    #print("here")
    mapping <- modifyList(aes_(y=~y), mapping)
    data <- mutate(data, y=tree_y(ggtree, data))
    gg <- ggplot(data=data, mapping = mapping, ...) +
      scale_y_continuous(limits=limits, expand=c(0,0))
  }
  return(gg)
}

no_legend <- function() theme(legend.position="none")

scale_y_tree <- function(expand=expand_scale(0, 0.6), ...){
  scale_y_continuous(expand=expand, ...)
}
## function to create summary plots for each genset for fig1:
get_cazy_summary <- function(df, rows_in_order) {
  print("getting summary")
  types <- c("AA", "GH", "CE", "PL", "GT", "CBM")
  #print(typeof(df))
  for (i in 1:length(types)) {
  colname <- paste(types[i], "_sum", sep="")
  #print(colname)
  trues <- grepl(types[i],colnames(df))
  trues <- length(trues[trues == TRUE])
  #print(trues) 
  if (trues == 0) {
    dd <- rep(0, length(rownames(df)))
    names(dd) <- rownames(df)
    dd <- as.data.frame(dd)
    colnames(dd) <- colname
    df <- cbind(df,dd)
    print("will skip")
    next
  }
  if (trues ==1) {
    dd <- df[,grepl(types[i],colnames(df))]
    names(dd) <- rownames(df)
    dd <- as.data.frame(dd)
    colnames(dd)<-colname
    df <- cbind(df,dd)
    next
  }
  if (trues > 1) {
    dd <- rowSums(df[,grepl(types[i],colnames(df))])
    names(dd) <- rownames(df)
    dd <- as.data.frame(dd)
    colnames(dd) <- colname
    df <- cbind(df,dd)
    next
  }
  }
  #print(df)
  summary_df <- df[,grepl("_sum", colnames(df))]
  colnames(summary_df) <- sub("_.*", "", colnames(summary_df))
  summary_df <- summary_df[types]

  # now do some reformating to make the df ready for plotting.	
  #summary_df$name <- rownames(summary_df)
  #summary_df <- rename_anc_labels(summary_df, rows_in_order)
  #summary_df <- melt(summary_df)
  #print(head(summary_df)) 
  #colnames(summary_df) <- c("label", "category", "value")
  #summary_df <- as_tibble(summary_df)
  #summary_df$label <- factor(summary_df$label, levels=c("root (R)", "Eurotiomycetes/Lecanoromycetes split (ELS)", "ancestral Eurotiomycetes (AE)", "ancestral Lecanoromycetes s.lat. (ALSL)", "ancestral Lecanoromycetes s.str. (ALSS)", "Lecanoromycetidae/Ostropomycetidae split (OLS)", "ancestral Lecanoromycetidae (AL)", "ancestral Ostropomycetidae (AO)", "ancestral Xylographa (AX)")) 
  #summary_df$category <- as.character(summary_df$category)
  #print("combine done")
  #print(summary_df)
  return(summary_df)
}


# function to rename ancestral state names for fig1 and others:
rename_anc_labels <- function(df, rows_in_order) {
  print("Renaming anc labels...")
  df$name <- factor(df$name, levels=rows_in_order)
  levels(df$name)[levels(df$name)=="root"] ="root (R)"
  levels(df$name)[levels(df$name)=="Eurotio_Lecanoro_split"] ="Eurotiomycetes/Lecanoromycetes split (ELS)"
  levels(df$name)[levels(df$name)=="ancestral_Eurotiomycetes"] ="ancestral Eurotiomycetes (AE)"
  levels(df$name)[levels(df$name)=="ancestral_Lecanoromycetes_slat"] ="ancestral Lecanoromycetes s.lat. (ALSL)"
  levels(df$name)[levels(df$name)=="ancestral_Lecanoromycetes_sstr"] ="ancestral Lecanoromycetes s.str. (ALSS)"
  levels(df$name)[levels(df$name)=="Lecanoro_Ostropo_split"] ="Lecanoromycetidae/Ostropomycetidae split (OLS)"
  levels(df$name)[levels(df$name)=="ancestral_Lecanoromycetidae"] ="ancestral Lecanoromycetidae (AL)"
  levels(df$name)[levels(df$name)=="ancestral_Ostropomycetidae"] ="ancestral Ostropomycetidae (AO)"
  levels(df$name)[levels(df$name)=="ancestral_Xylographa"] ="ancestral Xylographa (AX)"
  return(df)
}

create_anc_plot <- function(df){
  
  print("Creating anc plot")
  #print(df)
  p <- ggplot(df, aes(y=label, x=category)) + geom_tile(aes(fill = value)) + scale_fill_gradient2(low = "#EAF4F7",   high = "#F06449")+ theme_classic()+geom_text(aes(label = ifelse(round(df$value, 1)>0, round(df$value, 1), "")), size=2)
  p <- p + labs(x = "", y = "") + scale_x_discrete(position="top") +scale_y_discrete(expand = c(0, 0)) + theme(plot.margin = margin(0, 1, 1, 0, "cm"),legend.position = "none",axis.ticks = element_blank(), axis.text.x = element_text(size = base_size *0.8,angle=45, hjust = 0, colour = "grey50"))
  # make same changes to the plt margin
  p <- p + no_y_axis() + no_x_axis() + theme(plot.margin=margin(theme_grey()$plot.margin))

  return(p)

}

create_heatmap <- function(df) {
  print("create heatmap plot")
  #print(head(df))
  ptree <- get_tree()
  p <- ggtreeplot(ptree, df, aes(x=category)) + geom_tile(aes(fill=value)) +scale_fill_gradient2(low = "#EAF4F7", high = "#F06449") +theme_classic()+no_y_axis() +no_x_axis()+theme(axis.text.x = element_text(size = base_size *0.8, angle=45, hjust = 0, colour = "grey50")) + no_legend() +scale_x_discrete(position="top")
  p <- p + scale_y_continuous(expand=c(0,0)) + geom_text(aes(label = ifelse(df$value > 0, df$value, "")), colour="black", size=2)
  return(p)
}

## correct tip labels in tree:
tree$tip.label <- gsub("_", " ", tree$tip.label)
rownames(all_cazy) <- gsub("_", " ", rownames(all_cazy))

get_tree <- function() {
  ptree <- ggtree(tree, size=0.3) + geom_tiplab(align=TRUE, size=2.5) + geom_treescale(x=0.1,y=50, width=0.2, linesize=0.3) + scale_x_continuous(expand=expand_scale(0.2)) + scale_y_tree() +xlim(NA, 1.3)

  #highlight nodes in tree for plotting:
  bla <- c("R", "ELS", "ALSL", "ALSS", "OLS", "AO", "AX", "AL", "AE")
  #print(bla)
  #print(my_nodes)
  ptree <- ptree + geom_point2(aes(subset=(node %in% my_nodes)),color="black",size=2) + geom_text2(label=bla, aes(subset=(node %in% my_nodes)), nudge_x=-0.03,nudge_y=0.7, size=2)
  return(ptree)
}

# to expand empty dataframe
cbind.all <- function (...) 
{
    nm <- list(...)
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow))
    do.call(cbind, lapply(nm, function(x) rbind(x, matrix(, n - nrow(x), ncol(x)))))
}

## Run plotting routine for each interesting gene set:

sets <- c("cellulose", "hemicellulose", "phylosig", "cell_hemi", "lignin", "pectin")
sum_plots <- list()
anc_sum_plots <- list()
single_number_summary <- list()
single_number_anc_summary <- list()
for (foo in 1:length(sets)) {
  set <- sets[foo]
  if (set == "phylosig") {
    	cazies <- read.table(paste(wd,"/results/",prefix,"/phylosig/phylosig_sign_cazymes.txt", sep=""), header=T, sep=",")
    	cz <- cazies$cazyme
    	cellulose_df <- all_cazy[,colnames(all_cazy) %in% cazies$cazyme]
    	print("Will analyze phylosig cazy")
    	#print(colnames(cellulose_df))
  } else if (set == "cellulose") {
	cz <- c("GH3","GH5","GH6","GH7","GH10","GH11", "GH12", "GH28","GH43","GH45","GH61","GH74", "CBM1", "AA9", "CE1", "CE16", "CE5", "CE8", "CE12", "CE15")
	cellulose_df <- all_cazy[,colnames(all_cazy) %in% c("GH3","GH5","GH6","GH7","GH10","GH11", "GH12", "GH28","GH43","GH45","GH61","GH74", "CBM1", "AA9", "CE1", "CE16", "CE5", "CE8", "CE12", "CE15")]
    	print("Will analyze the cellulose set")
  } else if (set == "hemicellulose") {
	cz <- c("GH10", "GH11", "CE1", "CE12", "CE3", "CE2", "CE5", "CE7", "GH43_29", "GH43_34", "GH51", "GH45", "GH62")
	cellulose_df <- all_cazy[,colnames(all_cazy) %in% c("GH10", "GH11", "CE1", "CE12", "CE3", "CE2", "CE5", "CE7", "GH43_29", "GH43_34", "GH51", "GH45", "GH62")]
	print("Will analyze the hemicellulose set")
  } else if (set == "cell_hemi") {
	cz <- c("AA3","AA9","GH61","CBM1", "CBM35", "GH26", "GH1", "GH12", "GH45", "GH6", "GH7", "GH141", "GH29", "GH31", "GH5", "GH51", "GH55", "GH74", # this line are all cellulose
		"GH10","GH11","GH16","GH2", "GH26", "GH27", "CBM13", "CBM35", "GH3", "GH30", "GH35", "GH36", "GH39", "GH43", "CBM35", "CBM6", "CBM66", "GH67", "GH72", "GH95")
	cellulose_df <- all_cazy[,colnames(all_cazy) %in% cz]
	print("Will analyze combined cellulose + hemicellulose according to Miyauchi et al. 2020")
   } else if (set == "lignin") {
	cz <- c("AA1", "AA2", "AA5")
	cellulose_df <- all_cazy[,colnames(all_cazy) %in% cz]
	print("Will analyze lignin according to Miyauchi et al. 2020")
   } else if (set == "pectin") {
	cz <- c("CE8", "GH105", "GH28", "GH49", "GH53", "GH79", "GH88", "PL1", "PL3", "PL4", "PL9")
	cellulose_df <- all_cazy[,colnames(all_cazy) %in% cz]
	print("Will analyze pectin set according Miyauchi et al. 2020")	
   }
    else {
    print("No set specified")
    quit(save = "no", status = 1)
  }

  print("Finally create single number summary. This will give a single gene count for each set")
  single_number_summary[[set]] <- rowSums(cellulose_df) 

  #print(cellulose_df)
  print("First create heatmap plots...")
  base_size <- 11
  #reorder names for plot
  ancestral_states_plot <- ancestral_states[rows,colnames(ancestral_states) %in% cz]
  rows_in_order <- c("root", "Eurotio_Lecanoro_split", "ancestral_Eurotiomycetes","ancestral_Lecanoromycetes_slat", "ancestral_Lecanoromycetes_sstr","Lecanoro_Ostropo_split","ancestral_Lecanoromycetidae", "ancestral_Ostropomycetidae", "ancestral_Xylographa")
  ancestral_states_plot <- ancestral_states_plot[match(rows_in_order, rownames(ancestral_states_plot)),]
  print("Save single number for anc data")
  single_number_anc_summary[[set]] <- rowSums(ancestral_states_plot)
  ancestral_states_plot$name <- rownames(ancestral_states_plot)



  # this is where the summary plot for ancestral will be created:
  print("Now get the summary for ancestral states based on CAZyme categories")
  anc_sum <- get_cazy_summary(ancestral_states_plot, rows_in_order)
  anc_sum$name <- rownames(anc_sum)
  anc_sum <- rename_anc_labels(anc_sum, rows_in_order)
  anc_sum <- melt(anc_sum)
  colnames(anc_sum) <- c("label", "category", "value")
  anc_sum <- as_tibble(anc_sum)
  anc_sum$label <- factor(anc_sum$label, levels=c("root (R)", "Eurotiomycetes/Lecanoromycetes split (ELS)", "ancestral Eurotiomycetes (AE)", "ancestral Lecanoromycetes s.lat. (ALSL)", "ancestral Lecanoromycetes s.str. (ALSS)", "Lecanoromycetidae/Ostropomycetidae split (OLS)", "ancestral Lecanoromycetidae (AL)", "ancestral Ostropomycetidae (AO)", "ancestral Xylographa (AX)")) 
  anc_sum$category <- as.character(anc_sum$category)
  print("Create the summarized ancestral states heatmap plot") 
  anc_sum_plots[[foo]] <- create_anc_plot(anc_sum)
  
  # this is where the extended anc plot is created
  print("Create the complete ancestral states plot based on individual families...")
  plot_df1 <- ancestral_states_plot
  plot_df1$name <- rownames(ancestral_states_plot)
  plot_df1 <- rename_anc_labels(plot_df1, rows_in_order)
  plot_df1 <- melt(plot_df1)
  colnames(plot_df1) <- c("label", "category", "value")
  plot_df1 <- as_tibble(plot_df1)
  plot_df1$label <- factor(plot_df1$label, levels=c("root (R)", "Eurotiomycetes/Lecanoromycetes split (ELS)", "ancestral Eurotiomycetes (AE)", "ancestral Lecanoromycetes s.lat. (ALSL)", "ancestral Lecanoromycetes s.str. (ALSS)", "Lecanoromycetidae/Ostropomycetidae split (OLS)", "ancestral Lecanoromycetidae (AL)", "ancestral Ostropomycetidae (AO)", "ancestral Xylographa (AX)")) 
  plot_df1$category <- as.character(plot_df1$category)
  
  psummary_anc <- create_anc_plot(plot_df1) 
  
  # next we create the heatmap of cazy counts summary:
  print("Now calculate the summary for extant CAZyme numbers based on CAZyme categories...")
  sum <- get_cazy_summary(cellulose_df)
  sum$name <- rownames(sum)
  sum <- melt(sum)
  colnames(sum) <- c("label", "category", "value")
  sum <- as_tibble(sum)
  sum$label <- as.character(sum$label)
  sum$category <- as.character(sum$category)
  #print(head(sum))
  print("create plot for extant CAZy summary")
  sum_plots[[foo]] <- create_heatmap(sum)
  
  #now the heatmap for the cazies found in the genomes:
  print("Create complete heatmap...")
  plot_df_sum <- melt(cellulose_df)

  colnames(plot_df_sum) <- c("label", "category", "value")
  plot_df_sum <- as_tibble(plot_df_sum)
  plot_df_sum$label <- as.character(plot_df_sum$label)
  plot_df_sum$category <- as.character(plot_df_sum$category)
  #print(plot_df_sum$category)
  psummary <- create_heatmap(plot_df_sum)
  print("done")
 
 
  
  #pdf(file=paste("results/ancestral_states_cazy/plots/", prefix,"_",set,"_heatmap.pdf",sep=""), width=11.7, height=8.3)
  #print(psummary)
  #dev.off()

  #pdf(file=paste("results/ancestral_states_cazy/plots/", prefix,"_",set,"_anc_heatmap.pdf",sep=""), width=11.7, height=8.3)
  #print(psummary_anc)
  #dev.off()

  # Prepare the combined plot:
  print("Now plot combined plot")
  # Make some additional changes to the anc df before combining with the other plots:
  # To create the labels for the anc values while maintaining subplot alignment 
  # I had to create another empty plot with just the labels, but the extend into the plotting area:
  psummary_legend <- ggplot(plot_df1, aes(y=label, x=category)) + scale_y_discrete(position = "right") + theme(axis.text.y.right = element_text(hjust = 1, size=base_size*0.6))
  psummary_legend <- psummary_legend + theme(legend.position = "none", axis.title.y.right=element_blank(), panel.background = element_blank()) +no_y_axis() +no_x_axis()
  psummary_legend <- psummary_legend + theme(axis.text.y.right = element_text(margin = margin(0,0.5,0,-4.5, unit = 'cm'), colour="black"))
  
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
  AAABB
  CCCDD
  CCCDD
"
  
  print("combine output")
  pfiller <- ggplot() + geom_blank() + theme_minimal()
  pdf(file=paste("results/ancestral_states_cazy/plots/", prefix,"_",set,"_combined.pdf",sep=""), width=11.7, height=8.3)
  print(get_tree() + psummary + psummary_legend + psummary_anc + plot_layout(design=layout))
  dev.off()

}

# now plot the summarized cazymes for fig1:
layout <- "
  AAAAAAAAABBCCDD
  AAAAAAAAABBCCDD
  AAAAAAAAABBCCDD
  AAAAAAAAABBCCDD
  AAAAAAAAABBCCDD
  AAAAAAAAABBCCDD
  AAAAAAAAABBCCDD
  AAAAAAAAABBCCDD
  AAAAAAAAABBCCDD
  AAAAAAAAABBCCDD
  AAAAAAAAABBCCDD
  AAAAAAAAABBCCDD
  AAAAAAAAABBCCDD
  AAAAAAAAABBCCDD
  AAAAAAAAABBCCDD
  AAAAAAAAABBCCDD
  EEEEEEEEEFFGGHH
  EEEEEEEEEFFGGHH"

print("combine summary")
pdf(file=paste("results/ancestral_states_cazy/plots/", prefix,"_",set,"_combined_summary.pdf",sep=""), width=11.7, height=8.3)
print(get_tree() + sum_plots[[1]] +sum_plots[[2]] + sum_plots[[3]] + psummary_legend + anc_sum_plots[[1]] + anc_sum_plots[[2]] +anc_sum_plots[[3]] + plot_layout(design=layout))
dev.off()

# create single number summary plots for figure 1 plot.
plot_df <- as.data.frame(single_number_summary)
plot_df <- plot_df[,c("cell_hemi", "lignin", "pectin")]
plot_df$names <- rownames(plot_df)
print(head(plot_df))
plot_df <- melt(plot_df)
print(head(plot_df))
colnames(plot_df) <- c("label", "category", "value")
plot_df <- as_tibble(plot_df)
plot_df$label <- as.character(plot_df$label)
plot_df$category <- as.character(plot_df$category)

single_number_summary_plot <- create_heatmap(plot_df)

plot_df1 <- as.data.frame(single_number_anc_summary)
plot_df1 <- plot_df1[,c("cell_hemi", "lignin", "pectin")]
plot_df1$name <- rownames(plot_df1)
print(plot_df1)
print(rows_in_order)
plot_df1 <- rename_anc_labels(plot_df1, rows_in_order)
plot_df1 <- melt(plot_df1)
print(plot_df1)
colnames(plot_df1) <- c("label", "category", "value")
plot_df1 <- as_tibble(plot_df1)
plot_df1$label <- factor(plot_df1$label, levels=c("root (R)", "Eurotiomycetes/Lecanoromycetes split (ELS)", "ancestral Eurotiomycetes (AE)", "ancestral Lecanoromycetes s.lat. (ALSL)", "ancestral Lecanoromycetes s.str. (ALSS)", "Lecanoromycetidae/Ostropomycetidae split (OLS)", "ancestral Lecanoromycetidae (AL)", "ancestral Ostropomycetidae (AO)", "ancestral Xylographa (AX)"))
plot_df1$category <- as.character(plot_df1$category)

single_number_anc_summary_plot <- create_anc_plot(plot_df1)
layout <- "
  AAAAAAAB
  AAAAAAAB
  AAAAAAAB
  AAAAAAAB
  AAAAAAAB
  AAAAAAAB
  AAAAAAAB
  AAAAAAAB
  AAAAAAAB
  AAAAAAAB
  AAAAAAAB
  AAAAAAAB
  AAAAAAAB
  AAAAAAAB
  AAAAAAAB
  AAAAAAAB
  AAAAAAAB
  CCCCCCCD
  CCCCCCCD
"

treep <- get_tree()
treep <- treep +xlim(NA,1.2)
pdf(file=paste("results/ancestral_states_cazy/plots/", prefix,"_",set,"_single_number_summary.pdf",sep=""), width=11.7, height=8.3)
print(treep + single_number_summary_plot + psummary_legend + single_number_anc_summary_plot + plot_layout(design=layout))
dev.off()

df <- as.data.frame(single_number_anc_summary)
print(df)

print("saving environment")
save.image(envfile)


