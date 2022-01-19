#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

wd <- args[1]
cazy_file <- args[2]
treefile <- args[3]
prefix <- args[4]

print("Setting working directory to:")
print(wd)
setwd(wd)
#setwd('/Volumes/sinnafoch/sinnafoch/Dropbox/Philipp/Genomes/02_cazy_distribution/anc')

suppressWarnings(suppressMessages(library(phytools, quietly = TRUE)))
suppressWarnings(suppressMessages(library(ggplot2, quietly = TRUE)))
suppressWarnings(suppressMessages(library(ggpubr, quietly = TRUE)))
suppressWarnings(suppressMessages(library(ggtree, quietly = TRUE)))
suppressWarnings(suppressMessages(library(cowplot, quietly = TRUE)))
suppressWarnings(suppressMessages(library(reshape2, quietly = TRUE)))


# reformat the tree to ages between 0 and 1:
# Loading tree:

print("Loading tree...")
print(treefile)
tree <- read.tree(treefile)
mbt <- max(branching.times(tree))
tree$edge.length <- tree$edge.length / mbt
ptree <- ggtree(tree) +geom_tiplab() +xlim(NA, 2.2) +geom_treescale(x=0,y=45, width=0.5)
#ptree
d <- fortify(tree)
d <- subset(d, isTip)
tip_order <- with(d, label[order(y, decreasing=T)])




print("Read cazy summary input file...")
inputfile = paste(wd, "/", cazy_file,sep="")
print(inputfile)
summary_cazy <- read.csv(file=cazy_file, header=T)
rownames(summary_cazy) <- summary_cazy$X
summary_cazy$X <- NULL
summary_cazy <- t(summary_cazy)

# order dataframe according to tip labels in tree
summary_cazy <- summary_cazy[match(tip_order, rownames(summary_cazy)),]
summary_cazy <- summary_cazy [nrow(summary_cazy ):1,]
#summary_cazy
base_size <- 12

#change working directory for output
#print("Change WD for output...")
#setwd(paste(wd,"/results/",prefix,sep=""))

print("plotting tree...")
tree_out <- paste(prefix, "_tree.pdf", sep="")
pdf(file=tree_out, height=11.3, width=7.8) 
plot(tree, edge.width=2, cex=1.2, adj=0)
nodelabels()
add.scale.bar()
ptree
dev.off()

plot_df_sum <- melt(summary_cazy)
psummary <- ggplot(plot_df_sum, aes(Var2, Var1)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient(low = "#EAF4F7",   high = "#F06449")+ geom_text(aes(fill = plot_df_sum$value),label = round(plot_df_sum$value, 1), size=2)
psummary <- psummary+ theme_grey(base_size = base_size) + labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0),position="top") +scale_y_discrete(expand = c(0, 0)) + theme(plot.margin = margin(0, 1, 1, 0, "cm"),legend.position = "none",axis.ticks = element_blank(), axis.text.x = element_text(size = base_size *0.8,angle=45, hjust = 0, colour = "grey50"))
#psummary
pdf(file="overview_heatmap.pdf", width=4, height=10)
psummary
dev.off()


#ancestral state reconstruction of overview of cazymes
print("Reconstructing ancestral states")
fit_sum <- list()
obj_sum <- list()
for (i in 1:6){
  trait <- as.matrix(summary_cazy)[,i]
  fit_sum[[i]] <- anc.ML(tree, trait, model="OU")
  cat("convergence for ",i,"=",fit_sum[[i]]$convergence,"\n")
  obj_sum[[i]]<-contMap(tree,trait,plot=FALSE)
  obj_sum[[i]]<-setMap(obj_sum[[i]],colors=c("#ffffcc","#fc4e2a","#800026"), space="Lab")
}

names <- colnames(summary_cazy)
names

print("Reformatting ancestral state results")
rows <- c("root", "ancestral_Lecanoromycete", "ancestral_Eurotiomycete", "ancestral_Xylographa")
ancestral_states_summary <- data.frame(matrix(nrow=length(rows), ncol=length(colnames(summary_cazy))))
rownames(ancestral_states_summary) <- rows
colnames(ancestral_states_summary) <- colnames(summary_cazy)
for (i in 1:length(fit_sum)) {
  ancestral_states_summary["root", names[i]] <- round(fit_sum[[i]]$ace["50"])
  ancestral_states_summary["ancestral_Lecanoromycete", names[i]] <- round(fit_sum[[i]]$ace["51"])
  ancestral_states_summary["ancestral_Eurotiomycete", names[i]] <- round(fit_sum[[i]]$ace["80"])
  ancestral_states_summary["ancestral_Xylographa", names[i]] <- round(fit_sum[[i]]$ace["59"])
}

print("saving environment")
save.image(file="anc_cazy_summary.RData")

print("create anc heatmap")
base_size <- 11
ancestral_states_summary$name <- rownames(ancestral_states_summary)
plot_df1 <- melt(ancestral_states_summary)
psummary_anc <- ggplot(plot_df1, aes(variable, name)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient(low = "#EAF4F7",   high = "#F06449")+ geom_text(aes(fill = plot_df1$value),label = round(plot_df1$value, 1), size=2)
psummary_anc <- psummary_anc  + theme_grey(base_size = base_size) + labs(x = "", y = "") + scale_x_discrete(expand = c(0, 0),position="top") +scale_y_discrete(expand = c(0, 0)) + theme(plot.margin = margin(0, 1, 1, 0, "cm"), legend.position = "none", axis.ticks = element_blank(), axis.text.x = element_text(size = base_size *0.8,angle=45, hjust = 0, colour = "grey50"))
pdf(file="summary_anc_heatmap.pdf", width=11.4, height=8.3)
psummary_anc 
dev.off()


print("combine all output")
pfiller <- ggplot()+theme_minimal()
ggarrange(ptree, psummary, pfiller, psummary_anc)
library(gridExtra)
pdf(file="all_out.pdf", width=22.6, height=15.6)
#grid.arrange(ptree, psummary,pfiller, psummary_anc, ncol=2)  
plot_grid(ptree, psummary,NULL, psummary_anc,ncol=2, labels=c("A", "B", "", "C"))
dev.off()


