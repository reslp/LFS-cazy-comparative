library(ggplot2)
library(ggrepel)
library(tidyverse)
library(patchwork)
library(ggpubr)

args <- commandArgs(trailingOnly=TRUE)

rdata_to_load <- args[1]
sign_cazymes <- args[2]
genome_stats <- args[3]
outprefix <- args[4]

#load ancestral states data:
load("82_genomes_phylosig_anc_cazy_all.RData")

#load cazymes with sign K>1
sign_cazymes <- read.table("phylosig_cazymes_kgr1_pkl0.05.txt", header=T, sep=",")

#load taxonomic information for species:
add_info <- read.csv("stats_genomes.csv", sep=";", header=T) 
add_info["class"]
#all cazyme data needs to be ordered so that the taxonomy matches (add_info is ordered alphabetically)
all_cazy <- all_cazy[ order(row.names(all_cazy)), ]
#data
all_cazy
cazy_with_info <- cbind(all_cazy, add_info["class"])


#
cazy_subset <- all_cazy[,colnames(all_cazy) %in% sign_cazymes$cazyme]
#ancestral_states
#ancestral_states$name <- NULL
#cazy_subset <- rbind(cazy_subset, ancestral_states)
#rownames(cazy_subset)
#create phylogenetically informed PCA

# now create phylogenetic pca. The data will be log transformed and the results optimzied by Maximum Likelihood
a <- phyl.pca(tree, log(cazy_subset+1),opt="ML", mode="corr")
percent.var.phylo <-round(diag(a$Eval)/sum(diag(a$Eval)),4) * 100
percent.var.phylo

#pca <- prcomp(log(cazy_subset+1),center=T)
#pca$x
a$S

plot_data <- a$S
plot_data <- as.data.frame(plot_data)
plot_data$class <- cazy_with_info$class
plot_data$name <- str_replace(rownames(plot_data), "_", " ")

colors <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33")
p1 <- ggplot(plot_data, aes(PC1, PC2, colour=class)) + geom_point() +scale_color_manual(values=colors) +ggtitle("") + theme(legend.position = "none") + theme(plot.title = element_text(size = 10)) +xlab(paste("PC1 (", as.character(percent.var.phylo[1]),"%)", sep=""))+ylab(paste("PC2 (", as.character(percent.var.phylo[2]),"%)", sep=""))
p2 <- ggplot(plot_data, aes(PC1, PC2, label=rownames(plot_data))) + geom_point(color = ifelse(plot_data$class != "Lecanoromycetes", "grey50", "#984ea3")) + geom_text_repel(data = . %>% mutate(label = ifelse(class == "Lecanoromycetes" & PC2 > -1, name, "")), aes(label = label),segment.size=0.1, size=1.8) +ggtitle("")+ theme(plot.title = element_text(size = 10))+xlab(paste("PC1 (", as.character(percent.var.phylo[1]),"%)", sep=""))+ylab(paste("PC2 (", as.character(percent.var.phylo[2]),"%)", sep=""))
p3 <- ggplot(plot_data, aes(PC1, PC2, label=rownames(plot_data))) + geom_point(color = ifelse(plot_data$class != "Lecanoromycetes", "grey50", "#984ea3")) + geom_text_repel(data = . %>% mutate(label = ifelse(class == "Lecanoromycetes" & PC2 <= -1, name, "")), aes(label = label),segment.size=0.1, size=1.8) +ggtitle("")+ theme(plot.title = element_text(size = 10))+xlab(paste("PC1 (", as.character(percent.var.phylo[1]),"%)", sep=""))+ylab(paste("PC2 (", as.character(percent.var.phylo[2]),"%)", sep=""))

pp <- ggarrange(p1,p2,p3, ncol=3,nrow=1,common.legend = T, legend="bottom", labels="AUTO")
pp

pdf(file="phly_pca.pdf", width=11.7, height=5)
pp
dev.off()

# now create a plot with a normal PCA but with the reconstructed ancestral states on different nodes:
ancestral_states
ancestral_states$name <- NULL
cazy_subset2 <- rbind(cazy_subset, ancestral_states)
rownames(cazy_subset2)

pca <- prcomp(log(cazy_subset2+1),center=T)
pca_summary <-summary(pca)

#extract variance of PC1 and 2 for plot
PC1 <- as.data.frame(pca_summary$importance)$PC1[2]*100
PC2 <- as.data.frame(pca_summary$importance)$PC2[2]*100


plot_data2 <- pca$x
plot_data2 <- as.data.frame(plot_data2)
plot_data2$class <- c(cazy_with_info$class, "ancestral", "ancestral", "ancestral", "ancestral", "ancestral", "ancestral")
plot_data2$name <- str_replace(rownames(plot_data2), "_", " ")

colors2 <- c("#000000","#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33")
p1 <- ggplot(plot_data2, aes(PC1, PC2, colour=class)) + geom_point() +scale_color_manual(values=colors2) +ggtitle("") + theme(legend.position = "none") + theme(plot.title = element_text(size = 10)) +xlab(paste("PC1 (", as.character(PC1),"%)", sep=""))+ylab(paste("PC2 (", as.character(PC2),"%)", sep=""))+ geom_text_repel(data = . %>% mutate(label = ifelse(class == "ancestral", name, "")), aes(label = label),segment.size=0.1, size=1.8)
p2 <- ggplot(plot_data2, aes(PC1, PC2, label=rownames(plot_data2))) + geom_point(color = ifelse(plot_data2$class != "Lecanoromycetes", "grey50", "#984ea3")) + geom_text_repel(data = . %>% mutate(label = ifelse(class == "Lecanoromycetes" & PC2 > 0.6, name, "")), aes(label = label),segment.size=0.1, size=1.8) +ggtitle("")+ theme(plot.title = element_text(size = 10))+xlab(paste("PC1 (", as.character(PC1),"%)", sep=""))+ylab(paste("PC2 (", as.character(PC2),"%)", sep=""))
p3 <- ggplot(plot_data2, aes(PC1, PC2, label=rownames(plot_data2))) + geom_point(color = ifelse(plot_data2$class != "Lecanoromycetes", "grey50", "#984ea3")) + geom_text_repel(data = . %>% mutate(label = ifelse(class == "Lecanoromycetes" & PC2 <= 0.6, name, "")), aes(label = label),segment.size=0.1, size=1.8) +ggtitle("")+ theme(plot.title = element_text(size = 10))+xlab(paste("PC1 (", as.character(PC1),"%)", sep=""))+ylab(paste("PC2 (", as.character(PC2),"%)", sep=""))

pp <- ggarrange(p1,p2,p3, ncol=3,nrow=1,common.legend = T, legend="bottom", labels="AUTO")
pp


pdf(file="normal_pca.pdf", width=11.7, height=5)
pp
dev.off()
#<-round(diag(a$Eval)/sum(diag(a$Eval)),4) * 100
#create phylogenetically informed PCA
