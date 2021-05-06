library(ggplot2)
library(ggrepel)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(phytools)
args <- commandArgs(trailingOnly=TRUE)

rdata_to_load <- args[1]
#sign_cazymes <- args[2]
color_info <- args[2]
apriori_cazymes <- args[3]
genome_stats <- args[4]
outprefix <- args[5]



print("Reading input data")
#load ancestral states data:
#load("82_genomes_phylosig_anc_cazy_all.RData")
load(rdata_to_load)

#load cazymes with sign K>1

#sign_cazymes <- read.table("phylosig_cazymes_kgr1_pkl0.05.txt", header=T, sep=",")
#sign_cazymes <- read.table(sign_cazymes, header=T, sep=",")
apriori_cazymes  <- read.table(apriori_cazymes, header=T, sep=",", stringsAsFactors=F)
rownames(apriori_cazymes) <- apriori_cazymes$cazyme
#load taxonomic information for species:
#add_info <- read.csv("stats_genomes.csv", sep=";", header=T) 
add_info <- read.csv(genome_stats, sep=";", header=T) 
print("done")
print("Formatting data")
#all cazyme data needs to be ordered so that the taxonomy matches (add_info is ordered alphabetically)
all_cazy <- all_cazy[ order(row.names(all_cazy)), ]
#data
cazy_with_info <- cbind(all_cazy, add_info["class"])

setnames <- colnames(apriori_cazymes)
setnames <- setnames[2:length(setnames)]
#load and parse color info:
color_data <- read.csv(color_info,sep=",", header=T, stringsAsFactors=F)
colors <- color_data$color
names(colors) <- color_data$taxonomy
print(colors)
for (i in 1:length(setnames)) {
	set <- setnames[i]
	
	print("Subsetting data to only cazymes for:")
	print(set)
	if (set == "phylosig") {
		cazy_subset <- all_cazy[,colnames(all_cazy) %in% sign_cazymes$cazyme]
	}
	else {
		cazy_subset <- all_cazy[,colnames(all_cazy) %in% apriori_cazymes$cazyme[apriori_cazymes[,set] == 1]]
	}
	#ancestral_states
	#ancestral_states$name <- NULL
	#cazy_subset <- rbind(cazy_subset, ancestral_states)
	#rownames(cazy_subset)
	#create phylogenetically informed PCA

	# now create phylogenetic pca. The data will be log transformed and the results optimzied by Maximum Likelihood
	print("Calculating phylo pca")
	a <- phyl.pca(tree, log(cazy_subset+1),opt="ML", mode="corr")
	percent.var.phylo <-round(diag(a$Eval)/sum(diag(a$Eval)),4) * 100

	#pca <- prcomp(log(cazy_subset+1),center=T)
	#pca$x

	print("prepare plotting")
	plot_data <- a$S
	plot_data <- as.data.frame(plot_data)
	plot_data$class <- cazy_with_info$class
	plot_data$name <- str_replace(rownames(plot_data), "_", " ")

	print("creating plots")
	#colors <- c("#e41a1c", "#377eb8", "#4daf4a", "colors[1]", "#ff7f00", "#ffff33")
	print(colors)
	p1 <- ggplot(plot_data, aes(PC1, PC2, colour=class)) + geom_point() +scale_color_manual(values=colors) +ggtitle("") + theme(legend.position = "none") + theme(plot.title = element_text(size = 10)) +xlab(paste("PC1 (", as.character(percent.var.phylo[1]),"%)", sep=""))+ylab(paste("PC2 (", as.character(percent.var.phylo[2]),"%)", sep=""))
	# the cutoff is used to set which labels should be plotted in pcas A B and C according to cutoffs in PC1 (cellulose, hemicellulose) and PC2 (phylosig)
	if (set == "phylosig") {
		cutoff <- 1
		p2 <- ggplot(plot_data, aes(PC1, PC2, label=rownames(plot_data))) + geom_point(color = ifelse(plot_data$class != "Lecanoromycetes", "grey50", colors[1])) + geom_text_repel(data = . %>% mutate(label = ifelse(class == "Lecanoromycetes" & PC2 > cutoff, name, "")), aes(label = label),segment.size=0.1, size=1.8) +ggtitle("")+ theme(plot.title = element_text(size = 10))+xlab(paste("PC1 (", as.character(percent.var.phylo[1]),"%)", sep=""))+ylab(paste("PC2 (", as.character(percent.var.phylo[2]),"%)", sep=""))
	p3 <- ggplot(plot_data, aes(PC1, PC2, label=rownames(plot_data))) + geom_point(color = ifelse(plot_data$class != "Lecanoromycetes", "grey50", colors[1])) + geom_text_repel(data = . %>% mutate(label = ifelse(class == "Lecanoromycetes" & PC2 <= cutoff, name, "")), aes(label = label),segment.size=0.1, size=1.8) +ggtitle("")+ theme(plot.title = element_text(size = 10))+xlab(paste("PC1 (", as.character(percent.var.phylo[1]),"%)", sep=""))+ylab(paste("PC2 (", as.character(percent.var.phylo[2]),"%)", sep=""))
	} else if (set == "cellulose") {
		cutoff <- 3.125
                p2 <- ggplot(plot_data, aes(PC1, PC2, label=rownames(plot_data))) + geom_point(color = ifelse(plot_data$class != "Lecanoromycetes", "grey50", colors[1])) + geom_text_repel(data = . %>% mutate(label = ifelse(class == "Lecanoromycetes" & PC1 > cutoff, name, "")), aes(label = label),segment.size=0.1, size=1.8) +ggtitle("")+ theme(plot.title = element_text(size = 10))+xlab(paste("PC1 (", as.character(percent.var.phylo[1]),"%)", sep=""))+ylab(paste("PC2 (", as.character(percent.var.phylo[2]),"%)", sep=""))
        p3 <- ggplot(plot_data, aes(PC1, PC2, label=rownames(plot_data))) + geom_point(color = ifelse(plot_data$class != "Lecanoromycetes", "grey50", colors[1])) + geom_text_repel(data = . %>% mutate(label = ifelse(class == "Lecanoromycetes" & PC1 <= cutoff, name, "")), aes(label = label),segment.size=0.1, size=1.8) +ggtitle("")+ theme(plot.title = element_text(size = 10))+xlab(paste("PC1 (", as.character(percent.var.phylo[1]),"%)", sep=""))+ylab(paste("PC2 (", as.character(percent.var.phylo[2]),"%)", sep=""))
	} else if (set == "hemicellulose") {
		cutoff <- 2
		 p2 <- ggplot(plot_data, aes(PC1, PC2, label=rownames(plot_data))) + geom_point(color = ifelse(plot_data$class != "Lecanoromycetes", "grey50", colors[1])) + geom_text_repel(data = . %>% mutate(label = ifelse(class == "Lecanoromycetes" & PC1 > cutoff, name, "")), aes(label = label),segment.size=0.1, size=1.8) +ggtitle("")+ theme(plot.title = element_text(size = 10))+xlab(paste("PC1 (", as.character(percent.var.phylo[1]),"%)", sep=""))+ylab(paste("PC2 (", as.character(percent.var.phylo[2]),"%)", sep=""))
        p3 <- ggplot(plot_data, aes(PC1, PC2, label=rownames(plot_data))) + geom_point(color = ifelse(plot_data$class != "Lecanoromycetes", "grey50", colors[1])) + geom_text_repel(data = . %>% mutate(label = ifelse(class == "Lecanoromycetes" & PC1 <= cutoff, name, "")), aes(label = label),segment.size=0.1, size=1.8) +ggtitle("")+ theme(plot.title = element_text(size = 10))+xlab(paste("PC1 (", as.character(percent.var.phylo[1]),"%)", sep=""))+ylab(paste("PC2 (", as.character(percent.var.phylo[2]),"%)", sep=""))
	} else {
		cutoff <- 3.125
                p2 <- ggplot(plot_data, aes(PC1, PC2, label=rownames(plot_data))) + geom_point(color = ifelse(plot_data$class != "Lecanoromycetes", "grey50", colors[1])) + geom_text_repel(data = . %>% mutate(label = ifelse(class == "Lecanoromycetes" & PC1 > cutoff, name, "")), aes(label = label),segment.size=0.1, size=1.8) +ggtitle("")+ theme(plot.title = element_text(size = 10))+xlab(paste("PC1 (", as.character(percent.var.phylo[1]),"%)", sep=""))+ylab(paste("PC2 (", as.character(percent.var.phylo[2]),"%)", sep=""))
        p3 <- ggplot(plot_data, aes(PC1, PC2, label=rownames(plot_data))) + geom_point(color = ifelse(plot_data$class != "Lecanoromycetes", "grey50", colors[1])) + geom_text_repel(data = . %>% mutate(label = ifelse(class == "Lecanoromycetes" & PC1 <= cutoff, name, "")), aes(label = label),segment.size=0.1, size=1.8) +ggtitle("")+ theme(plot.title = element_text(size = 10))+xlab(paste("PC1 (", as.character(percent.var.phylo[1]),"%)", sep=""))+ylab(paste("PC2 (", as.character(percent.var.phylo[2]),"%)", sep=""))
	}
	pp <- ggarrange(p1,p2,p3, ncol=3,nrow=1,common.legend = T, legend="bottom", labels="AUTO")
	print("saving phylo.pca plots")
	pdf(file=paste(outprefix,"/",set,"_phly_pca.pdf", sep=""), width=11.7, height=5)
	print(pp)
	dev.off()

	print("Starting with normal PCA")
	# now create a plot with a normal PCA but with the reconstructed ancestral states on different nodes:
	ancestral_states$name <- NULL
	if (set == "phylosig") {
		ancestral_states_sub <- ancestral_states[,colnames(ancestral_states) %in% sign_cazymes$cazyme]
	} else {
		ancestral_states_sub <- ancestral_states[,colnames(ancestral_states) %in% apriori_cazymes$cazyme[apriori_cazymes[,set] == 1]]
	}
	cazy_subset2 <- rbind(cazy_subset, ancestral_states_sub)
	#rownames(cazy_subset2)

	print("Run normal PCA")
	pca <- prcomp(log(cazy_subset2+1),center=T)
	pca_summary <-summary(pca)
	print(pca_summary)
	#extract variance of PC1 and 2 for plot
	print("Extract explained variance for PC1 and PC2")
	PC1 <- as.data.frame(pca_summary$importance)$PC1[2]*100
	PC2 <- as.data.frame(pca_summary$importance)$PC2[2]*100


	print("prepare for plotting")
	plot_data2 <- pca$x
	plot_data2 <- as.data.frame(plot_data2)
	cazy_with_info$class
	plot_data2$class <- c(as.character(cazy_with_info$class), "ancestral", "ancestral", "ancestral", "ancestral", "ancestral", "ancestral", "ancestral", "ancestral", "ancestral")

	plot_data2$class
	plot_data2$name <- str_replace(rownames(plot_data2), "_", " ")

	print("Creating plots")
	ancestral_color <- "#000000"
	names(ancestral_color) <- "ancestral"
	colors2 <- c(ancestral_color, colors) #,"#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33")
	p1 <- ggplot(plot_data2, aes(PC1, PC2, colour=class)) + geom_point() +scale_color_manual(values=colors2) +ggtitle("") + theme(legend.position = "none") + theme(plot.title = element_text(size = 10)) +xlab(paste("PC1 (", as.character(PC1),"%)", sep=""))+ylab(paste("PC2 (", as.character(PC2),"%)", sep=""))+ geom_text_repel(data = . %>% mutate(label = ifelse(class == "ancestral", name, "")), aes(label = label),segment.size=0.1, size=1.8)
	p2 <- ggplot(plot_data2, aes(PC1, PC2, label=rownames(plot_data2))) + geom_point(color = ifelse(plot_data2$class != "Lecanoromycetes", "grey50", colors[1])) + geom_text_repel(data = . %>% mutate(label = ifelse(class == "Lecanoromycetes" & PC2 > 0.5, name, "")), aes(label = label),segment.size=0.1, size=1.8) +ggtitle("")+ theme(plot.title = element_text(size = 10))+xlab(paste("PC1 (", as.character(PC1),"%)", sep=""))+ylab(paste("PC2 (", as.character(PC2),"%)", sep=""))
	p3 <- ggplot(plot_data2, aes(PC1, PC2, label=rownames(plot_data2))) + geom_point(color = ifelse(plot_data2$class != "Lecanoromycetes", "grey50", colors[1])) + geom_text_repel(data = . %>% mutate(label = ifelse(class == "Lecanoromycetes" & PC2 <= 0.5, name, "")), aes(label = label),segment.size=0.1, size=1.8) +ggtitle("")+ theme(plot.title = element_text(size = 10))+xlab(paste("PC1 (", as.character(PC1),"%)", sep=""))+ylab(paste("PC2 (", as.character(PC2),"%)", sep=""))

	pp <- ggarrange(p1,p2,p3, ncol=3,nrow=1,common.legend = T, legend="bottom", labels="AUTO")

	print("Save normal PCA plots")
	pdf(file=paste(outprefix,"/",set,"_normal_pca.pdf",sep=""), width=11.7, height=5)
	print(pp)
	dev.off()
	print("all done")
	#<-round(diag(a$Eval)/sum(diag(a$Eval)),4) * 100
	#create phylogenetically informed PCA
}
