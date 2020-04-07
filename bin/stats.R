#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

data_file <- args[1]
output <- args[2]
data <- read.csv(data_file, sep=";", header = T)

#plots:
# genome size with taxonomic identity
# genome size vs. gc content
library(ggplot2)
library(scales)
library(wesanderson)
library(gridExtra)
library(ggpubr)

print("creating plots...")
colors <- wes_palette("FantasticFox1", length(unique(data$class)), type="continuous")
ggplot(data = data, aes(x = Assembly.Size, y=Percent.GC, color=class)) + geom_point()
textsize <- 8

print("genome size")
pSize <- ggplot(data=data, aes(x = isolate, y= Assembly.Size, fill=class))  +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(text = element_text(size=textsize)) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = colors) +
  xlab("")+
  ylab("Genome Size (bp)")+
  coord_flip() 

print("N50")
pN50 <- ggplot(data=data, aes(x = isolate, y= Scaffold.N50, fill=class))  +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(text = element_text(size=textsize)) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = colors) +
  xlab("")+
  ylab("N50 (bp)")+
  coord_flip() 

print("Largest")
pLargest <- ggplot(data=data, aes(x = isolate, y= Largest.Scaffold, fill=class))  +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(text = element_text(size=textsize)) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = colors) +
  xlab("")+
  ylab("Largest Scaffold (bp)")+
  coord_flip() 

print("Average")
pAverage <- ggplot(data=data, aes(x = isolate, y= Average.Scaffold, fill=class))  +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(text = element_text(size=textsize)) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = colors) +
  xlab("")+
  ylab("Average Scaffold (bp)")+
  coord_flip() 

print("NumScaf")
pNumScaf <- ggplot(data=data, aes(x = isolate, y= Num.Scaffolds, fill=class))  +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(text = element_text(size=textsize)) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = colors) +
  xlab("")+
  ylab("Largest Scaffold (bp)")+
  coord_flip()

print("GC")
pGC <- ggplot(data=data, aes(x = isolate, y= Percent.GC, fill=class))  +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(text = element_text(size=textsize)) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = colors) +
  xlab("")+
  ylab("GC Content (%)")+
  coord_flip() 

print("Creating first file")
pdf(file=paste(output, "genome_stats1.pdf", sep=""), width=11.8, height=8.3, onefile=FALSE)
p1 <- ggarrange(pSize, pN50, pLargest, pAverage, pNumScaf, pGC, nrow=2, ncol=3, common.legend = T, legend="bottom", labels="AUTO")
annotate_figure(p1, top = text_grob("Basic statistics of the studied genomes", color = "black", face = "bold", size = 14))
dev.off()

print("genes")
pGenes <- ggplot(data=data, aes(x = isolate, y= Num.Genes, fill=class))  +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(text = element_text(size=textsize)) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = colors) +
  xlab("")+
  ylab("Number of Genes")+
  coord_flip() 

print("proteins")
pProteins <- ggplot(data=data, aes(x = isolate, y= Num.Proteins, fill=class))  +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(text = element_text(size=textsize)) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = colors) +
  xlab("")+
  ylab("Number of Proteins")+
  coord_flip()

print("tRNA")
ptRNA <- ggplot(data=data, aes(x = isolate, y= Num.tRNA, fill=class))  +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(text = element_text(size=textsize)) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = colors) +
  xlab("")+
  ylab("Number of tRNAs")+
  coord_flip()

print("uniqe")
pUnique <- ggplot(data=data, aes(x = isolate, y= Unique.Proteins, fill=class))  +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(text = element_text(size=textsize)) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = colors) +
  xlab("")+
  ylab("Number of unqiue Proteins")+
  coord_flip()

print("protortho")
pProtOrth <- ggplot(data=data, aes(x = isolate, y= Prots.atleast.1.ortholog, fill=class))  +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(text = element_text(size=textsize)) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = colors) +
  xlab("")+
  ylab("No. of Proteins with at least one Ortholog")+
  coord_flip()

print("create second file")
pdf(file=paste(output, "genome_stats2.pdf", sep=""), width=11.8, height=8.3, onefile=FALSE)
p2 <- ggarrange(pGenes, pProteins, ptRNA, pUnique, pProtOrth, nrow=2, ncol=3, common.legend = T, legend="bottom", labels="AUTO")
annotate_figure(p2, top = text_grob("Statistics related to the gene content of the studied genomes", color = "black", face = "bold", size = 14))
dev.off()

print("saving environment")
save.image(file=paste(output, "statistics.RData",sep=""))
