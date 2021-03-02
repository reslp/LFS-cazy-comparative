library(ape)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(wesanderson)
library(ggpubr)
library(patchwork)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)
print(args)


wd <- args[1]
cogs_file <- args[2]
cazy_file <- args[3]
secmet_file <- args[4]
stats_file <- args[5]
gene2gene_file <- args[6]
genelength_file <- args[7]
tree_file <- args[8]
lifestyle_file <- args[9]

#setwd("/Users/sinnafoch/Dropbox/Philipp/xylographa_comparative_genomics/tmp/new_overview")
setwd(wd)

# load all input files:
#cogs <- read.csv("COGS.all.results.csv")
cogs <- read.csv(cogs_file)
cogs <- melt(cogs)
#cazy <- read.csv("CAZyme.summary.results.csv")
cazy <- read.csv(cazy_file)
rownames(cazy) <- cazy$X
cazy$X <- NULL
cazy <- t(cazy)
#pf00083 <- read.csv("PF00083_counts.txt", sep="\t", header=F)
#rownames(pf00083) <- pf00083$V1
#pf00083$V1 <- NULL
#colnames(pf00083) <- "PF00083"
#pf00083 <- as.data.frame(pf00083)
# load transporter counts based on phylogeny:
#transporters <- read.table("sugar_transporter_counts.txt", sep="\t", header=T)
#transporters$total <- NULL
#rownames(transporters) <- transporters$group
#transporters$group <- NULL
#transporters <- t(transporters)


# this is for later to plot also individual cazymes
cazy_individual <- cazy
#gene2gene <- read.table("gene2gene_species_medians.txt",header=F)
gene2gene <- read.table(gene2gene_file, header=F)
colnames(gene2gene) <- c("name", "gene2gene")
rownames(gene2gene) <- gene2gene$name
gene2gene$name <- NULL
#genelength <- read.table("gene_length_species_medians.txt",header=F)
genelength <- read.table(genelength_file, header=F)
colnames(genelength) <- c("name", "genelength")
rownames(genelength) <- genelength$name
genelength$name <- NULL

#secmet <- read.csv("SM.summary.results.csv")
#secmet <- read.csv(secmet_file)
#secmet <- melt(secmet)
#stats <- read.table("genome.stats.summary.csv", sep=",", header=T)
stats <- read.table(stats_file, sep=",", header=T)
lifestyles <- read.table(lifestyle_file, sep=",", header=T)
rownames(lifestyles) <- lifestyles$name
tree <- read.tree(tree_file)


#get species names and set current species.
# this will have to be a for loop later
rownames(stats) <- str_replace(stats$X, " ", "_")
stats$X <- NULL
species_names <- colnames(stats)
stats2 <- t(stats)
stats2 <- as.data.frame(stats2)
stats2$isolate <- NULL
rownames(stats2) <- str_replace(rownames(stats2),"\\.","_")

#reformat to numeric for plotting
stats2$Assembly_Size <- str_replace_all(str_replace(stats2$Assembly_Size," bp",""),",","")
stats2$Assembly_Size <- trunc(as.numeric(stats2$Assembly_Size)/1000000)
stats2$Percent_GC <- str_replace(stats2$Percent_GC,"%","")
stats2$Percent_GC <- as.numeric(stats2$Percent_GC)/100
stats2$Num_Genes <- as.numeric(str_replace(stats2$Num_Genes,",",""))
stats2$Num_tRNA <- as.numeric(stats2$Num_tRNA)

# order by other dataframe, maybe later by tree tip order
stats2 <- stats2[match(tree$tip.label, rownames(stats2)),]
#stats2
lifestyles <- lifestyles[match(tree$tip.label, rownames(lifestyles)),]
lifestyles
cazy <- cazy[match(tree$tip.label, rownames(cazy)),]
#cazy
cazy_individual <- cazy_individual[match(tree$tip.label, rownames(cazy_individual)),]
#cazy_individual
gene2gene <- gene2gene[match(tree$tip.label, rownames(gene2gene)),]
genelength <- genelength[match(tree$tip.label, rownames(genelength)),]
#genelength
#pf00083 <- pf00083[match(tree$tip.label, rownames(pf00083)),]
#transporters <- transporters[match(tree$tip.label, rownames(transporters)),]
#transporters<-as.data.frame(transporters)
#calculate total number of CAZymes
cazy <- as.data.frame(cazy)
cazy_names <- rownames(cazy)
cazy <- cazy %>% mutate(Total = rowSums(.))
cazy$name <- cazy_names

#add pf00083 counts 
#stats2$PF00083 <- pf00083
#stats2
#add lifestyles
stats2$lifestyle <- lifestyles$lifestyle2
stats2$Cazy <- cazy$Total
stats2$gene2gene <- gene2gene
stats2$genelength <- genelength

#add individual cazyme values:
cazy_individual <- as.data.frame(cazy_individual)
stats2$AA <- cazy_individual$AA
stats2$CE <- cazy_individual$CE
stats2$CBM <- cazy_individual$CBM
stats2$GH <- cazy_individual$GH
stats2$GT <- cazy_individual$GT
stats2$PL <- cazy_individual$PL

#add transporter groups:
#colnames(transporters)
#for (name in colnames(transporters)) {
#  stats2[name] <- transporters[name]
#}


#add name column
stats2$name <- rownames(stats2)
stats2$name <- factor(stats2$name, levels=unique(stats2$name))

#colors
stats2$lifestyle
lifestyle_colors <- c("#4b8ab8", "#922a19", "#c3b4a5", "#6e5854", "#31394a")
names(lifestyle_colors) <- unique(stats2$lifestyle)
stats2

create_barplot <- function(df, y_data, fungal_names, title) {
  ggplot(df, aes(x=name, y=y_data, fill=lifestyle)) +
    geom_bar(stat="identity", alpha = 0.8) +
    scale_fill_manual(values = lifestyle_colors) +
    scale_y_continuous(position= "right") +
    ylab(title)+
    geom_hline(yintercept=median(y_data), linetype="dotted", color = "grey40")+
    coord_flip(ylim = NULL) + 
    theme_minimal() +
    theme(legend.position= "none",legend.title=element_blank(),axis.title.x= element_text(colour= "grey30", size=9),axis.text.x = element_text(colour= "grey30", size=9))+
    theme(axis.text.y = fungal_names ,axis.title.y=element_blank())
}

create_dotplot <- function(df, y_data, fungal_names, title) {
  ggplot(df, aes(x=name, y=y_data, color=lifestyle), label=y_data) +
    geom_point(stat="identity", alpha = 0.8, size=2) +
    scale_color_manual(values = lifestyle_colors) +
    scale_y_continuous(position= "right") +
    ylab(title)+
    geom_hline(yintercept=median(y_data), linetype="dotted", color = "grey40")+
    coord_flip(ylim = NULL) + 
    theme_minimal() +
    theme(legend.position= "none",legend.title=element_blank(),axis.title.x= element_text(colour= "grey30", size=9),axis.text.x = element_text(colour= "grey30", size=9))+
    theme(axis.text.y = fungal_names ,axis.title.y=element_blank())
}

p_genome_size <- create_barplot(stats2, stats2$Assembly_Size, element_text(colour= "gray30", size=7, vjust=0.5, hjust=1), "Genome Size (MB)")
p_GC <- create_dotplot(stats2, stats2$Percent_GC, element_blank(), "GC Content")
p_genes <- create_barplot(stats2, stats2$Num_Genes, element_blank(), "Genes")
p_trna <- create_barplot(stats2, stats2$Num_tRNA, element_blank(), "tRNAs")
p_cazy <- create_dotplot(stats2, stats2$Cazy, element_blank(), "CAZymes")
p_gene2gene <- create_barplot(stats2, stats2$gene2gene, element_blank(), "gene2gene (median)")
p_genelength <- create_barplot(stats2, stats2$genelength, element_blank(), "gene length (median)")

layout <- 
"
  ABCDEFG
  ABCDEFG
"
pdf(file="genomes.pdf", width=11, height=11.7)
print(p_genome_size + p_GC + p_genes + p_genelength + p_gene2gene + p_trna + p_cazy + plot_layout(design=layout))
dev.off() 

ggarrange(p_genome_size, p_GC,p_genes,p_trna,p_cazy, ncol=5)


# plot individual cazyme groups the same way:



p_AA <- create_barplot(stats2, stats2$AA, element_text(colour= "gray30", size=7, vjust=0.5, hjust=1), "AA")
p_CBM <- create_barplot(stats2, stats2$CBM, element_blank(), "CBM")
p_CE <- create_barplot(stats2, stats2$CE, element_blank(), "CE")
p_GH <- create_barplot(stats2, stats2$GH, element_blank(), "GH")
p_GT <- create_barplot(stats2, stats2$GT, element_blank(), "GT")
p_PL <- create_barplot(stats2, stats2$PL, element_blank(), "PL")
#p_pf00083 <- create_barplot(stats2, stats2$PF00083, element_blank(), "PF00083")

layout2 <- 
  "
  ABCDEF
  ABCDEF
"
pdf(file="cazymes_overview_plus_transporters.pdf", width=11, height=11.7)
print(p_AA + p_CBM + p_CE + p_GH + p_GT + p_PL + plot_layout(design=layout2))
dev.off()

write.csv(file="median_values_lifestyles.csv", aggregate(.~lifestyle, data=stats2, median), row.names = F, quote=F)

# now plot different transporter groups the same way:
#colnames(transporters)

#p_pent_gluc <- create_barplot(stats2, stats2$group_pentose_glucose, element_text(colour= "gray30", size=7, vjust=0.5, hjust=1), "pentose_glucose")
#p_inositol <- create_barplot(stats2, stats2$group_inositol, element_text(colour= "gray30", size=7, vjust=0.5, hjust=1), "inositol")
#p_poly_fru_ino <- create_barplot(stats2, stats2$group_polyol_fructose_inositol, element_blank(), "poly_fru_ino")
#p_cellodex <- create_barplot(stats2, stats2$group_cellodextrin_lactose, element_blank(), "cellodextrin_lactose")
#p_malt <- create_barplot(stats2, stats2$group_maltose_sucrose, element_blank(), "maltose_sucrose")
#p_gluc_pent <- create_barplot(stats2, stats2$group_glucose_pentose, element_blank(), "glucose_pentose")
#p_put_vak <- create_barplot(stats2, stats2$group_putative_vakuolar, element_blank(), "put_vakuolar")
#p_hex <- create_barplot(stats2, stats2$group_hexose, element_blank(), "hexose")
#p_sugaralc <- create_barplot(stats2, stats2$group_additional_sugaralcohol, element_blank(), "add_sugaralc")
#p_xylose <- create_barplot(stats2, stats2$group_xylose, element_blank(), "xylose")
#p_dgala <- create_barplot(stats2, stats2$group_d_galacturonic, element_blank(), "d_galactu")
#p_quest_no_leca <- create_barplot(stats2, stats2$group_question_no_lecanoros, element_blank(), "unknown1")
#p_unknown2 <- create_barplot(stats2, stats2$group_unknown, element_blank(), "unknown2")

#layout2 <- 
#  "
#  ABCDE
#  ABCDE
#"
#pdf(file="overview_sugar_transporters1.pdf", width=11, height=11.7)
#print(p_inositol + p_poly_fru_ino + p_sugaralc + p_cellodex + p_xylose + plot_layout(design=layout2))
#dev.off()

#layout2 <- 
#  "
#  ABCDEFGH
#  ABCDEFGH
#"
#p_gluc_pent
#pdf(file="overview_sugar_transporters2.pdf", width=11, height=11.7)
#print(p_pent_gluc + p_malt + p_gluc_pent + p_put_vak + p_hex + p_dgala + p_quest_no_leca + p_unknown2 + plot_layout(design=layout2))
#dev.off()
