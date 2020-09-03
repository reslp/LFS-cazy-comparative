#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

wd <- args[1]
cazy_file <- args[2]
treefile <- args[3]
prefix <- args[4]
outdir <- args[5]

print(wd)
print(cazy_file)
print(treefile)
print(prefix)
print(outdir)

library(phytools)
library(reshape2)
library(ggtree)
library(ggplot2)

setwd(wd)

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
cellulose_df <- all_cazy

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


print("Reconstructing ancestral states...")
fit_cel <- list()
obj_cel <- list()
for (i in 1:length(colnames(cellulose_df))){
  print(i)
  trait <- as.matrix(cellulose_df)[,i]
  fit_cel[[i]] <- anc.ML(tree, trait, model="OU")
  #fit_cel[[i]] <- fastAnc(tree, trait) # this is only for faster prototyping
  cat("convergence for ",i,"=",fit_cel[[i]]$convergence,"\n")
  obj_cel[[i]]<-contMap(tree,trait,plot=FALSE)
  obj_cel[[i]]<-setMap(obj_cel[[i]],colors=c("#f7fcf5","#74c476","#00441b"), space="Lab")
}
names <- colnames(cellulose_df)

print(fit_cel[[1]])

# create dataframe to create heatmap for ancestral states:
print("Ancestral internal nodes as follows (for manual check):")

print("root:")
root_node <- findMRCA(tree, tree$tip.label)
print(root_node)

print("eurotio_lecanoro_split")
eurotio_lecanoro_split_node <- findMRCA(tree, c("Endocarpon_pusillum", "Onygena_corvina", "Evernia_prunastri", "Xylographa_parallela"))
print(eurotio_lecanoro_split_node)

eurotio_node <- findMRCA(tree, c("Pseudophaeomoniella_oleicola", "Penicillium_chrysogenum"))
print("Eurotiomycetes:")
print(eurotio_node)

lecanoro_slat_node <- findMRCA(tree, c("Peltigera_leucophlebia", "Xylographa_parallela", "Umbilicaria_muehlenbergii", "Acarospora_spec"))
print("Lecanoromycetes s.lat")
print(lecanoro_slat_node)

lecanoro_node <- findMRCA(tree, c("Peltigera_leucophlebia", "Xylographa_parallela", "Umbilicaria_muehlenbergii"))
print("Lecanoromycetes:")
print(lecanoro_node)

lecanoro_ostropo_split_node <- findMRCA(tree, c("Peltigera_leucophlebia", "Xylographa_parallela", "Toensbergia_leucococca"))
print("Lecanoro ostropo split:")
print(lecanoro_ostropo_split_node)

ostropo_node <- findMRCA(tree, c("Loxospora_cismonica", "Xylographa_parallela"))
print("Ostropomycetidae:")
print(ostropo_node)

lecanidae_node <- findMRCA(tree, c("Cladonia_metacorallifera", "Peltigera_leucophlebia", "Toensbergia_leucococca"))
print("Lecanoromycetidae:")
print(lecanidae_node)

xylo_node <- findMRCA(tree, c("Xylographa_parallela", "Xylographa_trunciseda"))
print("Xylographa:")
print(xylo_node)


print("Creating dataframe for ancestral states")
rows <- c("root", "Eurotio_Lecanoro_split", "ancestral_Eurotiomycetes","ancestral_Lecanoromycetes_slat", "ancestral_Lecanoromycetes_sstr","Lecanoro_Ostropo_split", "ancestral_Ostropomycetidae","ancestral_Lecanoromycetidae", "ancestral_Xylographa")
ancestral_states <- data.frame(matrix(nrow=length(rows), ncol=length(colnames(cellulose_df))))
rownames(ancestral_states) <- rows
colnames(ancestral_states) <- colnames(cellulose_df)

my_nodes <- c(root_node, eurotio_lecanoro_split_node, eurotio_node, lecanoro_slat_node, lecanoro_node,lecanoro_ostropo_split_node, ostropo_node, lecanidae_node, xylo_node)

for (i in 1:length(fit_cel)) {
  ancestral_states["ancestral_Xylographa", names[i]] <- round(fit_cel[[i]]$ace[toString(xylo_node)])
  ancestral_states["ancestral_Ostropomycetidae", names[i]] <- round(fit_cel[[i]]$ace[toString(ostropo_node)])
  ancestral_states["ancestral_Lecanoromycetidae", names[i]] <- round(fit_cel[[i]]$ace[toString(lecanidae_node)])
  ancestral_states["ancestral_Lecanoromycetes_sstr", names[i]] <- round(fit_cel[[i]]$ace[toString(lecanoro_node)])
  ancestral_states["Lecanoro_Ostropo_split", names[i]] <- round(fit_cel[[i]]$ace[toString(lecanoro_ostropo_split_node)])
  ancestral_states["ancestral_Eurotiomycetes", names[i]] <- round(fit_cel[[i]]$ace[toString(eurotio_node)])
  ancestral_states["Eurotio_Lecanoro_split", names[i]] <- round(fit_cel[[i]]$ace[toString(eurotio_lecanoro_split_node)])
  ancestral_states["ancestral_Lecanoromycetes_slat", names[i]] <- round(fit_cel[[i]]$ace[toString(lecanoro_slat_node)])
  ancestral_states["root", names[i]] <- round(fit_cel[[i]]$ace[toString(root_node)])
}
# uncomment this for fast prototyping with fastANC
#for (i in 1:length(fit_cel)) {
#  ancestral_states["ancestral_Xylographa", names[i]] <- round(fit_cel[[i]][toString(xylo_node)])
#  ancestral_states["ancestral_Ostropomycetidae", names[i]] <- round(fit_cel[[i]][toString(ostropo_node)])
#  ancestral_states["ancestral_Lecanoromycetidae", names[i]] <- round(fit_cel[[i]][toString(lecanidae_node)])
#  ancestral_states["ancestral_Lecanoromycetes_sstr", names[i]] <- round(fit_cel[[i]][toString(lecanoro_node)])
#  ancestral_states["Lecanoro_Ostropo_split", names[i]] <- round(fit_cel[[i]][toString(lecanoro_ostropo_split_node)])
#  ancestral_states["ancestral_Eurotiomycetes", names[i]] <- round(fit_cel[[i]][toString(eurotio_node)])
#  ancestral_states["Eurotio_Lecanoro_split", names[i]] <- round(fit_cel[[i]][toString(eurotio_lecanoro_split_node)])
#  ancestral_states["ancestral_Lecanoromycetes_slat", names[i]] <- round(fit_cel[[i]][toString(lecanoro_slat_node)])
#  ancestral_states["root", names[i]] <- round(fit_cel[[i]][toString(root_node)])
#}

print(ancestral_states)

print("saving environment")
save.image(file=paste(prefix,"_anc_cazy_all.RData",sep=""))
