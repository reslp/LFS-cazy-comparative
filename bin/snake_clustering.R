#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

wd <- args[1]
file_cazy <- args[2]
file_interpro <- args[3]
file_pfam <- args[4]
prefix <- args[5]

setwd(wd)
library(pvclust)
print("Read input files:")
all_cazy <- read.csv(file=file_cazy)
all_pfam <- read.csv(file=file_pfam)
all_interpro <- read.csv(file=file_interpro)

boots = 5000
variance = 0.1
widthp=15
heightp=10
method="mcquitty"

print("Set WD to output")

setwd(paste(wd,"/results/",prefix,"/similarity_clustering",sep=""))
rownames(all_cazy) <- all_cazy$X
all_cazy$X <- NULL

all_cazy <- t(all_cazy)
rownames(all_cazy)

#remove columns with zero variance, as they will not help with the clustering and pvclust does not like them
print("calculating Cazyme similarity..." )
all_cazy <- all_cazy + 1
all_cazy <- log(all_cazy)

variances<-apply(t(all_cazy), 1, var)
print(variances)
low_var_names <- names(variances[which(variances<=variance)])
all_cazy_subset <- all_cazy[, !colnames(all_cazy) %in% low_var_names]
length(colnames(all_cazy))
length(colnames(all_cazy_subset))

cazy.pv <- pvclust(t(all_cazy_subset),method.hclust=method,method.dist="correlation", nboot=boots,parallel=T)
pdf("cazy_clustering_all.pdf", width=widthp, height=heightp)
plot(cazy.pv, cex=1, cex.pv=0.9, main=paste("CAZY profile similarity clustering. method: ", method," variance: ", variance, sep=""))
dev.off()


######## PFAM
#print("calculating PFAM similarity..." )
#print(rownames(all_pfam))
#rownames(all_pfam) <- all_pfam$X
#all_pfam$X <- NULL
#all_pfam$description <- NULL
#all_pfam <- all_pfam[,-ncol(all_pfam)]
#all_pfam <- t(all_pfam)

#as.numeric(all_pfam)


#remove columns with zero variance, as they will not help with the clustering and pvclust does not like them
#variances<-apply(t(all_pfam), 1, var)
#low_var_names <- names(variances[which(variances<=variance)])
#all_pfam_subset <- all_pfam[, !colnames(all_pfam) %in% low_var_names]

#as.numeric(all_pfam_subset)

#pfam.pv <- pvclust(t(all_pfam_subset),method.hclust=method,method.dist="correlation", nboot=boots,parallel=T)
#pdf("pfam_clustering_all.pdf", width=widthp, height=heightp)
#plot(pfam.pv, cex=1, cex.pv=0.8, main=paste("PFAM profile similarity clustering. method: ", method,sep=""))
#dev.off()

######## interpro
#print("calculating Interpro similarity..." )
#rownames(all_interpro) <- all_interpro$X
#all_interpro$X <- NULL
#all_interpro$description <- NULL
#all_interpro <- all_interpro[,-ncol(all_interpro)]
#all_interpro <- t(all_interpro)
#rownames(all_interpro)


#remove columns with zero variance, as they will not help with the clustering and pvclust does not like them
#variances<-apply(t(all_interpro), 1, var)
#low_var_names <- names(variances[which(variances<=variance)])
#all_interpro_subset <- all_interpro[, !colnames(all_interpro) %in% low_var_names]

#interpro.pv <- pvclust(t(all_interpro_subset),method.hclust=method,method.dist="correlation", nboot=boots, parallel=T)
#interpro.pv <- pfam.pv
#pdf("interpro_clustering_all.pdf", width=widthp, height=heightp)
#plot(interpro.pv, cex=1, cex.pv=0.8, main=paste("InterPro profile similarity clustering. method: ", method,sep=""))
#dev.off()

