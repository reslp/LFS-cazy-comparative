library(phytools)
#library(tidyverse)
#library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

wd <- args[1]
treefile <- args[2]
discrete_chars <- args[3]
cazyme_chars <- args[4]
outdir <- args[5]




#setwd("/home/reslp/Dropbox/Philipp/xylographa_comparative_genomics/")
#setwd("/Users/sinnafoch/Dropbox/Philipp/xylographa_comparative_genomics/")
#treefile <- "tmp/correlation_test/82_genomes_ultra.tre"
#discrete_chars <- "tmp/correlation_test/character_information.csv"
#cazyme_chars <- "tmp/correlation_test/CAZyme.all.results.csv"
#outdir <- "tmp/correlation_test"

#load the tree
tree <- read.tree(treefile)

# reading raw character data and do some reformatting to be able to combine them into a single dataframe
data1 <- read.csv(discrete_chars)
rownames(data1) <- data1$species
data1$species <- NULL
data2 <- read.csv(cazyme_chars,stringsAsFactors=FALSE)

rownames(data2) <- data2$X
data2$X <- NULL
data2 <- t(data2)
data <- cbind(data1, data2)
all_cazy_groups <- colnames(data2)

combination <- list()
mcmcs <- list()
pdf(file=paste(outdir,"/r_values_overview.pdf", sep=""))
for (i in 1:length(all_cazy_groups)) {
  cat(format(Sys.time(), "%a %b %d %X %Y"))
  cat(paste(" - ", toString(i)," - ",all_cazy_groups[i],"\n", sep=""))
  subsampled_data <- data[,c("lichen", all_cazy_groups[i])]
  subsampled_data$lichen <- as.character(subsampled_data$lichen)
  combination[[i]] <- c("lichen", all_cazy_groups[i])
  mcmcs[[i]] <- threshBayes(tree, subsampled_data, types=c("discrete", "continuous"), control=list(sample=500, quiet=T),ngen=1000000)
  plot(mcmcs[[i]], bw=0.1)
  lines(rep(0.7,2),c(0,par()$usr[4]),lwd=2,lty="dashed",col="red")
  text(x=0.7,y=0.9*par()$usr[4],"simulated r",pos=4,cex=0.7)
  tt <- paste(combination[[i]], sep = " and ", collapse=" and ")
  title(tt)
}
dev.off()
save.image(paste(outdir, "/data.RData", sep=""))

