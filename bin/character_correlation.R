library(phytools)
#library(tidyverse)
#library(ggplot2)
library(parallel)
args <- commandArgs(trailingOnly=TRUE)

wd <- args[1]
treefile <- args[2]
discrete_chars <- args[3]
cazyme_chars <- args[4]
outdir <- args[5]
threads <- args[6]



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
ngen <- 2000000 # number of generations
burnin <- 0.2 * ngen # 20% burn-in
sample <- 200 # sample all n generations
quiet <- T
run_mcmc <- function(i, char) {
  subsampled_data <- data[,c(char, all_cazy_groups[i])]
  colnames(subsampled_data)[1] <- "my_char"
  subsampled_data$my_char <- as.character(subsampled_data$my_char)
  combination <- c(char, all_cazy_groups[i])
  cat(format(Sys.time(), "%a %b %d %X %Y"))
  cat(paste(" - ", toString(i)," - ",all_cazy_groups[i],"\n", sep=""))
  cat ("          Chain1\n")
  mcmc1 <- threshBayes(tree, subsampled_data, types=c("discrete", "continuous"), control=list(sample=sample, quiet=quiet),ngen=ngen)
  cat ("          Chain2\n")
  mcmc2 <- threshBayes(tree, subsampled_data, types=c("discrete", "continuous"), control=list(sample=sample, quiet=quiet),ngen=ngen)
  cat ("          Chain3\n")
  mcmc3 <- threshBayes(tree, subsampled_data, types=c("discrete", "continuous"), control=list(sample=sample, quiet=quiet),ngen=ngen)
  cat ("          Chain4\n")
  mcmc4 <- threshBayes(tree, subsampled_data, types=c("discrete", "continuous"), control=list(sample=sample, quiet=quiet),ngen=ngen)
  mcmc <- list()
  mcmc$par <- rbind(mcmc1$par, mcmc2$par, mcmc3$par, mcmc4$par)
  mcmc$liab <- rbind(mcmc1$liab, mcmc2$liab, mcmc3$liab, mcmc4$liab)
  mcmc$burnin <- mcmc1$burnin
  mcmc$levels <- mcmc1$levels
  mcmc$types <- mcmc1$types
  class(mcmc) <- "threshBayes"
  mcmc$combination <- combination
  return(mcmc)
}

num_cazy <- seq(1, length(all_cazy_groups))
print(detectCores())
print(detectCores(logical = FALSE))
#mcmcs <- mclapply(num_cazy,run_mcmc, mc.cores=threads)
#mcmcs <- mclapply(num_cazy,run_mcmc, mc.cores=4)
#mcmcs <- list()
x <- 1
for (char in colnames(data1)) {
	cat(paste("Running analysis for ", char, "\n", sep=""))
	mcmcs <- mclapply(num_cazy,run_mcmc,char, mc.cores=12)
	cat("Saving environment\n")
	save.image(paste(outdir, "/data_",char ,".RData", sep=""))
	#x <- x + 1
}
#mcmcs <- mclapply(num_cazy,run_mcmc, mc.cores=1)
#pdf(file=paste(outdir,"/r_values_overview.pdf", sep=""))
#for (i in 1:length(all_cazy_groups)) {
#  cat(format(Sys.time(), "%a %b %d %X %Y"))
#  cat(paste(" - ", toString(i)," - ",all_cazy_groups[i],"\n", sep=""))
#  subsampled_data <- data[,c("lichen", all_cazy_groups[i])]
#  subsampled_data$lichen <- as.character(subsampled_data$lichen)
#  combination[[i]] <- c("lichen", all_cazy_groups[i])
#  mcmcs[[i]] <- threshBayes(tree, subsampled_data, types=c("discrete", "continuous"), control=list(sample=500, quiet=T),ngen=1000000)
#  plot(mcmcs[[i]], bw=0.1)
#  lines(rep(0.7,2),c(0,par()$usr[4]),lwd=2,lty="dashed",col="red")
#  text(x=0.7,y=0.9*par()$usr[4],"simulated r",pos=4,cex=0.7)
#  tt <- paste(combination[[i]], sep = " and ", collapse=" and ")
#  title(tt)
#}
#dev.off()


#save.image(paste(outdir, "/data.RData", sep=""))

