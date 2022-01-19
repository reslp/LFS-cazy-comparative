library(phytools)

args <- commandArgs(trailingOnly=TRUE)

print("load input")
tree <- read.tree(args[1])
data <- read.csv(args[2], stringsAsFactors = F)
outdir <- args[3]

#parameters for cutoff:
kcutoff <- 1
pvalue <- 0.05

rownames(data) <- data$X
data$X <- NULL
data <- t(data)

options(scipen=999) # this is to omit scientific notation in values.
print("Calculating phylogenetic signal")
df <- data.frame(cazyme=character(),K=numeric(),pvalueK=numeric(),lambda=numeric(),pvalueL=numeric(),logL=numeric(), logL0=numeric(), stringsAsFactors=FALSE)
for (i in 1:length(colnames(data))) {
  print(i)
  which <- colnames(data)[i]
  sig <- phylosig(tree, log(data[,i]+1), method="K", test=TRUE, nsim=10000)
  sig2 <- phylosig(tree, log(data[,i]+1), method="lambda", test=TRUE)
  dat <- c(cazyme=which, K=sig$K, pvalueK=sig$P, lambda=sig2$lambda, pvalueL=sig2$P, logL=sig2$logL, logL0=sig2$logL0)
  df[i,] <- dat
}

print("Write output")
write.table(df, file=paste(outdir,"phylosig_cazymes.txt", sep=""), sep=",", row.names = FALSE, quote=FALSE)
write.table(df[df$K > kcutoff & df$pvalueK < pvalue, ], file=paste(outdir, "phylosig_sign_cazymes.txt", sep=""), sep=",", row.names = FALSE, quote=FALSE)

