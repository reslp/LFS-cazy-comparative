library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)

wd <- args[1]
taxinfo <- args[2]
cazyinfo <- args[3]

setwd(wd)

options(StringsAsFactors = F)
options(scipen=999)
cazy_sum_data <- read.table(cazyinfo, sep=",", header=T)
rownames(cazy_sum_data) <- cazy_sum_data$X
cazy_sum_data$X <- NULL
cazy_sum_data <- t(cazy_sum_data)



taxonomy_info <- read.table(taxinfo, sep=",", header=T)
taxonomy_info


lecanoromycetidae_names <- taxonomy_info[taxonomy_info$class == "Lecanoromycetidae",]$name
ostropomycetidae_names <- taxonomy_info[taxonomy_info$class == "Ostropomycetidae",]$name

lecanoromycetes_names <- c(ostropomycetidae_names, lecanoromycetidae_names)
non_lecanoromycetes_names <- taxonomy_info[!taxonomy_info$name %in% lecanoromycetes_names,]$name

total_number_lecanoromycetes <- sum(colSums(cazy_sum_data[rownames(cazy_sum_data) %in% lecanoromycetes_names,]))
total_number_nonlecanoromycetes <- sum(colSums(cazy_sum_data[rownames(cazy_sum_data) %in% non_lecanoromycetes_names,]))

percent_diff_lecanoro_non_lecanoro <- 100 - (100 / total_number_nonlecanoromycetes * total_number_lecanoromycetes)

total_number_lecanoromycetidae <- sum(colSums(cazy_sum_data[rownames(cazy_sum_data) %in% lecanoromycetidae_names,]))
total_number_ostropomycetidae <- sum(colSums(cazy_sum_data[rownames(cazy_sum_data) %in% ostropomycetidae_names,]))

percent_diff_lecanoro_vs_ostropo <- 100 - (100 / total_number_ostropomycetidae * total_number_lecanoromycetidae)
percent_diff_lecanoro_vs_ostropo



### GHs
total_GHs_lecanoromycetes <- colSums(cazy_sum_data[rownames(cazy_sum_data) %in% lecanoromycetes_names,])["GH"]
total_GHs_non_lecanoromycetes <- colSums(cazy_sum_data[rownames(cazy_sum_data) %in% non_lecanoromycetes_names,])["GH"]

gh_percent_diff_lecanoro_non_lecanoro <- 100 - (100 / total_GHs_non_lecanoromycetes * total_GHs_lecanoromycetes)


total_GHs_lecanoromycetidae <- colSums(cazy_sum_data[rownames(cazy_sum_data) %in% lecanoromycetidae_names,])["GH"]
total_GHs_ostropomycetidae <- colSums(cazy_sum_data[rownames(cazy_sum_data) %in% ostropomycetidae_names,])["GH"]

gh_percent_diff_lecanoro_vs_ostropo <- 100 - (100 / total_GHs_ostropomycetidae * total_GHs_lecanoromycetidae)

cazy_sum_data[rownames(cazy_sum_data) %in% lecanoromycetidae_names,][,"AA"]

cats <- c()
sum_lecanoro_counts <- c()
percent_dif <- c()
sum_other_asco_counts <- c()
pvalues <- c()
for (cat in colnames(cazy_sum_data)) {
  cats <- c(cats, cat)
  lecanoro_counts <- cazy_sum_data[rownames(cazy_sum_data) %in% lecanoromycetes_names,][,cat]
  sum_lecanoro_counts <- c(sum_lecanoro_counts,sum(lecanoro_counts))
  other_asco_counts <- cazy_sum_data[rownames(cazy_sum_data) %in% non_lecanoromycetes_names,][,cat]
  sum_other_asco_counts <- c(sum_other_asco_counts, sum(other_asco_counts))
  
  percent_dif <- c(percent_dif,100-(100/sum(other_asco_counts)*sum(lecanoro_counts)))
  test <- wilcox.test(lecanoro_counts, other_asco_counts)
  pvalues <- c(pvalues, test$p.value)
  #print(paste(cat, sum(lecanoro_counts), sum(other_asco_counts), percent_dif, test$p.value, sep=" "))
}


df <- data.frame(category = cats, Lecanoromycetes = sum_lecanoro_counts, Other_fungi = sum_other_asco_counts, Percent_Difference = percent_dif, Pvalue = pvalues)

write.csv(df,file="overview_cazyme_counts.csv", quote=F, row.names = F)
