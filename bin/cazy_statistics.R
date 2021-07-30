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
ostropales_names <- c("Graphis_scripta", "Stictis_urceolatum", "Gomphillus_americanus", "Thelotrema_lepadinum", "Cyanodermella_asteris")
ostropomycetiade_names <- ostropomycetidae_names[!ostropomycetidae_names %in% ostropales_names]


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
sum_ostropo_counts <- c()
sum_ostropales_counts <- c()
percent_dif_ostropo <- c()
percent_dif_ostropales <- c()
percent_dif_lecanoro <- c()
sum_other_asco_counts <- c()
pvalues_ostropo <- c()
pvalues_ostropales <- c()
pvalues_lecanoro <- c()

for (cat in colnames(cazy_sum_data)) {
  cats <- c(cats, cat)
  
  lecanoro_counts <- cazy_sum_data[rownames(cazy_sum_data) %in% lecanoromycetidae_names,][,cat]
  ostropo_counts <- cazy_sum_data[rownames(cazy_sum_data) %in% ostropomycetidae_names,][,cat]
  ostropales_counts <- cazy_sum_data[rownames(cazy_sum_data) %in% ostropales_names,][,cat]
  
  sum_lecanoro_counts <- c(sum_lecanoro_counts,mean(lecanoro_counts))
  sum_ostropo_counts <- c(sum_ostropo_counts,mean(ostropo_counts))
  sum_ostropales_counts <- c(sum_ostropales_counts,mean(ostropales_counts))
  
  other_asco_counts <- cazy_sum_data[rownames(cazy_sum_data) %in% non_lecanoromycetes_names,][,cat]
  sum_other_asco_counts <- c(sum_other_asco_counts, mean(other_asco_counts))
  
  percent_dif_ostropo <- c(percent_dif_ostropo,100-(100/mean(other_asco_counts)*mean(ostropo_counts)))
  percent_dif_ostropales <- c(percent_dif_ostropales,100-(100/mean(other_asco_counts)*mean(ostropales_counts)))
  percent_dif_lecanoro <- c(percent_dif_lecanoro,100-(100/mean(other_asco_counts)*mean(lecanoro_counts)))
  
  test1 <- wilcox.test(lecanoro_counts, other_asco_counts)
  pvalues_lecanoro <- c(pvalues_lecanoro, test1$p.value)
  
  test2 <- wilcox.test(ostropo_counts, other_asco_counts)
  pvalues_ostropo <- c(pvalues_ostropo, test2$p.value)
  
  test3 <- wilcox.test(ostropales_counts, other_asco_counts)
  pvalues_ostropales <- c(pvalues_ostropales, test3$p.value)
  #print(paste(cat, sum(lecanoro_counts), sum(other_asco_counts), percent_dif, test$p.value, sep=" "))
}

print(sum_lecanoro_counts)
print(sum_ostropo_counts)
print(sum_other_asco_counts)
print(percent_dif_lecanoro)
print(percent_dif_ostropo)
print(pvalues_ostropo)
print(pvalues_lecanoro)
print(cats)
df <- data.frame(category = cats, Lecanoromycetidae = sum_lecanoro_counts, Ostropomycetidae= sum_ostropo_counts, Ostropales=sum_ostropales_counts, Other_fungi = sum_other_asco_counts, Percent_Difference_Lecanoromycetidae = percent_dif_lecanoro, Percent_Difference_Ostropomycetidae= percent_dif_ostropo, Percent_Difference_Ostropales= percent_dif_ostropales, Pvalue_lecanoro = pvalues_lecanoro, Pvalue_ostropo = pvalues_ostropo, Pvalue_ostropales=pvalues_ostropales)

write.csv(df,file="overview_cazyme_counts.csv", quote=F, row.names = F)
