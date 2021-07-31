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
print(lecanoromycetes_names)

ostropomycetidae_all_names <- ostropomycetidae_names
ostropales_names <- c("Graphis_scripta", "Stictis_urceolatum", "Gomphillus_americanus", "Thelotrema_lepadinum", "Cyanodermella_asteris")
ostropomycetidae_names <- ostropomycetidae_names[!ostropomycetidae_names %in% ostropales_names]

non_lecanoromycetes_names <- taxonomy_info[!taxonomy_info$name %in% lecanoromycetes_names,]$name

# formula to calculate difference between numbers:
calc_diff <- function(V1, V2) {
	diff <- (abs(V1 - V2) / ((V1 + V2) / 2)) * 100
	return(diff)
}


cats <- c()
sum_lecanoromycetes_counts <- c()
sum_lecanoro_counts <- c()
sum_ostropo_all_counts <- c()
sum_ostropo_counts <- c()
sum_ostropales_counts <- c()
percent_dif_ostropo <- c()
percent_dif_ostropo_all <- c()
percent_dif_ostropo_lecanoro <- c()
percent_dif_ostropales <- c()
percent_dif_lecanoro <- c()
percent_dif_lecanoromycetes <- c()
sum_other_asco_counts <- c()
pvalues_ostropo <- c()
pvalues_ostropo_all <- c()
pvalues_ostropo_lecanoro <- c()
pvalues_ostropales <- c()
pvalues_lecanoro <- c()
pvalues_lecanoromycetes <- c()

for (cat in colnames(cazy_sum_data)) {
  cats <- c(cats, cat)
 
  # subset the values for calculations: 
  lecanoro_counts <- cazy_sum_data[rownames(cazy_sum_data) %in% lecanoromycetidae_names,][,cat]
  lecanoromycetes_counts <- cazy_sum_data[rownames(cazy_sum_data) %in% lecanoromycetes_names,][,cat]
  ostropo_counts <- cazy_sum_data[rownames(cazy_sum_data) %in% ostropomycetidae_names,][,cat]
  ostropo_all_counts <- cazy_sum_data[rownames(cazy_sum_data) %in% ostropomycetidae_all_names,][,cat]
  ostropales_counts <- cazy_sum_data[rownames(cazy_sum_data) %in% ostropales_names,][,cat]
  
  # calculate mean values for cazyme numbers
  # Lecanoromycetidae
  sum_lecanoro_counts <- c(sum_lecanoro_counts,mean(lecanoro_counts))
  # Lecanoromycetes
  sum_lecanoromycetes_counts <- c(sum_lecanoromycetes_counts,mean(lecanoromycetes_counts))
  # Ostropomycetidae without OG
  sum_ostropo_counts <- c(sum_ostropo_counts,mean(ostropo_counts))
  # Ostropomycetidae
  sum_ostropo_all_counts <- c(sum_ostropo_all_counts,mean(ostropo_all_counts))
  # OG clade
  sum_ostropales_counts <- c(sum_ostropales_counts,mean(ostropales_counts))
  
  # subset and calculate for other fungi
  other_asco_counts <- cazy_sum_data[rownames(cazy_sum_data) %in% non_lecanoromycetes_names,][,cat]
  sum_other_asco_counts <- c(sum_other_asco_counts, mean(other_asco_counts))
  
  # caculate the % difference between mean values found in different tax groups:
  # Ostropomycetidae without OG vs other fungi
  percent_dif_ostropo <- c(percent_dif_ostropo, calc_diff(mean(ostropo_counts), mean(other_asco_counts)))
  # Ostropomycetidae vs other fungi
  percent_dif_ostropo_all <- c(percent_dif_ostropo_all, calc_diff(mean(ostropo_all_counts), mean(other_asco_counts)))
  # Ostropomycetidae vs. other fungi
  percent_dif_ostropo_lecanoro <- c(percent_dif_ostropo_lecanoro, calc_diff(mean(ostropo_all_counts), mean(lecanoro_counts)))
  # OG clade vs. other fungi
  percent_dif_ostropales <- c(percent_dif_ostropales,calc_diff(mean(ostropales_counts),mean(other_asco_counts)))
  # Lecanoromycetidae vs. other fungi
  percent_dif_lecanoro <- c(percent_dif_lecanoro,calc_diff(mean(lecanoro_counts),mean(other_asco_counts)))
  # Lecanoromycetes vs. other fungi
  percent_dif_lecanoromycetes <- c(percent_dif_lecanoromycetes,calc_diff(mean(lecanoromycetes_counts),mean(other_asco_counts)))
  
  #compare Lecanoromycetes with other fungi
  test0 <- wilcox.test(lecanoromycetes_counts, other_asco_counts)
  pvalues_lecanoromycetes <- c(pvalues_lecanoromycetes, test0$p.value)
  
  #compare Lecanoromycetidae with other fungi
  test1 <- wilcox.test(lecanoro_counts, other_asco_counts)
  pvalues_lecanoro <- c(pvalues_lecanoro, test1$p.value)
  
  #compare Ostropomycetidae without OG with other fungi
  test2 <- wilcox.test(ostropo_counts, other_asco_counts)
  pvalues_ostropo <- c(pvalues_ostropo, test2$p.value)
  
  #compare OG clade counts with other fungi
  test3 <- wilcox.test(ostropales_counts, other_asco_counts)
  pvalues_ostropales <- c(pvalues_ostropales, test3$p.value)
  
  #compare the counts of Ostropomycetidae and Lecanoromycetidae
  test4 <- wilcox.test(ostropo_all_counts, lecanoro_counts)
  pvalues_ostropo_lecanoro <- c(pvalues_ostropo_lecanoro, test4$p.value)
  
  #compare the counts of Ostropomycetidae and other fungi
  test5 <- wilcox.test(ostropo_all_counts, other_asco_counts)
  pvalues_ostropo_all <- c(pvalues_ostropo_all, test5$p.value)
  #print(paste(cat, sum(lecanoro_counts), sum(other_asco_counts), percent_dif, test$p.value, sep=" "))
}

print(sum_lecanoro_counts)
print(sum_lecanoromycetes_counts)
print(sum_ostropo_counts)
print(sum_other_asco_counts)
print(percent_dif_lecanoro)
print(percent_dif_lecanoromycetes)
print(percent_dif_ostropo)
print(pvalues_ostropo)
print(pvalues_lecanoro)
print(pvalues_lecanoromycetes)
print(cats)

df <- data.frame(
	category = cats,
	Lecanoromycetes=sum_lecanoromycetes_counts,
	Lecanoromycetidae = sum_lecanoro_counts,
	Ostropomycetidae= sum_ostropo_all_counts,
	Ostropomycetidae_noOG = sum_ostropo_counts,
	Ostropales=sum_ostropales_counts,
	Other_fungi = sum_other_asco_counts,
	Percent_Difference_Lecanoromycetes=percent_dif_lecanoromycetes,
	Percent_Difference_Lecanoromycetidae = percent_dif_lecanoro,
	Percent_Difference_Ostropomycetidae= percent_dif_ostropo_all,
	Percent_Difference_Ostropomycetidae_noOG= percent_dif_ostropo,
	Percent_Difference_Ostropales= percent_dif_ostropales,
	Percent_Difference_OstroLecano= percent_dif_ostropo_lecanoro,
	Pvalue_lecanoromycetes = pvalues_lecanoromycetes,
	Pvalue_lecanoro = pvalues_lecanoro,
	Pvalue_ostropo = pvalues_ostropo_all,
	Pvalue_ostropo_noOG = pvalues_ostropo,
	Pvalue_ostropales=pvalues_ostropales,
	Pvalue_ostropo_lecanoro = pvalues_ostropo_lecanoro)

write.csv(df,file="overview_cazyme_counts.csv", quote=F, row.names = F)
