
library(ggtree)
library(tidyverse)
library(RColorBrewer)
library(cowplot)

args <- commandArgs(trailingOnly=TRUE)

wd <- args[1]
treefile <- args[2]
fam_data <- args[3]
additional_data <- args[4]
out_prefix <- args[5]
family <- args[6]
deeploc_file <- args[7]

prob <- 0.7
#setwd("/home/reslp/Dropbox/Philipp/xylographa_comparative_genomics/tmp/saccharis_trees")
#setwd("/Users/sinnafoch/Dropbox/Philipp/xylographa_comparative_genomics/tmp/saccharis_trees")


## load code from ggnewscale. The package is not available in conda so I could not include it in the environment.
# The code here was downloaded directly from github: https://github.com/eliocamp/ggnewscale
# Download was on June 2nd, 2020. The downloaded commit is: 16f0331

#first check if rData file already exists, in which case plotting will be skipped:
if (file.exists(paste(out_prefix,"/",family,".rData", sep="")) == TRUE) {
	print(paste("File for ", family, " already exists. Plot will be skipped"))
	quit()
}


new_scale <- function(new_aes) {
  structure(standardise_aes_names(new_aes), class = "new_aes")
}

new_scale_fill <- function() {
  new_scale("fill")
}


new_scale_color <- function() {
  new_scale("colour")
}


new_scale_colour <- function() {
  new_scale("colour")
}


ggplot_add.new_aes <- function(object, plot, object_name) {
  # To add default scales (I need to build the whole plot because they might be computed aesthetics)
  if (is.null(plot$scales$get_scales(object))) {
    plot$scales <- ggplot2::ggplot_build(plot)$plot$scales
  }
  # Global aes
  old_aes <- names(plot$mapping)[remove_new(names(plot$mapping)) %in% object]
  new_aes <- paste0(old_aes, "_new")
  names(plot$mapping)[names(plot$mapping) == old_aes] <- new_aes
  
  
  
  plot$layers <- bump_aes_layers(plot$layers, new_aes = object)
  plot$scales$scales <- bump_aes_scales(plot$scales$scales, new_aes = object)
  plot$labels <- bump_aes_labels(plot$labels, new_aes = object)
  plot
}

bump_aes_layers <- function(layers, new_aes) {
  lapply(layers, bump_aes_layer, new_aes = new_aes)
  
}

bump_aes_layer <- function(layer, new_aes) {
  original_aes <- new_aes
  
  new_layer <- ggplot2::ggproto(NULL, layer)
  
  # Get explicit mapping
  old_aes <- names(new_layer$mapping)[remove_new(names(new_layer$mapping)) %in% new_aes]
  
  # If not explicit, get the default
  if (length(old_aes) == 0) {
    old_aes <- names(new_layer$stat$default_aes)[remove_new(names(new_layer$stat$default_aes)) %in% new_aes]
    if (length(old_aes) == 0) {
      old_aes <- names(new_layer$geom$default_aes)[remove_new(names(new_layer$geom$default_aes)) %in% new_aes]
    }
  }
  new_aes <- paste0(old_aes, "_new")
  
  old_geom <- new_layer$geom
  
  old_setup <- old_geom$handle_na
  new_setup <- function(self, data, params) {
    colnames(data)[colnames(data) %in% new_aes] <- original_aes
    old_setup(data, params)
  }
  
  new_geom <- ggplot2::ggproto(paste0("New", class(old_geom)[1]), old_geom,
                               handle_na = new_setup)
  
  new_geom$default_aes <- change_name(new_geom$default_aes, old_aes, new_aes)
  new_geom$non_missing_aes <- change_name(new_geom$non_missing_aes, old_aes, new_aes)
  new_geom$required_aes <- change_name(new_geom$required_aes, old_aes, new_aes)
  new_geom$optional_aes <- change_name(new_geom$optional_aes, old_aes, new_aes)
  
  new_layer$geom <- new_geom
  
  old_stat <- new_layer$stat
  
  old_setup2 <- old_stat$handle_na
  new_setup <- function(self, data, params) {
    colnames(data)[colnames(data) %in% new_aes] <- original_aes
    old_setup2(data, params)
  }
  
  new_stat <- ggplot2::ggproto(paste0("New", class(old_stat)[1]), old_stat,
                               handle_na = new_setup)
  
  new_stat$default_aes <- change_name(new_stat$default_aes, old_aes, new_aes)
  new_stat$non_missing_aes <- change_name(new_stat$non_missing_aes, old_aes, new_aes)
  new_stat$required_aes <- change_name(new_stat$required_aes, old_aes, new_aes)
  new_stat$optional_aes <- change_name(new_stat$optional_aes, old_aes, new_aes)
  
  new_layer$stat <- new_stat
  
  new_layer$mapping <- change_name(new_layer$mapping, old_aes, new_aes)
  new_layer$aes_params <- change_name(new_layer$aes_params, old_aes, new_aes)
  new_layer
}

bump_aes_scales <- function(scales, new_aes) {
  lapply(scales, bump_aes_scale, new_aes = new_aes)
}


bump_aes_scale <- function(scale, new_aes) {
  old_aes <- scale$aesthetics[remove_new(scale$aesthetics) %in% new_aes]
  if (length(old_aes) != 0) {
    new_aes <- paste0(old_aes, "_new")
    
    scale$aesthetics[scale$aesthetics %in% old_aes] <- new_aes
    
    no_guide <- isFALSE(scale$guide) | isTRUE(scale$guide == "none")
    if (!no_guide) {
      if (is.character(scale$guide)) {
        scale$guide <- get(paste0("guide_", scale$guide), mode = "function")()
      }
      scale$guide$available_aes[scale$guide$available_aes %in% old_aes] <- new_aes
    }
  }
  
  scale
}

bump_aes_labels <- function(labels, new_aes) {
  old_aes <-  names(labels)[remove_new(names(labels)) %in% new_aes]
  new_aes <- paste0(old_aes, "_new")
  
  names(labels)[names(labels) %in% old_aes] <- new_aes
  labels
}


change_name <- function(list, old, new) {
  UseMethod("change_name")
}

change_name.character <- function(list, old, new) {
  list[list %in% old] <- new
  list
}

change_name.default <- function(list, old, new) {
  nam <- names(list)
  nam[nam %in% old] <- new
  names(list) <- nam
  list
}

change_name.NULL <- function(list, old, new) {
  NULL
}


remove_new <- function(aes) {
  gsub("(_new)*", "", aes, fixed = FALSE)
  # stringi::stri_replace_all(aes, "", regex = "(_new)*")
}


isTRUE <- function (x) {
  is.logical(x) && length(x) == 1L && !is.na(x) && x
}

isFALSE <- function (x) {
  is.logical(x) && length(x) == 1L && !is.na(x) && !x
}


## this is the end of ggnewscale code

#load data
tree <- read.tree(treefile)
print("Load cazy data")
cazy_data <- read.csv(fam_data, sep="\t", header=T, stringsAsFactors = F)
print("load additional mapping data")
additional_mapping <- read.csv(additional_data, header=T, stringsAsFactors = F, sep="\t")
deeploc <- read.csv(deeploc_file, header=T, sep="\t", stringsAsFactors = F)

#reformat deeploc data:
ids <- deeploc$ID
#strsplit(strsplit(ids[[147]], " ")[[1]][2]), "_(?=[^_]+$)", perl=TRUE))[[1]][1]
for (i in 1:length(deeploc$ID)){ # keep only first number if there are more than one, only the first number is downloaded
  if (startsWith(strsplit(ids[i]," ")[[1]][1], "00")) { # account for different naming in deeploc output resulting from different name in saccharis alignment
    ids[i] <- strsplit(strsplit(ids[[i]], " ")[[1]][2],"_(?=[^_]+$)", perl=T)[[1]][1] # split string at last occurence of _ (some NCBI names have _ in there names)
  }else {
    ids[i] <- strsplit(ids[i]," ")[[1]][1]
  }
  
}

rownames(deeploc) <- ids

## this is to make sure tip labels and dataframe rows match
## it is necessary because saccharis extracts only cazy motifs from unkown sequences and there could be more than 
## on such motive in a single sequence
## eg. U00647651 (a sequence from AA2) is represented twice in the tree as: U00647651a U00647651b
## most of the annotation is however sequence specific and not motif specific so it will be duplicated
## information about the function about the fragment can be derived from the pyhlogenetic context and sister sequences
## therefore the dataframe is extended with new rows using information from the gene based annotations:

#remove <- vector()
print("Duplicating correct dataframe rows containing deeploc location information")
for (i in 1:length(rownames(deeploc))) {
  for (j in 1:length(tree$tip.label)){
    if (grepl(rownames(deeploc)[i], tree$tip.label[j], fixed=T) && (rownames(deeploc)[i] != tree$tip.label[j])){
      print(tree$tip.label[j])
      deeploc <- rbind(deeploc,deeploc[i,])
      rownames(deeploc)[length(rownames(deeploc))] <- tree$tip.label[j]
      #remove <- c(remove, rownames(deeploc)[i])
    }
  }
}

#deeploc
rownames(deeploc) %in% tree$tip.label
# this should be empty:
tree$tip.label[!tree$tip.label %in% rownames(deeploc)]


## this is to subsample only ids above a certain probability. Some assignment probabilities are very low
best_value <- vector()
for (i in 1:length(rownames(deeploc))) {
  column <- gsub("/", ".",deeploc$Location[i]) # account for different labeling of columns introduced by r
  
  best_value <- c(best_value, deeploc[,column][i])
}
deeploc$best_value <- best_value
deeploc <- deeploc[deeploc$best_value >= prob,] # only keep labels with prob > 70%



#format mapping data
rownames(additional_mapping) <- additional_mapping$saccharis_name

## similar to what is done above for the deeploc data, taxonomy information needs to be transfered for duplicated 
## sequences as well. This is done here. The dataframe is huge in this case and this chunk can run for some time...

print("Duplicating correct dataframe rows containing additional mapping (eg. taxonomy)")
for (i in 1:length(rownames(additional_mapping))) {
  for (j in 1:length(tree$tip.label)){
    if (grepl(rownames(additional_mapping)[i], tree$tip.label[j], fixed=T) && (rownames(additional_mapping)[i] != tree$tip.label[j])){
      print(tree$tip.label[j])
      additional_mapping <- rbind(additional_mapping,additional_mapping[i,])
      rownames(additional_mapping)[length(rownames(additional_mapping))] <- tree$tip.label[j]
    }
  }
}

print("Some reformating of cazy data for easier handling later on")
# do some reformatting to cazy data to make it easier to handle later on
cl <- colnames(cazy_data)
cazy_data$sub
cl[1] <- "Family"
cl[2] <- "Domain"
colnames(cazy_data) <- cl
cazy_data <- cazy_data[cazy_data$GenBank != "",]

#format ncbi accession numbers correctly:
ncbi_names <- cazy_data$GenBank
ec_names <- cazy_data$EC.
for (i in 1:length(ncbi_names)){ # keep only first number if there are more than one, only the first number is downloaded
  ncbi_names[i] <- strsplit(ncbi_names[i]," ")[[1]][1]
  ec_names[i] <- strsplit(ec_names[i]," ")[[1]][1] 
}
if (length(rownames(cazy_data))>1) {
	cazy_data <- as.data.frame(sapply(cazy_data, as.character))
}

ncbi_names<-make.unique(ncbi_names,sep="_") # make accessions unique, this is not yet compatible with the tip labels (which contain letters)
cazy_data$GenBank <- ncbi_names
rownames(cazy_data) <- ncbi_names
cazy_data$EC. <- ec_names

#combine datasets before subsampling for plotting

#subsample data to plot

if ("Subf" %in% colnames(cazy_data)) {
	plot_data <- cazy_data %>% select("EC.", "Subf", "Domain")
} else {
	plot_data <- cazy_data %>% select("EC.", "Domain")
	}
#plot_data$Domain

names <- c(rownames(plot_data),rownames(additional_mapping)) 

taxonomy <- c(as.vector(plot_data$Domain),additional_mapping$class)
names(taxonomy) <- names
taxonomy_plotting <- taxonomy[names(taxonomy) %in% tree$tip.label]
taxonomy_plotting <-data.frame(t(data.frame(as.list(taxonomy_plotting))))
colnames(taxonomy_plotting) <- "taxonomy"
#taxonomy_plotting$taxonomy




#create number of distinct colors according to the plot data


# subfamily and EC code colors need to be created dynamically because the numers will always be different
if ("Subf" %in% colnames(cazy_data)) {
	n <- length(unique(plot_data$Subf))
	subf_colors <- colorRampPalette(brewer.pal(8, "BrBG"))(n)
	length(subf_colors)
} 

n <- length(unique(plot_data$EC.))
ec_colors <- colorRampPalette(brewer.pal(11, "Spectral"))(n)
length(ec_colors)

# colors for taxonomy and location will be hardcoded so that they are the same accross multiple plots
num_cat_tax <- c("Archaea", "Bacteria","Viruses", "Eukaryota", "Lecanoromycetes", "Leotiomycetes", "Sordariomycetes", "Arthoniomycetes", "Dothideomycetes", "Eurotiomycetes", "unclassified")
tax_colors <- c("#d9d9d9", "#969696", "#f0f0f0")
tax_colors <- c(tax_colors, brewer.pal(length(num_cat_tax)-4, "Spectral"))
tax_colors <- c(tax_colors, "#ffffff")
names(tax_colors) <- num_cat_tax


num_cat_deeploc <- c("Membrane", "Nucleus", "Cytoplasm", "Extracellular", "Mitochondrion", "Cell_membrane", "Endoplasmic_reticulum", "Plastid", "Golgi_apparatus", "Lysosome/Vacuole", "Peroxisome")
deeploc_colors <- brewer.pal(length(num_cat_deeploc), "Paired")
names(deeploc_colors)<- num_cat_deeploc


#if ("Subf" %in% colnames(cazy_data)) {
#	all_colors <- c(subf_colors, ec_colors, tax_colors)
#	names(all_colors) <- c(unique(plot_data$Subf), unique(plot_data$EC.), as.vector(unique(taxonomy_plotting$taxonomy))) 
#} else { all_colors <- c(ec_colors, tax_colors)
#	names(all_colors) <- c(unique(plot_data$EC.), as.vector(unique(taxonomy_plotting$taxonomy)))
#}

# plot tree
ggt <- ggtree(tree, layout="circular", size=0.2)
ggt <- ggt + xlim(-2, NA)
#ggt <- ggt + xlim(-2, NA)
ggt <- open_tree(ggt, 10)
ggt <- rotate_tree(ggt, 10)
#ggt +ggtitle("bla")+theme(plot.title = element_text(hjust = 0.5, vjust=-55))

if ("Subf" %in% colnames(cazy_data)) {

	subfp <- gheatmap(ggt, plot_data[,"Subf",drop=F], offset = 0.1, width=0.1, color=NULL, colnames_position="top", colnames_angle=90, colnames_offset_y=0, hjust=0, font.size=2)
	subfp  <- subfp  + scale_fill_manual("Subfamily", values=subf_colors) 

	ecp <- subfp + new_scale_fill()
	ecp <- gheatmap(ecp, plot_data[,"EC.",drop=F], offset = 1.1, width=0.1, color=NULL, colnames_position="top", colnames_angle=90, colnames_offset_y=0, hjust=0, font.size=2)
	ecp  <- ecp  + scale_fill_manual("EC number", values=ec_colors)
} else {
	ecp <- gheatmap(ggt, plot_data[,"EC.",drop=F], offset = 1.1, width=0.1, color=NULL, colnames_position="top", colnames_angle=90, colnames_offset_y=0, hjust=0, font.size=2)
        ecp  <- ecp  + scale_fill_manual("EC number", values=ec_colors)
	}

deeplocp <- ecp + new_scale_fill()
deeplocp <- gheatmap(deeplocp, deeploc[,"Location",drop=F], offset = 2.1, width=0.1, color=NULL, colnames_position="top", colnames_angle=90, colnames_offset_y=0, hjust=0, font.size=2)
deeplocp <- deeplocp + scale_fill_manual("Location", values=deeploc_colors)


domainp <- deeplocp + new_scale_fill()
domainp <- gheatmap(domainp, taxonomy_plotting[,"taxonomy", drop=F], offset = 3.1, width=0.1, color=NULL, colnames_position="top", colnames_angle=90, colnames_offset_y=0, hjust=0, font.size=2)
domainp <- domainp + scale_fill_manual("taxonomy", values=tax_colors) + ggtitle(family) + theme(plot.title = element_text(hjust = 0.5, vjust=-55))


# function to reduce size of legend
addSmallLegend <- function(myPlot, pointSize = 1, textSize = 7, spaceLegend = 0.5) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}
dd <- addSmallLegend(domainp)



#subfp <- gheatmap(ggt, plot_data[,"Subf",drop=F], offset = 0.1, width=0.1, color=NULL, colnames_position="top", colnames_angle=90, colnames_offset_y=0, hjust=0, font.size=2)
#subfp2 <- subfp + scale_fill_manual("Subfamily", values=subf_colors)
#subfp2 <- addSmallLegend(subfp2)
#leg_subfp <- get_legend(subfp2)

#ecp2 <- gheatmap(ggt, plot_data[,"EC.",drop=F], offset = 0.1, width=0.1, color=NULL, colnames_position="top", colnames_angle=90, colnames_offset_y=0, hjust=0, font.size=2)
#ecp2 <- ecp2 + scale_fill_manual("EC number", values=ec_colors)
#ecp2 <- addSmallLegend(ecp2)
#leg_ecp <- get_legend(ecp2)

#domainp2 <- gheatmap(ggt, taxonomy_plotting[,"taxonomy", drop=F], offset = 0.1, width=0.1, color=NULL, colnames_position="top", colnames_angle=90, colnames_offset_y=0, hjust=0, font.size=2)
#domainp2 <- domainp2 + scale_fill_manual("taxonomy", values=tax_colors)
#domainp2 <- addSmallLegend(domainp2)
#leg_domainp2 <- get_legend(domainp2)

#domainp <- domainp + theme(legend.position="none")
#p <- plot_grid(domainp, leg_subfp, leg_ecp, leg_domainp2, ncol=4, rel_widths=c(.7, .1, .1, .1))
p <- dd + theme(plot.margin = margin(-2, 0, -2, -2, "cm"))

print("Writing plots to file")
outfile <- paste(out_prefix, "/", family, "_tree.pdf", sep="")
pdf(file=outfile, width=11.8, height=7.3)
p
dev.off()

print("Saving R env for potential later use")
rfile <- paste(out_prefix,"/",family,".rData", sep="")
save.image(file=rfile)

