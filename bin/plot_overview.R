library(ape)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(wesanderson)
library(ggpubr)
library(patchwork)
options(show.error.locations = TRUE)

args <- commandArgs(trailingOnly=TRUE)
print(args)

wd <- args[1]
cogs_file <- args[2]
cazy_file <- args[3]
secmet_file <- args[4]
stats_file <- args[5]
outfile <- args[6]

#setwd("/Users/sinnafoch/Dropbox/Philipp/xylographa_comparative_genomics/tmp/new_overview")
setwd(wd)

# load all input files:
#cogs <- read.csv("COGS.all.results.csv")
cogs <- read.csv(cogs_file)
cogs <- melt(cogs)
#cazy <- read.csv("CAZyme.summary.results.csv")
cazy <- read.csv(cazy_file)
cazy <- melt(cazy)
#secmet <- read.csv("SM.summary.results.csv")
secmet <- read.csv(secmet_file)
secmet <- melt(secmet)
#stats <- read.table("genome.stats.summary.csv", sep=",", header=T)
stats <- read.table(stats_file, sep=",", header=T)


#get species names and set current species.
# this will have to be a for loop later
rownames(stats) <- stats$X
stats$X <- NULL
species_names <- colnames(stats)

all_plots <- list()
i <- 1
#generate plots
print(species_names)
for (species in species_names){
  print(species)
  genus <- substr(strsplit(species,"\\_")[[1]][1], 1,1 ) # get only first letter of genus name
  short_sp_name <- paste(genus, "..", strsplit(species,"\\_")[[1]][2], sep="")
  print(short_sp_name)
  stats_sp <- stats[species] # this needs to be adjusted because this file has different names
  sp <- short_sp_name
  sp_title <- str_replace(species, "\\.", " ")
  # get color palette for cogs:
  print("Get color palette: ")
  cogs_palette = colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(cogs$X)))
  
  
  # for here on, stuff will be created for every species:
  # calculate stuff for the pie chart
  print("Extracting COGS for species")
  current_sp_cogs <- cogs[cogs$variable==species, ]
  current_sp_cogs$fraction = current_sp_cogs$value / sum(current_sp_cogs$value)
  current_sp_cogs$ymax = cumsum(current_sp_cogs$fraction)
  current_sp_cogs$ymin = c(0, head(current_sp_cogs$ymax, n=-1))
  current_sp_cogs$label <- paste0(current_sp_cogs$X, "\n value: ", current_sp_cogs$count)
  current_sp_cogs$labelPosition <- (current_sp_cogs$ymax + current_sp_cogs$ymin) / 2
  
  # create plot for COGS
  print("Creating COGS plot")
  cog_plot <- ggplot(current_sp_cogs, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=X)) +
    geom_rect() +
    #geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
    scale_fill_manual(values=cogs_palette) +
    coord_polar(theta="y") +
    xlim(c(0.5, 4)) +
    theme_void() +
    theme(legend.position = "none") +ggtitle("COGs")
  
  # create plot for CAZymes
  print("Extract CAZY")
  cazy_palette <- colorRampPalette(wes_palette("Darjeeling1"))(6)
  current_sp_cazy <- cazy[cazy$variable==species,]
  cazy_plot <-ggplot(current_sp_cazy, aes(x=X, y=value, fill=as.factor(X))) + 
    geom_bar(stat = "identity") +theme_classic()+theme(axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank()) +
    scale_fill_manual(values=cazy_palette) +theme(legend.position = "none") + ggtitle("CAZymes")
  
  print("Extract SECMET")
  # create plot for secondary metabolite genes
  secmet_palette <- colorRampPalette(wes_palette("Cavalcanti1"))(length(unique(secmet$X)))
  
  current_sp_secmet <- secmet[secmet$variable==species,]
  secmet_plot <- ggplot(current_sp_secmet, aes(x=X, y=value, fill=as.factor(X))) + 
    geom_bar(stat = "identity") +theme_classic()+theme(axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank()) +
    scale_fill_manual(values=secmet_palette) +theme(legend.position = "none") + ggtitle("Seconday Metabolite genes")
  
  # create overview table:
  print("Overview table")
  all_rownames <- rownames(stats_sp)
  to_drop <- c("isolate", "locus_tag", "Prots atleast 1 ortholog", "Single-copy orthologs")
  stats_sp <- stats_sp %>% filter( ! rownames(stats_sp) %in% to_drop)
  new_rownames <- all_rownames[!all_rownames %in% to_drop]
  rownames(stats_sp) <- new_rownames
  
  
  tbody.style = tbody_style(color = "black", hjust=0, x=0.1, size=6) # this should left justify the text, but it does not work for some reason
  stats_plot <- ggtexttable(stats_sp, rows = rownames(stats_sp), 
                            theme = ttheme("blank", padding=unit(c(1,1),"mm"), base_size = 7, tbody.style = tbody.style))
  # create the assembly of plots:
  print("Plot assembly")
  sp_plot_assembly <- cog_plot + cazy_plot + secmet_plot + stats_plot + plot_layout(ncol=4) + plot_annotation(title = sp_title)
  sp_plot_assembly <- sp_plot_assembly & theme(text=element_text(size=6)) 
  all_plots[[i]] <- sp_plot_assembly
  i <- i + 1
}

how_many <- seq(1, length(all_plots), 7)
#save plots to file
pdf(width=8.3, height=11.7, file = outfile)
for (i in (1:length(how_many))){
 if (i == length(how_many)) { # special case for last element
   my_nums <- seq(how_many[i],length(all_plots))
   my_nums_length <- length(my_nums)
   if (my_nums_length != 7) { # if the last set of plots is < 7 add empty plots
     print("Will add empty plots")
     total_length <- length(all_plots)
     for (j in (1:(7-my_nums_length))) {
       all_plots[[total_length+j]] <- plot_spacer()
     }
     my_nums <- seq(how_many[i],length(all_plots))
     print(my_nums)
     print(all_plots[[my_nums[1]]] / all_plots[[my_nums[2]]] / all_plots[[my_nums[3]]] / all_plots[[my_nums[4]]] /all_plots[[my_nums[5]]] / all_plots[[my_nums[6]]] / all_plots[[my_nums[7]]])
     break
   }
   print(all_plots[[my_nums[1]]] / all_plots[[my_nums[2]]] / all_plots[[my_nums[3]]] / all_plots[[my_nums[4]]] /all_plots[[my_nums[5]]] / all_plots[[my_nums[6]]] / all_plots[[my_nums[7]]])
   break
 }
  #print(length(all_plots[how_many[i]:how_many[i+1]-1]))
  my_nums <- seq(how_many[i],how_many[i+1]-1)
  print(my_nums)
  print(all_plots[[my_nums[1]]] / all_plots[[my_nums[2]]] / all_plots[[my_nums[3]]] / all_plots[[my_nums[4]]] /all_plots[[my_nums[5]]] / all_plots[[my_nums[6]]] / all_plots[[my_nums[7]]])

}
dev.off()




