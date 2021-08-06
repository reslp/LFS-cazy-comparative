library(phytools)
library(ggplot2)
library(ggpubr)
library(ggtree)
library(cowplot)
library(reshape2)
library(patchwork)
library(tidyverse)
library(RColorBrewer)

options(StringsAsFactors=FALSE)

args <- commandArgs(trailingOnly=TRUE)

wd <- args[1]
treefile <- args[2]
countsfile <- args[3]

## function definitions. These functions have been modified from:
# https://thackl.github.io/ggtree-composite-plots
base_size <- 11
tree_y <-  function(ggtree, data){
  if(!inherits(ggtree, "ggtree"))
    stop("not a ggtree object")
  left_join(select(data, label), select(ggtree$data, label, y)) %>%
    pull(y)
}

# get the range of the ggtree y-axis data
tree_ylim <- function(ggtree){
  if(!inherits(ggtree, "ggtree"))
    stop("not a ggtree object")
  range(ggtree$data$y)
}

no_y_axis <- function ()
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
no_x_axis <- function ()
  theme(axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# plot data next to a ggtree aligned by shared labels
ggtreeplot <- function(ggtree, data = NULL, mapping = aes(), flip=FALSE,
                       expand_limits=expand_scale(0,.6), ...){
  #print(head(data))
  if(!inherits(ggtree, "ggtree"))
    stop("not a ggtree object")
  
  # match the tree limits
  limits <- tree_ylim(ggtree)
  limits[1] <- limits[1] + (limits[1] * expand_limits[1]) - expand_limits[2]
  limits[2] <- limits[2] + (limits[2] * expand_limits[3]) + expand_limits[4]

  if(flip){
    mapping <- modifyList(aes_(x=~x), mapping)
    data <- mutate(data, x=tree_y(ggtree, data))
    gg <- ggplot(data=data, mapping = mapping, ...) +
      scale_x_continuous(limits=limits, expand=c(0,0))
  }else{
    mapping <- modifyList(aes_(y=~y), mapping)
    data <- mutate(data, y=tree_y(ggtree, data))
    gg <- ggplot(data=data, mapping = mapping, ...) +
      scale_y_continuous(limits=limits, expand=c(0,0))
  }
  return(gg)
}
no_legend <- function() theme(legend.position="none")

scale_y_tree <- function(expand=expand_scale(0, 0.6), ...){
  scale_y_continuous(expand=expand, ...)
}

get_tree <- function() {
  ptree <- ggtree(tree, size=0.3) + geom_tiplab(align=TRUE, size=2.5) + scale_y_tree()+ xlim(0, 1500)
  
  
  #highlight nodes in tree for plotting:
  #bla <- c("R", "ELS", "ALSL", "ALSS", "OLS", "AO", "AX", "AL", "AE")
  #print(bla)
  #print(my_nodes)
  #ptree <- ptree + geom_point2(color="black",size=2) + geom_text2(nudge_x=-0.03,nudge_y=0.7, size=2)
  return(ptree)
}

create_heatmap <- function(df, heat_colors) {
  print("create heatmap plot")
  #print(head(df))
  ptree <- get_tree()
  
  p <- ggtreeplot(ptree, df, aes(x = category)) +
    geom_point(aes(x=category, color = value, size = value), alpha=.6) + 
    scale_color_gradientn(colors= heat_colors, na.value=NA) +
    geom_text(aes(label = ifelse(df$value > 0, df$value, "")), size = 2, color="grey50") +
    scale_size(range = c(1,max(df$value/3))) +theme_classic()+no_y_axis() +no_x_axis()+theme(axis.text.x = element_text(size = base_size *0.6, angle=45, hjust = 0, colour = "grey50")) + no_legend() +scale_x_discrete(position="top")
  
  
  #p <- ggtreeplot(ptree, df, aes(x=Orthogroup)) + geom_tile(aes(fill=value)) +scale_fill_gradient2(low = "#EAF4F7", high = "#F06449") +theme_classic()+no_y_axis() +no_x_axis()+theme(axis.text.x = element_text(size = base_size *0.8, angle=45, hjust = 0, colour = "grey50")) + no_legend() +scale_x_discrete(position="top")
  #p <- p + scale_y_continuous(expand=c(0,0)) + geom_text(aes(label = ifelse(df$value > 0, df$value, "")), colour="black", size=2)
  return(p)
}

create_heatmap2 <- function(df, heat_colors) {
  print("create heatmap plot")
  #print(head(df))
  #ptree <- get_tree()
  
  p <- ggplot(df, aes(x=category, y=label,width.x=100)) + geom_point(aes(x=category, color = value, size = value,y=label), alpha=.6) + 
    scale_color_gradientn(colors= heat_colors, na.value=NA) +
    geom_text(aes(label = ifelse(df$value > 0, df$value, "")), size = 2, color="grey50") +
    scale_size(range = c(1,max(df$value/3))) +theme_classic() +no_y_axis()+no_x_axis()+theme(axis.text.y = element_text(colour = "grey50", size = base_size*0.5, hjust=1), axis.text.x.top = element_text(size = base_size *0.6, angle=90, hjust = 0,vjust=0.5, colour = "grey50")) + no_legend() +
    scale_x_discrete(position="top")+scale_y_discrete()
  gt <- cowplot::as_gtable(p)
  #p <- ggtreeplot(ptree, df, aes(x=Orthogroup)) + geom_tile(aes(fill=value)) +scale_fill_gradient2(low = "#EAF4F7", high = "#F06449") +theme_classic()+no_y_axis() +no_x_axis()+theme(axis.text.x = element_text(size = base_size *0.8, angle=45, hjust = 0, colour = "grey50")) + no_legend() +scale_x_discrete(position="top")
  #p <- p + scale_y_continuous(expand=c(0,0)) + geom_text(aes(label = ifelse(df$value > 0, df$value, "")), colour="black", size=2)
  return(p)
}


Blues <- function(n) {
  cols<-colorRampPalette(brewer.pal(9,"Blues"))(n)
  cols[1:n]
}
# get colors

### end of function definitions

tree <- read.tree(treefile)
tree
data <- read.csv(countsfile, sep="\t")
data[is.na(data)] <- 0
data$Total <-NULL
data$Characterized_transporters <- NULL
rownames(data) <- make.names(data$Orthogroup, unique=TRUE)

data$Orthogroup <- NULL
data$Orthogroup_old <- NULL
data <- t(data)

data <- data[match(rev(tree$tip.label), rownames(data)), ]
colnames(data) <- gsub(";", "/", colnames(data))
colnames(data)

## correct tip labels in tree:
tree$tip.label <- gsub("_", " ", tree$tip.label)
rownames(data) <- gsub("_", " ", rownames(data))

#subsample to plot only a few:
#data[,c(2,3,4,11,15)]


data2 <- melt(data)
#logtransform for testing
#data2$value <-log(data2$value+1)
head(data2)
colnames(data2) <- c("label", "category", "value")

data2 <- as_tibble(data2)
heat_colors <- Blues(30)



create_heatmap(data2, heat_colors)

treep <- get_tree()
#treep + xlim(0, 1500)

#range(treep$data$x)

hm <- create_heatmap(data2, heat_colors)
hm2 <- create_heatmap2(data2, heat_colors)
hm+theme(axis.text.x = element_text(size = base_size *0.6, angle=90, hjust = 0, colour = "grey50"))
data[,c(2,3)]

pdf(file=paste0(wd,"/transporter_heatmap.pdf"), width=11.7, height=8.3)
print(hm2)
dev.off()

#pdf(file="transporter_heatmap_combined.pdf", width=11.7, height=8.3)
#print(treep+hm)
#dev.off()



