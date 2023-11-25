######load packages
library(ohchibi)
library(dplyr)
library(paletteer)
library(reshape2)
library(tidyverse)
library(egg)
library(ggrepel)
library(ggtree)

######set seed
set.seed(12356)

######set work directory
setwd("C:/Users/Castrillo Lab/OneDrive - The University of Nottingham/Publication 2019/Paper preparation/manuscripts/scripts")

######helper functions
source('0-Clean_up_plots.R')
#source('0_2-create_dataset_object.R')
#source('0_3-save_pdf.R')
#source('0_4-subset_Dataset.R')

######Open Dat object
Dat <- readRDS(file = "../cleandata/dat_leafionome_with_without_syncom.RDS")

######subset only the sterile condition
Dat_sub <- Dat %>%
  subset(Treatment != 'Full syncom', drop = T)

#####Exploring the data
Tab <- Dat_sub$Tab  %>% t

###Transform in z-score matrix 
#A z-score is a measure of how many standard deviations 
#(how far my observation is from the mean) below or above
#the populatio mean a raw score is.
Tab_z <- Tab %>% scale

#####melt the matrix
melted <- Tab_z %>% melt

#####name the columns
colnames(melted) <- c('UiD', 'ion', 'zscore')

######merge the melted data with the metadata
melted <- merge(melted, Dat_sub$Map, by='UiD')

######create a group column
melted$group <- paste(melted$Nutrient, melted$NumLeaf, 
                      sep = '_')

######average the Ions
Tab_sum <- acast(data = melted,formula = ion~group,
                 fun.aggregate = mean,value.var = "zscore")

######cluster the ions
mclust_ions <- hclust(d = as.dist(1-cor(Tab_sum %>% t)),
                      method = "complete")

######plot cluster
plot(mclust_ions)

######order ions
order_ions <- mclust_ions$order %>% mclust_ions$labels[.]

######melted the tab_sum
melted_sum <- Tab_sum %>% melt
colnames(melted_sum) <- c("ion","group","zscore")

######separate the different group in the melted_sum
melted_sum$Nutrient <- sapply(strsplit(as.character(melted_sum$group),
                                   '_'), `[`, 1) 

melted_sum$NumLeaf <- sapply(strsplit(as.character(melted_sum$group),
                                   '_'), `[`, 2)

#To plot with ggplot you wanan order the factor 
melted_sum$ion <- melted_sum$ion %>% factor(levels = order_ions)

######organize the nutrient 
melted_sum$Nutrient <- melted_sum$Nutrient %>%
  factor(levels = c('Full Hoagland', '1/4 Hoagland'))

######plot 
p_heatmap <- ggplot(data =melted_sum,aes(x = Nutrient,y = ion)) +
  geom_raster(aes(fill = zscore)) +
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-1.5,1.5),oob = squish)+
  facet_grid(.~NumLeaf)+
  ylab(label = NULL)+
  xlab(label = NULL)+
  clean+
  guides(color = guide_legend(override.aes = list(size = 10)))+
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title.y =element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 15, angle = 45, hjust=1, vjust=1),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        legend.text = element_text(size=10),
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.position = 'bottom')

p_heatmap

######create a tree for the ion clustering
tree_ion <- mclust_ions %>% as.phylo

######design the tree
p_tree_ion <- ggtree(tree_ion, ladderize = F, branch.length='none',size=1.2)+
  #The expand parth is fundamental to make the tree fit the composition with the other plots
  scale_y_continuous(expand = c(0.05,0.05))

####composition
composition <- egg::ggarrange(
  #Row 1
  p_tree_ion, p_heatmap,
  #Hoe many rows and columns
  nrow = 1,ncol = 2,byrow = T,
  #Here you define the width of the rows and columns
  #We wanna make the heatmap enormous in comparison to the dendrograms
  widths = c(0.2,1))

######save the plot as pdf
oh.save.pdf(p = composition,outname = "figS3G_Heatmap_leaf_ionome_sterile_condition.pdf",
            outdir = "../figures/",width = 10,height = 8)
