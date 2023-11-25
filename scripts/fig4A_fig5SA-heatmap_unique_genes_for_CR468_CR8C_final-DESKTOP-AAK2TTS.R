# Setting up environment ===================================================
# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
#options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

######load package
library(ohchibi)
library(dplyr)
library(scales)
library(paletteer)
library(ggtree)
library(tidyverse)

######set seed
set.seed(130816)

######set work directory
setwd("C:/Users/Castrillo Lab/OneDrive - The University of Nottingham/Publication 2019/Paper preparation/manuscripts/scripts")

######helper functions
source('0-Clean_up_plots.R')
source('0_5-heatmap_ohchibi.R')
source('0_6-theme_ohchibi.R')
source('0_7-chibi.heatmap.R')

######Read the contrasts and dataset
Dat_av  <- readRDS(file = "../cleandata/res_tab_av_rnaseq.RDS")
Tab_av <- Dat_av$Tab_av

######open the set of genes
df_gene <- read.csv('../cleandata/gene_setCR468_8c.csv')

#####select the unique genes
sig_unique <- df_gene$IdRows %>%
  unique

######select the genes on the Tab dataset
Tab_sub <- Tab_av[,colnames(Tab_av) %>% grep(pattern = "^L3")]

Tab_sub <- match(sig_unique,rownames(Tab_sub)) %>%
  Tab_sub[.,]

######create the clustering for the genes
mclust <- hclust(d = (1-cor(Tab_sub %>%t)) %>% as.dist,method = "ward.D2")
order <- mclust$order %>% mclust$labels[.]

######Visualise the dendogram
mclust %>% plot

######melt the dataset
melted_sub <- Tab_sub %>% melt %>%
  dplyr::rename(.data = .,gene_id = Var1,group = Var2)

###Order genes according to the cluster dendogram
melted_sub$gene_id <- melted_sub$gene_id %>% factor(levels = order)

#####Divide the dendogram in cluster
df_clustgene <- mclust %>% dendextend::cutree(tree=., k=4) %>%
  data.frame(gene_id=names(.), Clustgene=paste0("CR", .),
             row.names = NULL) 

#remove first column
df_clustgene <- df_clustgene[,-1]

######create a dataframe for the Tab_sub
df_new <- Tab_sub %>% data.frame

######create the IdRows column
df_new$gene_id <- rownames(df_new)

#####remove the rowname
rownames(df_new) <- NULL

######merge with the df_clustgene
df_new_final <- merge(df_clustgene, df_new, by='gene_id')

######save the csv
write.table(df_new_final, file='../cleandata/TableS7_1.csv', sep=',',
            row.names = F, col.names = T)

######prepare the dataset
melted_sub$group <- melted_sub$group %>% 
  factor(levels = c("L3_NB","L3_FullSynCom","L3_L3","L3_L4_up","L3_L5"))

######create the NumLeaf column
melted_sub$NumLeaf <- "L3"

######create the Treatment column
melted_sub$Treatment <- melted_sub$group %>%
  gsub(pattern = "^L[345]_",replacement = "") %>% as.character

######rename the information on the tReatment column
melted_sub$Treatment[which(melted_sub$Treatment == 'FullSynCom')] <- '+SynCom'
melted_sub$Treatment[which(melted_sub$Treatment == 'L3')] <- '+SynCom-L3'
melted_sub$Treatment[which(melted_sub$Treatment == 'L4_up')] <- '+SynCom-L4'
melted_sub$Treatment[which(melted_sub$Treatment == 'L5')] <- '+SynCom-L5'

######merge the melted and df_clustgene
melted_sum <- merge(melted_sub, df_clustgene, by='gene_id')

######order the Treatment columns
melted_sum$Treatment <- melted_sum$Treatment %>% 
  factor(levels=c('NB', '+SynCom', '+SynCom-L3','+SynCom-L4','+SynCom-L5'))

######plot the heatmap
p_heatmap <- ggplot(data =melted_sum,aes(x = Treatment,y = gene_id)) +
  geom_raster(aes(fill = value)) +
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-1.0,1.0),oob = squish)+
  facet_grid(Clustgene~Treatment,space = "free",scales = "free")+
  ylab(label = NULL)+
  xlab(label = NULL)+
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  labs(fill = "Standardized\nexpression")+
  clean+
  guides(color = guide_legend(override.aes = list(size = 10)))+
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title.y =element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text= element_text(size = 18),
        strip.text.y = element_text(size=20, angle = 360),
        #strip.text.y = element_blank(),
        legend.text = element_text(size=10),
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.position = 'bottom')

p_heatmap

######save the p_heatmap
oh.save.pdf(p = p_heatmap,
            outname = "fig4B-heatmap_L3_RNAseq_unique_four_clusters_gene_legend.pdf",
            outdir = "../figures/",width = 10,height = 20)

######list of genes selected in arabidopsis mutants
df_mutants <- read.csv(file = '../rawdata/arabidopsis_mutants.csv')

######mutants id
unique_mutants <- df_mutants$gene_id %>% unique

######add the column mutants in the object
melted_sum$mutants <- 'No'
melted_sum$mutants[which(melted_sum$gene_id %in% unique_mutants)] <- 'Yes'

######add color 
melted_sum$color <- 'white'
melted_sum$color[which(melted_sum$mutants == 'Yes')] <- 'black'

######define mutants as factor
melted_sum$mutants <- melted_sum$mutants %>%
  factor(levels = c('No', 'Yes'))

######add a column for the x
melted_sum$bar <- 'bar'

#######select only the yes mutants
df_yes <- melted_sum %>%
  subset(mutants == 'yes')
unique_yes <- df_yes$gene_id %>% unique

######plot the heatmap
p_heatmap1 <- ggplot(data =melted_sum,aes(x = Treatment,y = gene_id)) +
  geom_raster(aes(fill = value)) +
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-1.0,1.0),oob = squish)+
  facet_grid(Clustgene~Treatment,space = "free",scales = "free")+
  ylab(label = NULL)+
  xlab(label = NULL)+
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  labs(fill = "Standardized\nexpression")+
  clean+
  guides(color = guide_legend(override.aes = list(size = 10)))+
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title.y =element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text= element_text(size = 15),
        strip.text.y = element_text(size=20, angle = 360),
        #strip.text.y = element_blank(),
        legend.text = element_text(size=10),
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.position = 'bottom')

p_heatmap1

######mutants
p_heatmap2 <- ggplot(data =melted_sum,aes(x = bar,y = gene_id)) +
  geom_raster(aes(fill = mutants)) +
  scale_fill_manual(name = 'Arabidopsis\nmutants',values = c('white', 'black'))+
  facet_grid(Clustgene~.,space = "free",scales = "free")+
  clean+
  guides(color = guide_legend(override.aes = list(size = 10)))+
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title.y =element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text= element_text(size = 10),
        #strip.text.y = element_text(size=20, angle = 360),
        strip.text.y = element_blank(),
        legend.text = element_text(size=15),
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.position = 'left')

p_heatmap2

######tree composition
composition <- egg::ggarrange(p_heatmap2,p_heatmap1,
                              nrow=1, ncol=2, widths = c( 0.02,1))
######save the p_heatmap
oh.save.pdf(p = composition,
            outname = "fig5A-heatmap_L3_RNAseq_unique_gene_four_clusters_legend_selected_finalTEST.pdf",
 