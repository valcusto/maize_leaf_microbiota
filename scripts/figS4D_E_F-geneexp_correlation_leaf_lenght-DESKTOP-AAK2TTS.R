########load packages
library(ohchibi)
library(paletteer)
library(scales)
library(car)
library(dplyr)
library(ggtree)
library(biomaRt)

######set seed
set.seed(130816)

######set work directory
setwd("C:/Users/Castrillo Lab/OneDrive - The University of Nottingham/Publication 2019/Paper preparation/manuscripts/scripts")

######helper functions
source('0_5-heatmap_ohchibi.R')
source('0-Clean_up_plots.R')
size_strip_text_row = 9
size_strip_text_col = 0
size_axis_text_col = 10
size_axis_text_row = 0
axis_ticks_row = F
size_legend_text = 10
size_legend_title = 11
size_border_tile = 0.5
width_border_tile = 0.85
height_border_tile = 0.85
palette_border = c("black")
legend_proportion_size = 0.75
size_dendrogram = 0.3
panel_border_heatmap = 0.01
panel_spacing = 0
font_family = "Helvetica"
font_face = "plain"
legend_position = "right"
tree_scale_factor_cols = c(0.005,0.005)
tree_scale_factor_rows = c(0.001,0.001)
egg_heights = c(0.2,1)
egg_widths = c(0.1,1)

### Read the contrasts and dataset
Res_contrasts <- readRDS(file = "../cleandata/res_rnaseq_contrasts_syncom_do.RDS")
Dat_av  <- readRDS(file = "../cleandata/res_tab_av_rnaseq.RDS")
Tab_av <- Dat_av$Tab_av
#df_trans <- Dat_av$df_trans
#df_dict_trans <- df_trans[,c("ensembl_gene_id","NCBI_gene_id")] %>% unique

#####################L3
### Correlation 
df_cor <- read.table(file = "../rawdata/res_cor_geneexp_length.tsv")
df_l3 <- df_cor %>% subset(V4 < 0.01) %>% subset(V2 == "L3")
up_l3 <- df_l3 %>% subset(V3 > 0)
down_l3 <- df_l3 %>% subset(V3 < 0)

######find the range
max(up_l3$V3)#0.70
min(up_l3$V3)#0.46

max(down_l3$V3)#-0.46
min(down_l3$V3)#-0.68

######select the significants for L3
sig_genes_l3 <- df_l3$V1 %>% as.character %>%
  unique

### Map to entrez id 
Tab_sub_l3 <- Tab_av[,colnames(Tab_av) %>% grep(pattern = "^L3")]

Tab_sub_l3 <- match(sig_genes_l3,rownames(Tab_sub_l3)) %>%
  Tab_sub_l3[.,]

######rename the Tab_sub_l3 column names
colnames(Tab_sub_l3) <- c('+SynCom', '+SynCom-L3', '+SynCom-L4', 
                          '+SynCom-L5', 'NB')

#######calculate the distance
dist_rows_l3 <- dist(Tab_sub_l3,method = 'euclidean')

######mclust
mclust_rows_l3 <- hclust(d = dist_rows_l3,method = 'ward.D')

######plot the mclust
mclust_rows_l3 %>% plot

######order according to mclust
order_rows_l3 <- mclust_rows_l3$order %>% mclust_rows_l3$labels[.]

######define the clusters
df_clust_rows_l3 <- mclust_rows_l3 %>% cutree(k = 9) %>%
  data.frame(IdRows = names(.), ClusterRows = paste0("CR",.,"C"),row.names = NULL)

######remove the first column
df_clust_rows_l3 <- df_clust_rows_l3[,-1]

######melt the tab
melted_l3 <- Tab_sub_l3 %>% melt
colnames(melted_l3) <- c("IdRows","IdCols","value")

######merge the two datasets
melted_l3 <- merge(melted_l3,df_clust_rows_l3, by = "IdRows")

######order the rows
melted_l3$IdRows <- melted_l3$IdRows %>%
  factor(levels = order_rows_l3)

######order groups
order_groups_rows_l3 <- with(melted_l3,order(IdRows)) %>%
  melted_l3$ClusterRows[.] %>% as.character %>%
  unique

######add it to the ClusterRows
melted_l3$ClusterRows <- melted_l3$ClusterRows %>%
  factor(levels = order_groups_rows_l3 %>% rev)

######helper functions
palette_heatmap = "pals::kovesi.diverging_bwr_55_98_c37"
range_fill_heatmap = c(-2,2)

######plot L3
p_l3 <- ggplot(data = melted_l3,mapping = aes(IdCols,IdRows)) +
  geom_tile(aes(fill = value),color = '#000000') +
  facet_grid(ClusterRows~.,scales = "free",space = "free") +
  scale_x_discrete(expand = c(0,0), position = 'top')+
  scale_y_discrete(expand =c(0,0)) +
  ggtitle('L3')+
  scale_fill_paletteer_c(palette_heatmap,
                         limits = range_fill_heatmap,
                         oob = squish,name = "Standardized\nexpression") +
  clean +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title.y =element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=20, 
                                   angle=360,
                                   hjust = 0.5),
        panel.background = element_blank(),
        panel.spacing = unit (0.5, 'lines'),
        strip.background = element_blank(),
        strip.text = element_text(size = 30),
        strip.text.y = element_text(size = 15),
        legend.text = element_text(size=10),
        legend.key.width = unit(1, 'cm'))

p_l3

######select only NB and FullSynCom
melted_sub_l3 <- melted_l3 %>% subset((IdCols == 'NB') | 
                                        (IdCols == '+SynCom')) %>% droplevels

######order the IdCols
melted_sub_l3$IdCols <- melted_sub_l3$IdCols %>%
  factor(levels = c('NB', '+SynCom'))

######plot
p_sub_l3 <- ggplot(data = melted_sub_l3,mapping = aes(IdCols,IdRows)) +
  geom_tile(aes(fill = value)) +
  facet_grid(ClusterRows~.,scales = "free",space = "free") +
  scale_x_discrete(expand = c(0,0), position = 'top')+
  scale_y_discrete(expand =c(0,0)) +
  ggtitle('L3')+
  scale_fill_paletteer_c(palette_heatmap,
                         limits = range_fill_heatmap,
                         oob = squish,name = "Standardized\nexpression") +
  clean +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.y =element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 20),
        panel.background = element_blank(),
        panel.spacing = unit (0.2, 'lines'),
        panel.border = element_rect(color="black", size=0.1, linetype="solid"),
        strip.background = element_blank(),
        strip.text = element_text(size = 30),
        strip.text.y = element_blank(),
        legend.position = 'bottom',
        legend.title = element_text(size=15),
        legend.text = element_text(size=10),
        legend.key.width = unit(1, 'cm'))

p_sub_l3

######composition
tree_rows_l3 <- mclust_rows_l3 %>% as.phylo
p_tree_rows_l3 <- ggtree(tree_rows_l3,ladderize = F,size = size_dendrogram) +
  #geom_tiplab()+
  #The expand parth is fundamental to make the tree fit the composition with the other plots
  scale_y_continuous(expand =c(0.000001,0.1))

######order the correlation
df_l3 <- with(df_l3,order(-V3)) %>%
  df_l3[.,]
df_l3$V1 <- df_l3$V1 %>% factor(levels = order_rows_l3)

######change the first column name on df_l3
colnames(df_l3)[1] <- 'IdRows'

######add the clusterRows information
df_l3 <- merge(df_l3, df_clust_rows_l3, by='IdRows')

#######order the clusterRows
df_l3$ClusterRows <- df_l3$ClusterRows %>%
  factor(levels = order_groups_rows_l3 %>% rev)

######create the label
m_text <- 'Spearman p against leaf length'

#######plot the correlation with leaf length
p_l3 <- ggplot(data = df_l3,aes(IdRows,V3)) +
  geom_bar(stat = "identity", color='grey') +
  geom_hline(yintercept = 0.46)+
  geom_hline(yintercept = -0.46)+
  facet_grid(ClusterRows~.,scales = "free",space = "free")+
  coord_flip() +
  ylab(label = m_text)+
  ylim(-1,1)+
  clean +
  theme(axis.title.y =element_blank(),
        axis.title.x = element_text(size=20),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit (0.2, 'lines'),
        panel.border = element_rect(color="black", size=0.1, linetype="solid"),
        strip.background = element_blank(),
        strip.text = element_text(size = 20),
        strip.text.y = element_text(size = 15, angle = 360))

p_l3

######compsotion
composition_l3 <- egg::ggarrange(
  p_tree_rows_l3,p_sub_l3, p_l3,
  ncol = 3,nrow = 1,
  widths = c(0.1,1,1))

######save as pdf
oh.save.pdf(p = composition_l3,
            outname = "figS4D-heatmap_L3_RNAseq_correlation_leaf_with_CRs.pdf",
            outdir = "../figures/",width = 16,height = 12)

######check in which cluster the selection of the genes belongs
df_negative <- read.csv('../rawdata/negative_selection.csv')

#####merge the df_negative with df_clust_rows_l3
merged_selection <- merge(df_negative, df_clust_rows_l3, by='IdRows')

#####################L4
### Correlation 
df_l4 <- df_cor %>% subset(V4 < 0.01) %>% subset(V2 == "L4")
up_l4 <- df_l4 %>% subset(V3 > 0)
down_l4 <- df_l4 %>% subset(V3 < 0) #%$% V1

######find the range
max(up_l4$V3)#0.60
min(up_l4$V3)#0.46

max(down_l4$V3)#-0.46
min(down_l4$V3)#-0.63

######select the significants for l4
sig_genes_l4 <- df_l4$V1 %>% as.character %>%
  unique

### Map to entrez id 
Tab_sub_l4 <- Tab_av[,colnames(Tab_av) %>% grep(pattern = "^L4")]

Tab_sub_l4 <- match(sig_genes_l4,rownames(Tab_sub_l4)) %>%
  Tab_sub_l4[.,]

######rename the Tab_sub_l4 column names
colnames(Tab_sub_l4) <- c('+SynCom', '+SynCom-L3', '+SynCom-L4', 
                          '+SynCom-L5', 'NB')

#######calculate the distance
dist_rows_l4 <- dist(Tab_sub_l4,method = 'euclidean')

######mclust
mclust_rows_l4 <- hclust(d = dist_rows_l4,method = 'ward.D')

######plot the mclust
mclust_rows_l4 %>% plot

######order according to mclust
order_rows_l4 <- mclust_rows_l4$order %>% mclust_rows_l4$labels[.]

######define the clusters
df_clust_rows_l4 <- mclust_rows_l4 %>% cutree(k = 7) %>%
  data.frame(IdRows = names(.), ClusterRows = paste0("CR",.,"C"),row.names = NULL)

######remove the first column
df_clust_rows_l4 <- df_clust_rows_l4[,-1]

######trensform CR8 in CR7
df_clust_rows_l4$ClusterRows <- df_clust_rows_l4$ClusterRows %>% 
  as.character %>%
  gsub(pattern = "CR7C|CR4C",replacement = "CR4C")

######melt the tab
melted_l4 <- Tab_sub_l4 %>% melt
colnames(melted_l4) <- c("IdRows","IdCols","value")

######merge the two datasets
melted_l4 <- merge(melted_l4,df_clust_rows_l4, by = "IdRows")

######order the rows
melted_l4$IdRows <- melted_l4$IdRows %>%
  factor(levels = order_rows_l4)

######order groups
order_groups_rows_l4 <- with(melted_l4,order(IdRows)) %>%
  melted_l4$ClusterRows[.] %>% as.character %>%
  unique

######add it to the ClusterRows
melted_l4$ClusterRows <- melted_l4$ClusterRows %>%
  factor(levels = order_groups_rows_l4 %>% rev)

######helper functions
palette_heatmap = "pals::kovesi.diverging_bwr_55_98_c37"
range_fill_heatmap = c(-2,2)

######plot l4
p_sub_l4 <- ggplot(data = melted_l4,mapping = aes(IdCols,IdRows)) +
  geom_tile(aes(fill = value)) +
  facet_grid(ClusterRows~.,scales = "free",space = "free") +
  scale_x_discrete(expand = c(0,0), position = 'top')+
  scale_y_discrete(expand =c(0,0)) +
  ggtitle('L4')+
  scale_fill_paletteer_c(palette_heatmap,
                         limits = range_fill_heatmap,
                         oob = squish,name = "Standardized\nexpression") +
  clean +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.y =element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 20),
        panel.background = element_blank(),
        panel.spacing = unit (0.2, 'lines'),
        panel.border = element_rect(color="black", size=0.1, linetype="solid"),
        strip.background = element_blank(),
        strip.text = element_text(size = 30),
        strip.text.y = element_blank(),
        legend.position = 'bottom',
        legend.title = element_text(size=15),
        legend.text = element_text(size=10),
        legend.key.width = unit(1, 'cm'))

p_sub_l4

######select only NB and FullSynCom
melted_sub_l4 <- melted_l4 %>% subset((IdCols == 'NB') | 
                                        (IdCols == '+SynCom')) %>% droplevels

######order the IdCols
melted_sub_l4$IdCols <- melted_sub_l4$IdCols %>%
  factor(levels = c('NB', '+SynCom'))

######plot
p_sub_l4 <- ggplot(data = melted_sub_l4,mapping = aes(IdCols,IdRows)) +
  geom_tile(aes(fill = value)) +
  facet_grid(ClusterRows~.,scales = "free",space = "free") +
  scale_x_discrete(expand = c(0,0), position = 'top')+
  scale_y_discrete(expand =c(0,0)) +
  ggtitle('L4')+
  scale_fill_paletteer_c(palette_heatmap,
                         limits = range_fill_heatmap,
                         oob = squish,name = "Standardized\nexpression") +
  clean +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.y =element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 20),
        panel.background = element_blank(),
        panel.spacing = unit (0.2, 'lines'),
        panel.border = element_rect(color="black", size=0.1, linetype="solid"),
        strip.background = element_blank(),
        strip.text = element_text(size = 30),
        strip.text.y = element_blank(),
        legend.position = 'bottom',
        legend.title = element_text(size=15),
        legend.text = element_text(size=10),
        legend.key.width = unit(1, 'cm'))

p_sub_l4

######composition
tree_rows_l4 <- mclust_rows_l4 %>% as.phylo
p_tree_rows_l4 <- ggtree(tree_rows_l4,ladderize = F,size = size_dendrogram) +
  #geom_tiplab()+
  #The expand parth is fundamental to make the tree fit the composition with the other plots
  scale_y_continuous(expand =c(0.001,0.001))

######order the correlation
df_l4 <- with(df_l4,order(-V3)) %>%
  df_l4[.,]
df_l4$V1 <- df_l4$V1 %>% factor(levels = order_rows_l4)

######change the first column name on df_l3
colnames(df_l4)[1] <- 'IdRows'

######add the clusterRows information
df_l4 <- merge(df_l4, df_clust_rows_l4, by='IdRows')

#######order the clusterRows
df_l4$ClusterRows <- df_l4$ClusterRows %>%
  factor(levels = order_groups_rows_l4 %>% rev)

######create the label
m_text <- 'Spearman ?? against leaf length'

#######plot the correlation with leaf length
p_l4 <- ggplot(data = df_l4,aes(IdRows,V3)) +
  geom_bar(stat = "identity", color='grey') +
  geom_hline(yintercept = 0.46)+
  geom_hline(yintercept = -0.46)+
  facet_grid(ClusterRows~.,scales = "free",space = "free")+
  coord_flip() +
  ylab(label = m_text)+
  ggtitle('L4')+
  ylim(-1,1)+
  clean +
  theme(axis.title.y =element_blank(),
        axis.title.x = element_text(size=20),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit (0.2, 'lines'),
        panel.border = element_rect(color="black", size=0.1, linetype="solid"),
        strip.background = element_blank(),
        strip.text = element_text(size = 20),
        strip.text.y = element_text(size = 15, angle = 360))

p_l4

######compsotion
composition_l4 <- egg::ggarrange(
  p_tree_rows_l4,p_sub_l4, p_l4,
  ncol = 3,nrow = 1,
  widths = c(0.1,1,1))

######save as pdf
oh.save.pdf(p = composition_l4,
            outname = "figS4D-heatmap_L4_RNAseq_correlation_leaf_with_CRs.pdf",
            outdir = "../figures/",width = 16,height = 12)

#####################L5
### Correlation 
df_l5 <- df_cor %>% subset(V4 < 0.01) %>% subset(V2 == "L5")
up_l5 <- df_l5 %>% subset(V3 > 0) #%$% V1
down_l5 <- df_l5 %>% subset(V3 < 0) #%$% V1

######find the range
max(up_l5$V3)#0.69
min(up_l5$V3)#0.46

max(down_l5$V3)#-0.46
min(down_l5$V3)#-0.77

######select the significants for l5
sig_genes_l5 <- df_l5$V1 %>% as.character %>%
  unique

### Map to entrez id 
Tab_sub_l5 <- Tab_av[,colnames(Tab_av) %>% grep(pattern = "^L5")]

Tab_sub_l5 <- match(sig_genes_l5,rownames(Tab_sub_l5)) %>%
  Tab_sub_l5[.,]

######rename the Tab_sub_l5 column names
colnames(Tab_sub_l5) <- c('+SynCom', '+SynCom-L3', '+SynCom-L4', 
                          '+SynCom-L5', 'NB')

#######calculate the distance
dist_rows_l5 <- dist(Tab_sub_l5,method = 'euclidean')

######mclust
mclust_rows_l5 <- hclust(d = dist_rows_l5,method = 'ward.D')

######plot the mclust
mclust_rows_l5 %>% plot

######order according to mclust
order_rows_l5 <- mclust_rows_l5$order %>% mclust_rows_l5$labels[.]

######define the clusters
df_clust_rows_l5 <- mclust_rows_l5 %>% cutree(k = 8) %>%
  data.frame(IdRows = names(.), ClusterRows = paste0("CR",.,"C"),row.names = NULL)

######remove the first column
df_clust_rows_l5 <- df_clust_rows_l5[,-1]

######trensform CR8 in CR7
df_clust_rows_l5$ClusterRows <- df_clust_rows_l5$ClusterRows %>% 
  as.character %>%
  gsub(pattern = "CR7C|CR8C",replacement = "CR7C")

######melt the tab
melted_l5 <- Tab_sub_l5 %>% melt
colnames(melted_l5) <- c("IdRows","IdCols","value")

######merge the two datasets
melted_l5 <- merge(melted_l5,df_clust_rows_l5, by = "IdRows")

######order the rows
melted_l5$IdRows <- melted_l5$IdRows %>%
  factor(levels = order_rows_l5)

######order groups
order_groups_rows_l5 <- with(melted_l5,order(IdRows)) %>%
  melted_l5$ClusterRows[.] %>% as.character %>%
  unique

######add it to the ClusterRows
melted_l5$ClusterRows <- melted_l5$ClusterRows %>%
  factor(levels = order_groups_rows_l5 %>% rev)

######helper functions
palette_heatmap = "pals::kovesi.diverging_bwr_55_98_c37"
range_fill_heatmap = c(-2,2)

######plot l5
p_l5 <- ggplot(data = melted_l5,mapping = aes(IdCols,IdRows)) +
  geom_tile(aes(fill = value),color = '#000000') +
  facet_grid(ClusterRows~IdCols,scales = "free",space = "free") +
  scale_x_discrete(expand = c(0,0), position = 'top')+
  scale_y_discrete(expand =c(0,0)) +
  ggtitle('l5')+
  scale_fill_paletteer_c(palette_heatmap,
                         limits = range_fill_heatmap,
                         oob = squish,name = "Standardized\nexpression") +
  clean +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title.y =element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=20, 
                                   angle=360,
                                   hjust = 0.5),
        panel.background = element_blank(),
        panel.spacing = unit (0.5, 'lines'),
        strip.background = element_blank(),
        strip.text = element_text(size = 30),
        strip.text.y = element_text(size = 15),
        legend.text = element_text(size=10),
        legend.key.width = unit(1, 'cm'))

p_l5

######select only NB and FullSynCom
melted_sub_l5 <- melted_l5 %>% subset((IdCols == 'NB') | 
                                        (IdCols == '+SynCom')) %>% droplevels

######order the IdCols
melted_sub_l5$IdCols <- melted_sub_l5$IdCols %>%
  factor(levels = c('NB', '+SynCom'))

######plot
p_sub_l5 <- ggplot(data = melted_sub_l5,mapping = aes(IdCols,IdRows)) +
  geom_tile(aes(fill = value)) +
  facet_grid(ClusterRows~.,scales = "free",space = "free") +
  scale_x_discrete(expand = c(0,0), position = 'top')+
  scale_y_discrete(expand =c(0,0)) +
  ggtitle('L5')+
  scale_fill_paletteer_c(palette_heatmap,
                         limits = range_fill_heatmap,
                         oob = squish,name = "Standardized\nexpression") +
  clean +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.y =element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 20),
        panel.background = element_blank(),
        panel.spacing = unit (0.2, 'lines'),
        panel.border = element_rect(color="black", size=0.1, linetype="solid"),
        strip.background = element_blank(),
        strip.text = element_text(size = 30),
        strip.text.y = element_blank(),
        legend.position = 'bottom',
        legend.title = element_text(size=15),
        legend.text = element_text(size=10),
        legend.key.width = unit(1, 'cm'))

p_sub_l5

######composition
tree_rows_l5 <- mclust_rows_l5 %>% as.phylo
p_tree_rows_l5 <- ggtree(tree_rows_l5,ladderize = F,size = size_dendrogram) +
  #geom_tiplab()+
  #The expand parth is fundamental to make the tree fit the composition with the other plots
  scale_y_continuous(expand =c(0.001,0.001))

######order the correlation
df_l5 <- with(df_l5,order(-V3)) %>%
  df_l5[.,]
df_l5$V1 <- df_l5$V1 %>% factor(levels = order_rows_l5)

######change the first column name on df_l3
colnames(df_l5)[1] <- 'IdRows'

######add the clusterRows information
df_l5 <- merge(df_l5, df_clust_rows_l5, by='IdRows')

#######order the clusterRows
df_l5$ClusterRows <- df_l5$ClusterRows %>%
  factor(levels = order_groups_rows_l5 %>% rev)

######create the label
m_text <- 'Spearman ?? against leaf length'

#######plot the correlation with leaf length
p_l5 <- ggplot(data = df_l5,aes(IdRows,V3)) +
  geom_bar(stat = "identity", color='grey') +
  geom_hline(yintercept = 0.46)+
  geom_hline(yintercept = -0.46)+
  facet_grid(ClusterRows~.,scales = "free",space = "free")+
  coord_flip() +
  ylab(label = m_text)+
  ylim(-1,1)+
  #theme_ohchibi() +
  ggtitle('L5')+
  clean +
  theme(axis.title.y =element_blank(),
        axis.title.x = element_text(size=20),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit (0.2, 'lines'),
        panel.border = element_rect(color="black", size=0.1, linetype="solid"),
        strip.background = element_blank(),
        strip.text = element_text(size = 20),
        strip.text.y = element_text(size = 15, angle = 360))

p_l5

######compsotion
composition_l5 <- egg::ggarrange(
  p_tree_rows_l5,p_sub_l5, p_l5,
  ncol = 3,nrow = 1,
  widths = c(0.1,1,1))

######save as pdf
oh.save.pdf(p = composition_l5,
            outname = "figS4D-heatmap_L5_RNAseq_correlation_leaf_with_CRs.pdf",
            outdir = "../figures/",width = 16,height = 12)
