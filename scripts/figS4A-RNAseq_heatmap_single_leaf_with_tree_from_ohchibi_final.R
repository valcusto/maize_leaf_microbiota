######load packages
library(ohchibi)
library(dplyr)
library(paletteer)
library(scales)
library(ggtree)
library(biomaRt)
library(egg)

######set seed
set.seed(130816)

######set work directory
setwd("C:/Users/Castrillo Lab/OneDrive - The University of Nottingham/Publication 2019/Paper preparation/manuscripts/scripts")

######helper functions
source('0-Clean_up_plots.R')
source('0_7-chibi.heatmap.R')

### Read the contrasts and dataset
Res_contrasts <- readRDS(file = "../cleandata/res_rnaseq_contrasts_syncom_do.RDS")
Dat_av  <- readRDS(file = "../cleandata/res_tab_av_rnaseq.RDS")
Tab_av <- Dat_av$Tab_av

#######select only the significant genes from L3 and NB
Res_sig_l3  <- Res_contrasts %>% 
  subset(padj < 0.1 )  %>%
  dplyr::filter(.data = .,Contrast %>% grepl(pattern = "_L3_FullSynCom")) %>%
  droplevels

sig_genes_l3 <- Res_sig_l3$gene_id %>% as.character %>%
  unique

### Map to entrez id 
Tab_sub_l3 <- Tab_av[,colnames(Tab_av) %>% grep(pattern = "^L3")]

Tab_sub_l3 <- match(sig_genes_l3,rownames(Tab_sub_l3)) %>%
  Tab_sub_l3[.,]

res_heatmap_l3 <- chibi.heatmap(Tab = Tab_sub_l3,
                                dist_method_rows = "euclidean",
                                hclust_method_rows = "ward.D",
                                range_fill_heatmap = c(-2,2),
                                k_rows = 11,panel_spacing = 0.2,
                                k_cols = 5,size_axis_text_row = 6,axis_ticks_row = T)

a <- res_heatmap_l3$heatmap
a

#######calculate the distance
dist_rows_l3 <- dist(Tab_sub_l3,method = 'euclidean')

######mclust
mclust_rows_l3 <- hclust(d = dist_rows_l3,method = 'ward.D')

######plot the mclust
mclust_rows_l3 %>% plot

######order according to mclust
order_rows_l3 <- mclust_rows_l3$order %>% mclust_rows_l3$labels[.]

######define the clusters
df_clust_rows_l3 <- mclust_rows_l3 %>% cutree(k = 11) %>%
  data.frame(IdRows = names(.), ClusterRows = paste0("CR",.),row.names = NULL)

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

######change the columns
melted_l3$IdCols <- melted_l3$IdCols %>%
  gsub(pattern = "^L[345]_",replacement = "") %>% 
  factor(levels = c("NB","FullSynCom","L3","L4_up","L5"))

######rename the columns
melted_l3$IdCols <- melted_l3$IdCols %>% as.character
melted_l3$IdCols[which(melted_l3$IdCols == 'FullSynCom')] <- '+SynCom'
melted_l3$IdCols[which(melted_l3$IdCols == 'L3')] <- '+SynCom-L3'
melted_l3$IdCols[which(melted_l3$IdCols == 'L4_up')] <- '+SynCom-L4'
melted_l3$IdCols[which(melted_l3$IdCols == 'L5')] <- '+SynCom-L5'

######order the factors
melted_l3$IdCols <- melted_l3$IdCols %>%
  factor(levels = c("NB","+SynCom","+SynCom-L3","+SynCom-L4","+SynCom-L5"))

######plot L3
p_l3 <- ggplot(data = melted_l3,mapping = aes(IdCols,IdRows)) +
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
        strip.text.y = element_text(size = 20, angle=360),
        legend.position = 'bottom',
        legend.title = element_text(size=15),
        legend.text = element_text(size=10),
        legend.key.width = unit(1, 'cm'))
p_l3

size_dendrogram = 0.3
######composition
tree_rows_l3 <- mclust_rows_l3 %>% as.phylo
p_tree_rows_l3 <- ggtree(tree_rows_l3,ladderize = F,size = size_dendrogram) +
  #geom_tiplab()+
  #The expand parth is fundamental to make the tree fit the composition with the other plots
  scale_y_continuous(expand =c(0.000001,0.1))


######compsotion
composition_l3 <- egg::ggarrange(
  p_tree_rows_l3,p_l3,
  ncol = 2,nrow = 1,
  widths = c(0.1,1))


######save as pdf
oh.save.pdf(p = composition_l3,
            outname = "figS4A-heatmap_L3_RNAseq_ohchibi_with_all_CRs.pdf",
            outdir = "../figures/",width = 12,height = 20)

#######select only the significant genes from L4 and NB
Res_sig_l4  <- Res_contrasts %>% 
  subset(padj < 0.1 )  %>%
  dplyr::filter(.data = .,Contrast %>% grepl(pattern = "_L4_FullSynCom")) %>%
  droplevels

sig_genes_l4 <- Res_sig_l4$gene_id %>% as.character %>%
  unique

### Map to entrez id 
Tab_sub_l4 <- Tab_av[,colnames(Tab_av) %>% grep(pattern = "^L4")]

Tab_sub_l4 <- match(sig_genes_l4,rownames(Tab_sub_l4)) %>%
  Tab_sub_l4[.,]

res_heatmap_l4 <- chibi.heatmap(Tab = Tab_sub_l4,
                                dist_method_rows = "euclidean",
                                hclust_method_rows = "ward.D",
                                range_fill_heatmap = c(-2,2),
                                k_rows = 8,panel_spacing = 0.2,
                                k_cols = 5,size_axis_text_row = 6,axis_ticks_row = T)

b <- res_heatmap_l4$heatmap
b

#######calculate the distance
dist_rows_l4 <- dist(Tab_sub_l4,method = 'euclidean')

######mclust
mclust_rows_l4 <- hclust(d = dist_rows_l4,method = 'ward.D')

######plot the mclust
mclust_rows_l4 %>% plot

######order according to mclust
order_rows_l4 <- mclust_rows_l4$order %>% mclust_rows_l4$labels[.]

######define the clusters
df_clust_rows_l4 <- mclust_rows_l4 %>% cutree(k = 8) %>%
  data.frame(IdRows = names(.), ClusterRows = paste0("CR",.),row.names = NULL)

######remove the first column
df_clust_rows_l4 <- df_clust_rows_l4[,-1]

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

######change the columns
melted_l4$IdCols <- melted_l4$IdCols %>%
  gsub(pattern = "^L[345]_",replacement = "") %>% 
  factor(levels = c("NB","FullSynCom","L3","L4_up","L5"))

######rename the columns
melted_l4$IdCols <- melted_l4$IdCols %>% as.character
melted_l4$IdCols[which(melted_l4$IdCols == 'FullSynCom')] <- '+SynCom'
melted_l4$IdCols[which(melted_l4$IdCols == 'L3')] <- '+SynCom-L3'
melted_l4$IdCols[which(melted_l4$IdCols == 'L4_up')] <- '+SynCom-L4'
melted_l4$IdCols[which(melted_l4$IdCols == 'L5')] <- '+SynCom-L5'

######order the factors
melted_l4$IdCols <- melted_l4$IdCols %>%
  factor(levels = c("NB","+SynCom","+SynCom-L3","+SynCom-L4","+SynCom-L5"))

######plot l4
p_l4 <- ggplot(data = melted_l4,mapping = aes(IdCols,IdRows)) +
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
        strip.text.y = element_text(size = 20, angle=360),
        legend.position = 'bottom',
        legend.title = element_text(size=15),
        legend.text = element_text(size=10),
        legend.key.width = unit(1, 'cm'))
p_l4

size_dendrogram = 0.3
######composition
tree_rows_l4 <- mclust_rows_l4 %>% as.phylo
p_tree_rows_l4 <- ggtree(tree_rows_l4,ladderize = F,size = size_dendrogram) +
  #geom_tiplab()+
  #The expand parth is fundamental to make the tree fit the composition with the other plots
  scale_y_continuous(expand =c(0.000001,0.1))


######compsotion
composition_l4 <- egg::ggarrange(
  p_tree_rows_l4,p_l4,
  ncol = 2,nrow = 1,
  widths = c(0.1,1))


######save as pdf
oh.save.pdf(p = composition_l4,
            outname = "figS4A-heatmap_l4_RNAseq_ohchibi_with_all_CRs.pdf",
            outdir = "../figures/",width = 12,height = 20)

