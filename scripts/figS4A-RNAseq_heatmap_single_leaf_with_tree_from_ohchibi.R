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
#source('0_7-chibi.heatmap.R')
source('0_5-heatmap_ohchibi.R')

### Read the contrasts and dataset
Res_contrasts <- readRDS(file = "../cleandata/res_rnaseq_contrasts_syncom_do.RDS")
Dat_av  <- readRDS(file = "../cleandata/res_tab_av_rnaseq.RDS")
Tab_av <- Dat_av$Tab_av

######select all the significant genes from the dataset
Res_sig_all  <- Res_contrasts %>% 
  subset(padj < 0.1 )  %>%
  droplevels

######select the significant ids
sig_genes_all <- Res_sig_all$gene_id %>% as.character %>%
  unique

### create a tab with the signifcant ones 
Tab_sub_all <- Tab_av[,colnames(Tab_av)]

Tab_sub_all <- match(sig_genes_all,rownames(Tab_sub_all)) %>%
  Tab_sub_all[.,]

#####transform Tab_sub_all in data frame
df <- Tab_sub_all %>% data.frame

######create a column
df$IdRows <- row.names(df)

#######select only the significant genes from L3 and NB
Res_sig  <- Res_contrasts %>% 
  subset(padj < 0.1 )  %>%
  dplyr::filter(.data = .,Contrast %>% grepl(pattern = "_L3_FullSynCom")) %>%
  droplevels

sig_genes <- Res_sig$gene_id %>% as.character %>%
  unique

### Map to entrez id 
Tab_sub <- Tab_av[,colnames(Tab_av) %>% grep(pattern = "^L3")]

Tab_sub <- match(sig_genes,rownames(Tab_sub)) %>%
  Tab_sub[.,]

######rearrange the columns
Tab_L3 <- Tab_sub[,c(5,1,2,3,4)]

######change the columns
colnames(Tab_L3) <- colnames(Tab_L3) %>%
  gsub(pattern = "^L[345]_",replacement = "") %>% 
  factor(levels = c("NB","FullSynCom","L3","L4_up","L5"))

######rename the columns
colnames(Tab_L3)[which(colnames(Tab_L3) == 'FullSynCom')] <- '+SynCom'
colnames(Tab_L3)[which(colnames(Tab_L3) == 'L3')] <- '+SynCom-L3'
colnames(Tab_L3)[which(colnames(Tab_L3) == 'L4_up')] <- '+SynCom-L4'
colnames(Tab_L3)[which(colnames(Tab_L3) == 'L5')] <- '+SynCom-L5'

######Heatmap
res_heatmap <- chibi.heatmap_l3(Tab = Tab_L3,
                             dist_method_rows = "euclidean",
                             hclust_method_rows = "ward.D",
                             range_fill_heatmap = c(-1,1),
                             k_rows = 11,panel_spacing = 0.1,
                             k_cols = 1,size_axis_text_row = 6,axis_ticks_row = F,
                             mtheme = theme(plot.title = element_text(hjust = 0.5, size = 20),
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
                                            legend.key.width = unit(1, 'cm')))

a <- res_heatmap$heatmap
a

######select the clusters from heatmap
df_clust_rows <- res_heatmap$df_clust_rows

######select the dataset
df <- res_heatmap$melted

######merge with df
df_final <- merge(df_clust_rows, df, by='IdRows') 

######save the Tab_sub_all
#write.table(df_final, file='../cleandata/TableS3.csv', sep=',',
            #row.names = F, col.names = T)

######select the L3_NB and L3_FullSynCom
Tab_L3 <- Tab_sub[,c(5,1)]

######change the column name
colnames(Tab_L3) <- c('NB', '+SynCom')

#####plot the heatmap
res_heatmap_l3 <- chibi.heatmap_l3(Tab = Tab_L3,
                             dist_method_rows = "euclidean",
                             hclust_method_rows = "ward.D",
                             range_fill_heatmap = c(-2,2),
                             k_rows = 11,panel_spacing = 0.1,
                             k_cols = 1,size_axis_text_row = 6,
                             axis_ticks_row = F,
                             legend_position = "none",
                             mtheme = theme(plot.title = element_text(hjust = 0.5, 
                                                                      size = 20),
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
                                            legend.key.width = unit(1, 'cm')))

b <- res_heatmap_l3$heatmap
b

######save as pdf
oh.save.pdf(p = b,
            outname = "fig4-heatmap_L3_RNAseq_ohchibi_with_CRs.pdf",
            outdir = "../figures/",width = 12,height = 20)

#######select only the significant genes from L4
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

c <- res_heatmap_l4$heatmap
c

######select the clusters from heatmap
df_clust_rows_l4 <- res_heatmap_l4$df_clust_rows

######merge with df
df_final_l4 <- merge(df_clust_rows_l4, df, by='IdRows') 

######save the Tab_sub_all
write.table(df_final_l4, file='../cleandata/TableS3_L4.csv', sep=',',
            row.names = F, col.names = T)


######select the L3_NB and L3_FullSynCom
Tab_L4 <- Tab_sub_l4[,c(5,1)]

######change column name
colnames(Tab_L4) <- c('NB', '+SynCom')

#####plot the heatmap
res_heatmap_l4 <- chibi.heatmap_l4(Tab = Tab_L4,
                                dist_method_rows = "euclidean",
                                hclust_method_rows = "ward.D",
                                range_fill_heatmap = c(-2,2),
                                k_rows = 8,panel_spacing = 0.2,
                                k_cols = 2,size_axis_text_row = 6,
                                axis_ticks_row = T,
                                legend_position = "bottom",
                                mtheme = theme(plot.title = element_text(hjust = 0.5, size = 20),
                                               axis.title.y =element_blank(),
                                               axis.ticks.y=element_blank(),
                                               axis.ticks.x = element_blank(),
                                               axis.text.y = element_blank(),
                                               axis.text.x = element_text(size=20, 
                                                                          angle=360,
                                                                          hjust = 0.5),
                                               panel.background = element_blank(),
                                               panel.spacing = unit (0.5, 'lines'),
                                               strip.background = element_blank(),
                                               strip.text = element_text(size = 30),
                                               strip.text.y = element_text(size=15),
                                               legend.text = element_text(size=10),
                                               legend.key.width = unit(1, 'cm')))

d <- res_heatmap_l4$heatmap
d

#######save the heatmap
oh.save.pdf(p = d,
            outname = "fig4-heatmap_L4_RNAseq_ohchibi_CRs.pdf",
            outdir = "../figures/",width = 12,height = 20)

#######select only the significant genes from L3 and NB
Res_sig_l5  <- Res_contrasts %>% 
  subset(padj < 0.1 )  %>%
  dplyr::filter(.data = .,Contrast %>% grepl(pattern = "_L5_FullSynCom")) %>%
  droplevels

sig_genes_l5 <- Res_sig_l5$gene_id %>% as.character %>%
  unique

### Map to entrez id 
Tab_sub_l5 <- Tab_av[,colnames(Tab_av) %>% grep(pattern = "^L5")]

Tab_sub_l5 <- match(sig_genes_l5,rownames(Tab_sub_l5)) %>%
  Tab_sub_l5[.,]

######select the L3_NB and L3_FullSynCom
Tab_l5 <- Tab_sub_l5[,c(5,1)]

######change the colnames
colnames(Tab_l5) <- c('NB', '+SynCom')

######test the optimal number of cluster
library(NbClust)
NbClust(data = Tab_l5, diss = NULL, distance = "euclidean",
        min.nc = 2, max.nc = 11, method = 'ward.D')

#####plot the heatmap
res_heatmap_l5 <- chibi.heatmap_l5(Tab = Tab_l5,
                                dist_method_rows = "euclidean",
                                hclust_method_rows = "ward.D",
                                range_fill_heatmap = c(-2,2),
                                k_rows = 4,panel_spacing = 0.2,
                                k_cols = 2,size_axis_text_row = 6,
                                axis_ticks_row = T,
                                legend_position = "none",
                                tree_scale_factor_rows = c(0.009,0.009),
                                egg_widths = c(0.1,1),
                                mtheme = theme(plot.title = element_text(hjust = 0.5, size = 20),
                                               axis.title.y =element_blank(),
                                               axis.ticks.y=element_blank(),
                                               axis.ticks.x = element_blank(),
                                               axis.text.y = element_blank(),
                                               axis.text.x = element_text(size=20, 
                                                                          angle=360,
                                                                          hjust = 0.5),
                                               panel.background = element_blank(),
                                               panel.spacing = unit (0.5, 'lines'),
                                               strip.background = element_blank(),
                                               strip.text = element_text(size = 30),
                                               strip.text.y = element_text(size=15),
                                               legend.text = element_text(size=10),
                                               legend.key.width = unit(1, 'cm')))

f <- res_heatmap_l5$heatmap
f


######select the clusters from heatmap
df_clust_rows_l5 <- res_heatmap_l5$df_clust_rows

######merge with df
df_final_l5 <- merge(df_clust_rows_l5, df, by='IdRows') 

######save the Tab_sub_all
write.table(df_final_l5, file='../cleandata/TableS3_L5.csv', sep=',',
            row.names = F, col.names = T)

#######save the heatmap
oh.save.pdf(p = f,
            outname = "fig4-heatmap_l5_RNAseq_ohchibi_CRs.pdf",
            outdir = "../figures/",width = 12,height = 20)

######composition
composition <- grid.arrange(b,d,f, nrow = 1)

#save
g <- arrangeGrob(b, d,f, nrow = 1) #generates g
ggsave(file="../figures/fig4-heatmap_RNAseq_all_leaves_nb_syncom_ohchibi_CRs.pdf", g)

pdf("../figures/fig4-heatmap_RNAseq_all_leaves_nb_syncom_ohchibi.pdf", 
    width = 12, height = 20) # Open a new pdf file
grid.arrange(b, d,f, nrow = 1) # Write the grid.arrange in the file
dev.off()

#######save the heatmap
oh.save.pdf(p = composition,
            outname = "fig4-heatmap_RNAseq_all_leaves_nb_syncom_ohchibi_CRs.pdf",
            outdir = "../figures/",width = 12,height = 20)
