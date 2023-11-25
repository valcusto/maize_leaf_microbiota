######load packages
library(ohchibi)
library(dplyr)
library(ggh4x)
library(emmeans)
library(ggtree)
library(paletteer)
library(scales)
library(egg)

######set seed
set.seed(130816)

######set work directory
setwd("C:/Users/Castrillo Lab/OneDrive - The University of Nottingham/Publication 2019/Paper preparation/manuscripts/scripts")

######helper functions
source('0-Clean_up_plots.R')
source('0_7-chibi.heatmap.R')
source('0_11-theme_ohchibi_2.R')
######paleta leaf
paleta_leaf <- c('#7b3294','#c2a5cf','#a6dba0')
names(paleta_leaf) <- c('L3', 'L4', 'L5')

### Read the contrasts and dataset
Res_contrasts <- readRDS(file = "../cleandata/res_rnaseq_contrasts_syncom_do.RDS")
Dat_av  <- readRDS(file = "../cleandata/res_tab_av_rnaseq.RDS")
Tab_av <- Dat_av$Tab_av

#######select only the significant genes from L3 and NB
Res_sig  <- Res_contrasts %>% 
  subset(padj < 0.1 )  %>%
  dplyr::filter(.data = .,Contrast %>% grepl(pattern = "_L3_FullSynCom")) %>%
  droplevels

######defense genes
df_gene <- read.csv('../cleandata/defense_genes_in_CR468_CR8C.csv')

######select only defense genes
df_defense <- df_gene %>%
  subset(category == 'defense') %>% droplevels

sig_genes <- df_defense$IdRows %>% as.character %>% unique

######select the defense genes
Tab_av <- match(sig_genes,rownames(Tab_av)) %>%
  Tab_av[.,]

res_heatmap <- chibi.heatmap(Tab = Tab_av,
                             dist_method_rows = "euclidean",
                             hclust_method_rows = "ward.D",
                             range_fill_heatmap = c(-2,2),
                             k_rows = 11,panel_spacing = 0.2,
                             k_cols = 5,size_axis_text_row = 6,axis_ticks_row = T)

a <- res_heatmap$heatmap
a

######select the L3_NB and L3_FullSynCom
Tab_all <- Tab_av[,c(1,5,6,10,11,15)]

distfun <- function(x,method) vegan::vegdist(x = x, method = "euclidean")
dist_mun <- distfun(Tab_all)

#######Cluster and order the family patterns
mclust <- hclust(d = dist_mun,method = "ward.D")
order <- mclust$order %>% mclust$labels[.]

###### View the dendogram
mclust %>%plot

#####plot the heatmap
res_heatmap_all <- chibi.heatmap(Tab = Tab_all,
                             dist_method_rows = "euclidean",
                             hclust_method_rows = "ward.D",
                             range_fill_heatmap = c(-2,2),
                             k_rows = 3,panel_spacing = 0.2,
                             k_cols = 6,size_axis_text_row = 6,axis_ticks_row = T)

b <- res_heatmap_all$heatmap
b

######clustrows
df_clustrows <- res_heatmap_all$df_clust_rows

######melt the table
melted <- Tab_all %>% melt

######name the columns
colnames(melted) <- c('IdRows','group','zscore')

######separate the NumLeaf and Treatment
melted$Treatment <- melted$group %>%
  gsub(pattern = '^L[345]_', replacement = '') %>% as.character
melted$NumLeaf <- 'L3'
melted$NumLeaf[melted$group %>% grep(pattern = "^L4")] <- "L4"
melted$NumLeaf[melted$group %>% grep(pattern = "^L5")] <- "L5"

######change the name of fullsyncom
melted$Treatment[which(melted$Treatment == 'FullSynCom')] <- '+SynCom'

######order the Treatment columns
melted$Treatment <- melted$Treatment %>% 
  factor(levels=c('NB', '+SynCom'))

melted$NumLeaf <- melted$NumLeaf %>% 
  factor(levels=c('L3', 'L4', 'L5'))

#####merge with the clusters
df_new_final <- merge(df_clustrows, melted, by='IdRows')

######plot
######plot the heatmap
p_heatmap <- ggplot(data =df_new_final,aes(x = Treatment,y = IdRows)) +
  geom_raster(aes(fill = zscore)) +
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-1,1),oob = squish)+
  facet_nested(ClusterRows~NumLeaf+Treatment,space = "free",scales = "free")+
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
        axis.text.y = element_text(size = 8),
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
            outname = "fig4C-heatmap_defense_genes_Jia_clusters_final.pdf",
            outdir = "../figures/",width = 10,height = 12)


######calculate the mean 
df_final <- aggregate(melted,zscore ~ group, mean)

######order the factor
df_final$group <- df_final$group %>%
  factor(levels = c('L3_NB', 'L3_FullSynCom',
                    'L4_NB', 'L4_FullSynCom',
                    'L5_NB', 'L5_FullSynCom'))

######separate the number leaf
df_final$Treatment <- df_final$group %>%
  gsub(pattern = '^L[345]_', replacement = '') %>% as.character
df_final$NumLeaf[df_final$group %>% grep(pattern = "^L3")] <- "L3"
df_final$NumLeaf[df_final$group %>% grep(pattern = "^L4")] <- "L4"
df_final$NumLeaf[df_final$group %>% grep(pattern = "^L5")] <- "L5"

######order the factor
df_final$Treatment <- df_final$Treatment %>%
  factor(levels = c('NB', 'FullSynCom'))

###create the unique names for the for loop
mleaf <- melted$NumLeaf %>% unique

###Empty dataset
Res_em <- NULL
Res_pval <- NULL

###perform the for loop for ttest
for (i in mleaf) {
  df_leaf <- melted %>% subset(NumLeaf == i) %>% droplevels
    #linear model
   m1 <- lm(formula = zscore ~ Treatment,data = df_leaf)
    ######
  m1_em <- emmeans(m1,pairwise ~ Treatment,
                    adjust = "none")
  df_em <- m1_em$emmeans %>% as.data.frame
  df_pval <- m1_em$contrasts %>% as.data.frame
  df_pval$NumLeaf<- i
  df_em$NumLeaf  <- i
    #Extract pvalue
  Res_em <- rbind(Res_em,df_em)
  Res_pval <- rbind(Res_pval,df_pval)
}

#####add significance on the column
pthres = 0.05
Res_pval$Significance <- "NS"
Res_pval$Significance[which(Res_pval$p.value < pthres)] <- 'p<0.05'

######select the column 7 to 10
pval <- Res_pval[,c(7:8)]

#####merge Res_em and pval
Res <- merge(Res_em, pval, by='NumLeaf')

######factor
Res$Treatment <- Res$Treatment %>% as.character
Res$Treatment[which(Res$Treatment == 'FullSynCom')] <- '+SynCom'

#####order the Treatment factor
Res$Treatment <- Res$Treatment %>%
  factor(levels = c('NB', '+SynCom'))

######Step10:Plot the informations
p1 <- ggplot(Res, aes(x=Treatment, y=emmean)) + 
  geom_line(aes(group = NumLeaf,color = Significance),size = 1) +
  geom_point(shape = 21,aes(fill = NumLeaf),size = 10) +
  theme_ohchibi(size_panel_border = 2) +
  scale_fill_manual(name = 'NumLeaf',values = paleta_leaf) +
  scale_color_manual(values = c("#D9D9D9","black")) +
  ylim(-0.4,0.4)+
  ggtitle('Defence genes')+
  xlab(label = NULL) + ylab(label = "Standardized expression") +
  theme_ohchibi( size_axis_title.x = 22,
                 size_axis_title.y = 22,
                 legend_proportion_size = 2,
                 size_title_text = 30,
                 size_legend_text = 20,
                 size_panel_border = 1.5,
                 size_lines_panel = 0,
                 size_axis_text.y = 22,
                 size_axis_text.x = 22) +
  clean +
  theme(legend.position = "right",
    panel.background = element_rect(fill = "white"),
    strip.text.x = element_text(family = "Helvetica",face = "bold",size = 20),
    plot.title = element_text(size=30, hjust = 0.5, vjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=30),
    axis.text.x = element_text(size=20),
    axis.text.y = element_text(size=20)) 

p1

#######save the heatmap
oh.save.pdf(p = p1,
            outname = "fig4D-expression_defence_genes_from_RNAseq_clusters_Jia.pdf",
            outdir = "../figures/",width = 10,height = 12)
