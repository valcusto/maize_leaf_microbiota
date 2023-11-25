#######load packages
library(ohchibi)
library(ggtree)
library(dplyr)
library(multcomp)
library(emmeans)
library(egg)
library(tidyr)

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

######select the significant expression
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

######heatmap
res_heatmap <- chibi.heatmap(Tab = Tab_sub,
                             dist_method_rows = "euclidean",hclust_method_rows = "ward.D",
                             range_fill_heatmap = c(-2,2),
                             k_rows = 11,panel_spacing = 0.2,
                             k_cols = 5,size_axis_text_row = 6,axis_ticks_row = T)

a <- res_heatmap$heatmap
a


#oh.save.pdf(p = a,outname = "heatmap_l3.pdf",outdir = "./figures/",width = 12,height = 48)
df_clust_rows <- res_heatmap$df_clust_rows

sig_genes <- df_clust_rows$IdRows

######select only the CR4, CR6 and CR8
mlist_query <- list(
  CR8 = df_clust_rows %>% subset(ClusterRows == "CR8") %$% IdRows %>% as.character,
  CR6 = df_clust_rows %>% subset(ClusterRows == "CR6") %$% IdRows %>% as.character,
  CR4 = df_clust_rows %>% subset(ClusterRows == "CR4") %$% IdRows %>% as.character
)

######select the genes significant genes in L3
toselect <- mlist_query %>% unlist %>% as.character %>% unique

######select the cluster 4,6,8 for the different leaves
melted_all <- Tab_av %>% melt %>%
  dplyr::rename(.data = .,gene_id = Var1,group = Var2)

######create the list with the gene id for the clusters
#toselect_all <- mlist_query %>% unlist %>% as.character %>% unique

######select the genes
melted_sub_all <- melted_all %>%
  dplyr::filter(.data = .,gene_id %in% toselect) %>% droplevels

######add the cluster information
melted_sub_all$ClusterRows <- match(melted_sub_all$gene_id,df_clust_rows$IdRows) %>%
  df_clust_rows$ClusterRows[.]

######create the group
melted_sub_all$group <- melted_sub_all$group %>% 
  factor(levels = c("L3_NB","L3_FullSynCom","L3_L3","L3_L4_up","L3_L5",
                    "L4_NB","L4_FullSynCom","L4_L3","L4_L4_up","L4_L5",
                    "L5_NB","L5_FullSynCom","L5_L3","L5_L4_up","L5_L5"))

melted_sub_all$NumLeaf <- "L3"
melted_sub_all$NumLeaf[melted_sub_all$group %>% grep(pattern = "^L4")] <- "L4"
melted_sub_all$NumLeaf[melted_sub_all$group %>% grep(pattern = "^L5")] <- "L5"
melted_sub_all$NumLeaf <- melted_sub_all$NumLeaf %>% factor

melted_sub_all$Treatment <- melted_sub_all$group %>%
  gsub(pattern = "^L[345]_",replacement = "") %>% 
  factor(levels = c("NB","FullSynCom","L3","L4_up","L5"))

######select only NB and FullSynCom
melted_full <- melted_sub_all %>%
  subset((Treatment == 'FullSynCom')) %>% droplevels
melted_nb <- melted_sub_all %>%
  subset(Treatment == 'NB') %>% droplevels
#melted_nb <- melted_nb %>%
  #subset(NumLeaf == 'L3') %>% droplevels

######combine the two datasets
df_final <- rbind(melted_full, melted_nb)

######rename the information on the tReatment column
df_final$Treatment <- df_final$Treatment %>% as.character
df_final$Treatment[which(df_final$Treatment == 'FullSynCom')] <- '+SynCom'

######order the Treatment columns
df_final$Treatment <- df_final$Treatment %>% 
  factor(levels=c('NB', '+SynCom'))

#######create boxplot
ggplot(data=df_final, mapping = aes(x=group,y=value))+
  geom_boxplot()+
  geom_point()+
  facet_grid(.~ClusterRows)

df_final$group <- df_final$group %>% factor(levels = c('L3_NB', 'L3_FullSynCom',
                                                       'L4_NB', 'L4_FullSynCom',
                                                       'L5_NB', 'L5_FullSynCom'))

######plot the common cluster
ggplot(data=df_final, aes(x=group, y=value, color=group))+
  geom_boxplot()+
  geom_point(shape=21, size=1.5, alpha=1)+
  facet_grid(.~ClusterRows, scales = 'free')

###create the unique names for the for loop
mleaf <- df_final$NumLeaf %>% unique
mcluster <- df_final$ClusterRows %>% unique

###Empty dataset
Res <- NULL

###perform the for loop
for (cluster in mcluster) {
  df_clust <- df_final %>% subset(ClusterRows == cluster) %>% droplevels
  for (leaf in mleaf){
    df_leaf <- df_clust %>% subset(NumLeaf == leaf) %>% droplevels
    #linear model
    m1 <- lm(formula = value ~ Treatment,data = df_leaf)
    ######
    df_em <- emmeans(object = m1,specs = "Treatment")
    cld <- multcomp::cld(df_em, Letters=c("abcdefghi")) %>% data.frame
    cld$ClusterRows <- cluster
    cld$NumLeaf <- leaf
    #combine the results
    Res <- rbind(Res, cld)
  }
}

######chanhe the .group name
colnames(Res)[7] <- 'letters'

######remove space from letters
Res$letters <- Res$letters %>%
  gsub(pattern = ' ', replacement = '')

######create color for the dataset
paleta_syncom <- c('#756bb1','#54278f',
                   '#bcbddc','#9e9ac8',
                   '#99d8c9', '#2ca25f')
names(paleta_syncom) <- c('L3\nNB','L3\n+SynCom',
                          'L4\nNB','L4\n+SynCom',
                          'L5\nNB','L5\n+SynCom')

######transform the group in character
df_final$group <- df_final$group %>% as.character

######correct the group name
df_final$group[which(df_final$group == 'L3_FullSynCom')] <- 'L3\n+SynCom'
df_final$group[which(df_final$group == 'L3_NB')] <- 'L3\nNB'
df_final$group[which(df_final$group == 'L4_FullSynCom')] <- 'L4\n+SynCom'
df_final$group[which(df_final$group == 'L4_NB')] <- 'L4\nNB'
df_final$group[which(df_final$group == 'L5_FullSynCom')] <- 'L5\n+SynCom'
df_final$group[which(df_final$group == 'L5_NB')] <- 'L5\nNB'

######create a group in the Res dataset
Res$group <- paste(Res$NumLeaf, Res$Treatment, sep='_')

######correct the group names
Res$group[which(Res$group == 'L3_+SynCom')] <- 'L3\n+SynCom'
Res$group[which(Res$group == 'L3_NB')] <- 'L3\nNB'
Res$group[which(Res$group == 'L4_+SynCom')] <- 'L4\n+SynCom'
Res$group[which(Res$group == 'L4_NB')] <- 'L4\nNB'
Res$group[which(Res$group == 'L5_+SynCom')] <- 'L5\n+SynCom'
Res$group[which(Res$group == 'L5_NB')] <- 'L5\nNB'

######organize the group
df_final$group <- df_final$group %>%
  factor(levels = c('L3\nNB', 'L3\n+SynCom',
                    'L4\nNB','L4\n+SynCom', 
                    'L5\nNB','L5\n+SynCom'))

######plot the result
p_all <- ggplot(data=df_final, aes(x=group, y=value, color=group))+
  geom_point(shape=21, size=1.5, alpha=1)+
  geom_hline(yintercept = 0)+
  facet_grid(.~ClusterRows, scales = 'free')+
  geom_pointrange(data=Res, aes(y=emmean, ymin=lower.CL, ymax=upper.CL,
                                color=group), 
                  size=1)+
  #geom_text(data = Res, aes(x = group,y = 2,
                            #label = letters),
            #inherit.aes = F,size = 6,family ="Arial",color = "black") +
  scale_color_manual(values = paleta_syncom)+
  ylim(-2,2)+
  xlab(NULL)+ 
  ylab('Standardized expression')+
  clean +
  theme(axis.text.x = element_text(size=12),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.position = 'none')

p_all

######add the significance on the plot
p <- p_all +
  geom_line(data = tibble(x=c(1,2), y=c(1.9, 1.9)), aes(x=x, y=y),
            inherit.aes = FALSE) +
  geom_line(data = tibble(x=c(3,4), y=c(1.9, 1.9)), aes(x=x, y=y),
            inherit.aes = FALSE) +
  geom_line(data = tibble(x=c(5,6), y=c(1.9, 1.9)), aes(x=x, y=y),
            inherit.aes = FALSE) +
  geom_text(data = tibble(x=1.5, y=1.95), aes(x=x, y=y, label = '*'),
            size=6,
            inherit.aes = FALSE) +
  geom_text(data = tibble(x=3.5, y=2), aes(x=x, y=y, label = 'NS'),
            inherit.aes = FALSE) +
  geom_text(data = tibble(x=5.5, y=2), aes(x=x, y=y, label = 'NS'),
            inherit.aes = FALSE)

p

######save figures
oh.save.pdf(p = p,outname = "figS4L_gene_expression_per_cluster_all_leaves_ttest.pdf",
            outdir = "../figures/",width = 20,height = 10)
