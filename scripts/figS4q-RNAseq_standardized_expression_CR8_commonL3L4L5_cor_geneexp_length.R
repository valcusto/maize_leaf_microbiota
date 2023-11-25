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

### Read the contrasts and dataset
Res_contrasts <- readRDS(file = "../cleandata/res_rnaseq_contrasts_syncom_do.RDS")
Dat_av  <- readRDS(file = "../cleandata/res_tab_av_rnaseq.RDS")
Tab_av <- Dat_av$Tab_av

#####################L3
### Correlation 
df_cor <- read.table(file = "../rawdata/res_cor_geneexp_length.tsv")
df_l3 <- df_cor %>% subset(V4 < 0.01) %>% subset(V2 == "L3")

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

######select only the CR8C
mlist_query <- list(
  CR8 = df_clust_rows_l3 %>% subset(ClusterRows == "CR8C") %$% IdRows %>% as.character
)

######define the toselect
toselect <- mlist_query %>% unlist %>% as.character %>% unique

######melt the dataset
melted_all <- Tab_av %>% melt %>%
  dplyr::rename(.data = .,gene_id = Var1,group = Var2)

######select genes belong to the cluster 8C
melted_sub_all <- melted_all %>%
  dplyr::filter(.data = .,gene_id %in% toselect) %>% droplevels

######add the cluster information
melted_sub_all$ClusterRows <- match(melted_sub_all$gene_id,
                                    df_clust_rows_l3$IdRows) %>%
  df_clust_rows_l3$ClusterRows[.]

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
  geom_hline(yintercept = 0,linetype="dashed", color = "black", size = 0.5)+
  facet_grid(.~ClusterRows, scales = 'free')+
  geom_pointrange(data=Res, aes(y=emmean, ymin=lower.CL, ymax=upper.CL,
                                color=group), 
                  size=1)+
  scale_color_manual(values = paleta_syncom)+
  ylim(-1,1)+
  xlab(NULL)+ 
  ylab('Standardized expression')+
  clean +
  theme(axis.text.x = element_text(size=15),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.position = 'none')

p_all

######add the significance on the plot
p <- p_all +
  geom_line(data = tibble(x=c(1,2), y=c(0.9, 0.9)), aes(x=x, y=y),
            inherit.aes = FALSE) +
  geom_line(data = tibble(x=c(3,4), y=c(0.9, 0.9)), aes(x=x, y=y),
            inherit.aes = FALSE) +
  geom_line(data = tibble(x=c(5,6), y=c(0.9, 0.9)), aes(x=x, y=y),
            inherit.aes = FALSE) +
  geom_text(data = tibble(x=1.5, y=0.95), aes(x=x, y=y, label = '*'),
            size=6,
            inherit.aes = FALSE) +
  geom_text(data = tibble(x=3.5, y=1), aes(x=x, y=y, label = 'NS'),
            inherit.aes = FALSE) +
  geom_text(data = tibble(x=5.5, y=1), aes(x=x, y=y, label = 'NS'),
            inherit.aes = FALSE)

p

######save figures
oh.save.pdf(p = p,outname = "figS4q_gene_expression_cor_geneexp_length_all_leaves_ttest.pdf",
            outdir = "../figures/",width = 6,height = 8)
