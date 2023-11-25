#######load packages
library(ohchibi)
library(paletteer)
library(scales)
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

######remove the gene Zm00001eb032600 from the dataset as it is an outlier
df_final <- df_final %>%
  subset(gene_id != 'Zm00001eb032600') %>% droplevels

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

######perform the t test
cr4_l3 <- df_final %>%
  subset(ClusterRows == 'CR4' & NumLeaf == 'L3') %>% droplevels
cr4_l4 <- df_final %>%
  subset(ClusterRows == 'CR4' & NumLeaf == 'L4') %>% droplevels
cr4_l5 <- df_final %>%
  subset(ClusterRows == 'CR4' & NumLeaf == 'L5') %>% droplevels

cr6_l3 <- df_final %>%
  subset(ClusterRows == 'CR6' & NumLeaf == 'L3') %>% droplevels
cr6_l4 <- df_final %>%
  subset(ClusterRows == 'CR6' & NumLeaf == 'L4') %>% droplevels
cr6_l5 <- df_final %>%
  subset(ClusterRows == 'CR6' & NumLeaf == 'L5') %>% droplevels

cr8_l3 <- df_final %>%
  subset(ClusterRows == 'CR8' & NumLeaf == 'L3') %>% droplevels
cr8_l4 <- df_final %>%
  subset(ClusterRows == 'CR8' & NumLeaf == 'L4') %>% droplevels
cr8_l5 <- df_final %>%
  subset(ClusterRows == 'CR8' & NumLeaf == 'L5') %>% droplevels

######remove the column group
cr4_l3 <- cr4_l3[,-2]
cr4_l4 <- cr4_l4[,-2]
cr4_l5 <- cr4_l5[,-2]

cr6_l3 <- cr6_l3[,-2]
cr6_l4 <- cr6_l4[,-2]
cr6_l5 <- cr6_l5[,-2]

cr8_l3 <- cr8_l3[,-2]
cr8_l4 <- cr8_l4[,-2]
cr8_l5 <- cr8_l5[,-2]

######create the columns with the +SynCom and NB
cr4_l3_ttest <- cr4_l3 %>%
  pivot_wider(names_from = Treatment, values_from = value)
cr4_l4_ttest <- cr4_l4 %>%
  pivot_wider(names_from = Treatment, values_from = value)
cr4_l5_ttest <- cr4_l5 %>%
  pivot_wider(names_from = Treatment, values_from = value)

cr6_l3_ttest <- cr6_l3 %>%
  pivot_wider(names_from = Treatment, values_from = value)
cr6_l4_ttest <- cr6_l4 %>%
  pivot_wider(names_from = Treatment, values_from = value)
cr6_l5_ttest <- cr6_l5 %>%
  pivot_wider(names_from = Treatment, values_from = value)

cr8_l3_ttest <- cr8_l3 %>%
  pivot_wider(names_from = Treatment, values_from = value)
cr8_l4_ttest <- cr8_l4 %>%
  pivot_wider(names_from = Treatment, values_from = value)
cr8_l5_ttest <- cr8_l5 %>%
  pivot_wider(names_from = Treatment, values_from = value)

######calculate the t test
Res_cr4_l3_ttest <- t.test(cr4_l3_ttest$NB, cr4_l3_ttest$`+SynCom`, 
                           alternative = "two.sided")
Res_cr4_l4_ttest <- t.test(cr4_l4_ttest$NB, cr4_l4_ttest$`+SynCom`, 
                           alternative = "two.sided")
Res_cr4_l5_ttest <- t.test(cr4_l5_ttest$NB, cr4_l5_ttest$`+SynCom`, 
                           alternative = "two.sided")

Res_cr6_l3_ttest <- t.test(cr6_l3_ttest$NB, cr6_l3_ttest$`+SynCom`, 
                           alternative = "two.sided")
Res_cr6_l4_ttest <- t.test(cr6_l4_ttest$NB, cr6_l4_ttest$`+SynCom`, 
                           alternative = "two.sided")
Res_cr6_l5_ttest <- t.test(cr6_l5_ttest$NB, cr6_l5_ttest$`+SynCom`, 
                           alternative = "two.sided")

Res_cr8_l3_ttest <- t.test(cr8_l3_ttest$NB, cr8_l3_ttest$`+SynCom`, 
                           alternative = "two.sided")
Res_cr8_l4_ttest <- t.test(cr8_l4_ttest$NB, cr8_l4_ttest$`+SynCom`, 
                           alternative = "two.sided")
Res_cr8_l5_ttest <- t.test(cr8_l5_ttest$NB, cr8_l5_ttest$`+SynCom`, 
                           alternative = "two.sided")

######create a dataframe with the results
df_ttest_l3 <- data.frame(ClusterRows = c('CR4', 'CR6', 'CR8'),
                       NumLeaf = c('L3', 'L3', 'L3'),
                       pvalue= c(Res_cr4_l3_ttest$p.value,
                                 Res_cr6_l3_ttest$p.value,
                                 Res_cr8_l3_ttest$p.value))

df_ttest_l4 <- data.frame(ClusterRows = c('CR4', 'CR6', 'CR8'),
                          NumLeaf = c('L4', 'L4', 'L4'),
                          pvalue= c(Res_cr4_l4_ttest$p.value,
                                    Res_cr6_l4_ttest$p.value,
                                    Res_cr8_l4_ttest$p.value))

df_ttest_l5 <- data.frame(ClusterRows = c('CR4', 'CR6', 'CR8'),
                          NumLeaf = c('L5', 'L5', 'L5'),
                          pvalue= c(Res_cr4_l5_ttest$p.value,
                                    Res_cr6_l5_ttest$p.value,
                                    Res_cr8_l5_ttest$p.value))

######bind the datasets
df_ttest <- rbind(df_ttest_l3, df_ttest_l4, df_ttest_l5)

#####add significance column
df_ttest$significance <- 'NS'
df_ttest$significance[which(df_ttest$pvalue < 0.05)] <- 'Significant'

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

######select only CR4
cr4 <- df_final %>%
  subset(ClusterRows == 'CR4') %>% droplevels
cr6 <- df_final %>%
  subset(ClusterRows == 'CR6') %>% droplevels
cr8 <- df_final %>%
  subset(ClusterRows == 'CR8') %>% droplevels

######select CRs for the Res
Res_cr4 <- Res %>%
  subset(ClusterRows == 'CR4') %>% droplevels
Res_cr6 <- Res %>%
  subset(ClusterRows == 'CR6') %>% droplevels
Res_cr8 <- Res %>%
  subset(ClusterRows == 'CR8') %>% droplevels

######plot only CR4
p_cr4 <- ggplot(data=cr4, aes(x=group, y=value, color=group))+
  geom_point(shape=21, size=1.5, alpha=1)+
  geom_hline(yintercept = 0,linetype="dashed", color = "black", size = 0.5)+
  facet_grid(.~ClusterRows, scales = 'free')+
  geom_pointrange(data=Res_cr4, aes(y=emmean, ymin=lower.CL, ymax=upper.CL,
                                color=group), 
                  size=1)+
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

p_cr4

p_cr6 <- ggplot(data=cr6, aes(x=group, y=value, color=group))+
  geom_point(shape=21, size=1.5, alpha=1)+
  geom_hline(yintercept = 0,linetype="dashed", color = "black", size = 0.5)+
  facet_grid(.~ClusterRows, scales = 'free')+
  geom_pointrange(data=Res_cr6, aes(y=emmean, ymin=lower.CL, ymax=upper.CL,
                                    color=group), 
                  size=1)+
  scale_color_manual(values = paleta_syncom)+
  ylim(-2,2)+
  xlab(NULL)+ 
  ylab('Standardized expression')+
  clean +
  theme(axis.text.x = element_text(size=12),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'none')

p_cr6

p_cr8 <- ggplot(data=cr8, aes(x=group, y=value, color=group))+
  geom_point(shape=21, size=1.5, alpha=1)+
  geom_hline(yintercept = 0,linetype="dashed", color = "black", size = 0.5)+
  facet_grid(.~ClusterRows, scales = 'free')+
  geom_pointrange(data=Res_cr8, aes(y=emmean, ymin=lower.CL, ymax=upper.CL,
                                    color=group), 
                  size=1)+
  scale_color_manual(values = paleta_syncom)+
  ylim(-2,2)+
  xlab(NULL)+ 
  ylab('Standardized expression')+
  clean +
  theme(axis.text.x = element_text(size=12),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'none')

p_cr8

######add the significance on the plot
p_cr4_final <- p_cr4 +
  geom_line(data = tibble(x=c(1,2), y=c(1.9, 1.9)), aes(x=x, y=y),
            inherit.aes = FALSE) +
  geom_line(data = tibble(x=c(3,4), y=c(1.9, 1.9)), aes(x=x, y=y),
            inherit.aes = FALSE) +
  geom_line(data = tibble(x=c(5,6), y=c(1.9, 1.9)), aes(x=x, y=y),
            inherit.aes = FALSE) +
  geom_text(data = tibble(x=1.5, y=1.95), aes(x=x, y=y, label = '*'),
            size=6,
            inherit.aes = FALSE) +
  geom_text(data = tibble(x=3.5, y=1.95), aes(x=x, y=y, label = '*'),
            size=6,
            inherit.aes = FALSE) +
  geom_text(data = tibble(x=5.5, y=2), aes(x=x, y=y, label = 'NS'),
            inherit.aes = FALSE)

p_cr4_final

p_cr6_final <- p_cr6 +
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

p_cr6_final

p_cr8_final <- p_cr8 +
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

p_cr8_final

######arrange the panels
composition <- egg::ggarrange(p_cr4_final, p_cr6_final, p_cr8_final,
                             nrow=1)

######save figures
oh.save.pdf(p = composition,outname = "figS4L_gene_expression_per_cluster_all_leaves_ttest_final.pdf",
            outdir = "../figures/",width = 15,height = 10)
