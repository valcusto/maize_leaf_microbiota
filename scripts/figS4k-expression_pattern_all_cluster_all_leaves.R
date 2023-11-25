#######load packages
library(ohchibi)
library(car)
library(FSA)
library(rcompanion)
library(boot)
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
source('0-Clean_up_plots.R')
source('0_6-theme_ohchibi.R')
source('0_7-chibi.heatmap.R')
Bmean <- function(data, indices) {
  d <- data[indices] # allows boot to select sample 
  return(mean(d))
}


### Read the contrasts and dataset
Res_contrasts <- readRDS(file = "../cleandata/res_rnaseq_contrasts_syncom_do.RDS")
Dat_av  <- readRDS(file = "../cleandata/res_tab_av_rnaseq.RDS")
Tab_av <- Dat_av$Tab_av

######select the significant expression
Res_sig  <- Res_contrasts %>% 
  subset(padj < 0.1)  %>%
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

######select the clusters
df_clust_rows <- res_heatmap$df_clust_rows

######transform the Tab_sub_l3 in data frame
df_l3 <- Tab_sub %>% data.frame

######create the gene_id column
df_l3$gene_id <- row.names(df_l3)

######clean the rows
row.names(df_l3) <- NULL

########select the clusters
mlist_query <- list(
  CR1 = df_clust_rows %>% subset(ClusterRows == "CR1") %$% IdRows %>% as.character,
  CR2 = df_clust_rows %>% subset(ClusterRows == "CR2") %$% IdRows %>% as.character,
  CR3 = df_clust_rows %>% subset(ClusterRows == "CR3") %$% IdRows %>% as.character,
  CR4 = df_clust_rows %>% subset(ClusterRows == "CR4") %$% IdRows %>% as.character,
  CR5 = df_clust_rows %>% subset(ClusterRows == "CR5") %$% IdRows %>% as.character,
  CR6 = df_clust_rows %>% subset(ClusterRows == "CR6") %$% IdRows %>% as.character,
  CR7 = df_clust_rows %>% subset(ClusterRows == "CR7") %$% IdRows %>% as.character,
  CR8 = df_clust_rows %>% subset(ClusterRows == "CR8") %$% IdRows %>% as.character,
  CR9 = df_clust_rows %>% subset(ClusterRows == "CR9") %$% IdRows %>% as.character,
  CR10 = df_clust_rows %>% subset(ClusterRows == "CR10") %$% IdRows %>% as.character,
  CR11 = df_clust_rows %>% subset(ClusterRows == "CR11") %$% IdRows %>% as.character
)

######select the genes significant genes in L3
toselect <- mlist_query %>% unlist %>% as.character %>% unique

######select the cluster 4,6,8 for the different leaves
melted_all <- Tab_av %>% melt %>%
  dplyr::rename(.data = .,gene_id = Var1,group = Var2)

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

######rename the information on the tReatment column
melted_sub_all$Treatment <- melted_sub_all$Treatment %>% as.character
melted_sub_all$Treatment[which(melted_sub_all$Treatment == 'FullSynCom')] <- '+SynCom'
melted_sub_all$Treatment[which(melted_sub_all$Treatment == 'L3')] <- '+SynCom-L3'
melted_sub_all$Treatment[which(melted_sub_all$Treatment == 'L4_up')] <- '+SynCom-L4'
melted_sub_all$Treatment[which(melted_sub_all$Treatment == 'L5')] <- '+SynCom-L5'

######order the Treatment columns
melted_sub_all$Treatment <- melted_sub_all$Treatment %>% 
  factor(levels=c('NB', '+SynCom', '+SynCom-L3', '+SynCom-L4', '+SynCom-L5'))

#######create boxplot
ggplot(data=melted_sub_all, mapping = aes(x=group,y=value))+
  geom_boxplot()+
  geom_point()+
  facet_grid(NumLeaf~ClusterRows)

######select only L3
df_final_l3 <- melted_sub_all %>% subset(NumLeaf == 'L3') %>% droplevels

######check the expression
ggplot(data=df_final_l3, mapping = aes(x=group,y=value))+
  geom_boxplot()+
  geom_point()+
  facet_grid(NumLeaf~ClusterRows)

#######create a group with group and ClusterRows
df_final_l3$group1 <- paste(df_final_l3$ClusterRows, df_final_l3$group, sep='_')

######transfor the group as factor
df_final_l3$group1 <- df_final_l3$group1 %>% 
  factor()

######statistical assumptions
#1.interdependence of the variables

#2.variance homogeneity
bartlett.test(value ~ group1, data = df_final_l3)
#p-value<0.05. Reject the null hypothesis. The variance are different

#3.normality
qqnorm(df_final_l3$value)
qqline(df_final_l3$value, col='red')

###create the unique names for the for loop
mcluster <- df_final_l3$ClusterRows %>% unique

###Empty dataset
Res_letter <- NULL

###perform the for loop
for (cluster in mcluster){
  df_clust <- df_final_l3 %>% subset(ClusterRows == cluster) %>% droplevels
  #kruskal wallis
  m1 <- kruskal.test(value ~ group1,data = df_clust)
  ######
  DT = dunnTest(value ~ group1,data = df_clust,method="bh")
  PT = DT$res
  PT
  df_letter <- cldList(P.adj ~ Comparison,
                       data = PT, remove.zero = F,
                       threshold = 0.05)
  df_letter$ClusterRows <- cluster
  Res_letter <- rbind(Res_letter, df_letter)
}

######calculate the mean and the confidence interval for each analysis
###create the unique names for the for loop
mgroup1 <- df_final_l3$group1 %>% unique
###Empty dataset
Res_mean <- NULL
Res_ci <- NULL

###perform the for loop
for (i in mgroup1){
  df_group1 <- df_final_l3 %>% subset(group1 == i) %>% droplevels
  #calculate means
  df_mean <- mean(df_group1$value) %>% as.data.frame
  ###calculate ci
  # bootstrapping with 1000 replications 
  results <- boot(data=df_group1$value, statistic=Bmean, R=1000)
  # get 95% confidence interval 
  results_boot <- boot.ci(results, type=c("norm", "basic", "perc", "bca"))
  ###create a data frame for ci
  df_ci <- results_boot$bca %>% as.data.frame
  #add the leaf name on data frame
  df_mean$group1 <- i
  df_ci$group1 <- i
  #combine the results
  Res_mean <- rbind(Res_mean, df_mean)
  Res_ci <- rbind(Res_ci, df_ci)
}

######create the dataset formean
colnames(Res_mean)[1]<-'mean'
###select column from 4to6
Res_ci <- Res_ci[,4:6]
###rename the first column
colnames(Res_ci)[1]<-'lower'
colnames(Res_ci)[2]<-'upper'

######merge the two datasets
res <- merge(Res_mean, Res_ci, by='group1')

######remove res between group names
res$group1 <- res$group1 %>% gsub(pattern = ' ', replacement = '')

######split group1
res$group1 <- res$group1 %>% as.character
res$ClusterRows <- lapply(strsplit(res$group1, '_'), `[`,1) %>% as.character
res$NumLeaf <- lapply(strsplit(res$group1, '_'), `[`, 2) %>% as.character
res$Treatment <- lapply(strsplit(res$group1, '_'), `[`, 3) %>% as.character

######rename the information on the tReatment column
res$Treatment[which(res$Treatment == 'FullSynCom')] <- '+SynCom'
res$Treatment[which(res$Treatment == 'L3')] <- '+SynCom-L3'
res$Treatment[which(res$Treatment == 'L4')] <- '+SynCom-L4'
res$Treatment[which(res$Treatment == 'L5')] <- '+SynCom-L5'

######change the first column name 
colnames(Res_letter)[1] <- 'group1'

######merge the dataset
res <- merge(res, Res_letter, by= c('group1', 'ClusterRows'))

######organize the CRs
df_final_l3$ClusterRows <- df_final_l3$ClusterRows %>%
  factor(levels = c('CR1', 'CR2', 'CR3', 'CR4', 'CR5',
                    'CR6', 'CR7', 'CR8', 'CR9', 'CR10', 'CR11'))

res$ClusterRows <- res$ClusterRows %>%
  factor(levels = c('CR1', 'CR2', 'CR3', 'CR4', 'CR5',
                    'CR6', 'CR7', 'CR8', 'CR9', 'CR10', 'CR11'))

######define the color for the dropout
paleta_syncom <- c('#d0d1e6','#a6bddb','#67a9cf',
                   '#3690c0','#02818a','#016450')
names(paleta_syncom) <- c('NB','+SynCom','+SynCom-L1&L2',
                          '+SynCom-L3','+SynCom-L4','+SynCom-L5')

######plot the result
p_l3 <- ggplot(data=df_final_l3, aes(x=Treatment, y=value, color=Treatment))+
  geom_point(shape=21, size=1.5, alpha=1)+
  facet_grid(.~ClusterRows, scales = 'free')+
  geom_pointrange(data=res, aes(y=mean, ymin=lower, ymax=upper,
                                color=Treatment), 
                  size=1)+
  geom_hline(yintercept=0, linetype="dashed", color = "black", size = 0.5)+
  geom_text(data = res, aes(x = Treatment,y = 2.5,
                            label = Letter),
            inherit.aes = F,size = 6,family ="Arial",color = "black") +
  #scale_color_manual(values = paleta_syncom)+
  xlab(NULL)+ 
  ylab('Standardized expression')+
  scale_color_manual(values = paleta_syncom)+
  clean +
  theme(axis.text.x = element_text(size=12, angle=45, hjust=1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.position = 'none')

p_l3

######save the figure
oh.save.pdf(p = p_l3,
            outname = "figS4B_gene_expression_L3_all_clusters.pdf",
            outdir = "../figures/",width = 30,height = 10)


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

######heatmap
res_heatmap_l4 <- chibi.heatmap(Tab = Tab_sub_l4,
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

b <- res_heatmap_l4$heatmap
b

######select the clusters
df_clust_rows_l4 <- res_heatmap_l4$df_clust_rows

######transform the Tab_sub_l3 in data frame
df_l4 <- Tab_sub_l4 %>% data.frame

######create the gene_id column
df_l4$gene_id <- row.names(df_l4)

######clean the rows
row.names(df_l4) <- NULL

########select the clusters
mlist_query_l4 <- list(
  CR1 = df_clust_rows_l4 %>% subset(ClusterRows == "CR1") %$% IdRows %>% as.character,
  CR2 = df_clust_rows_l4 %>% subset(ClusterRows == "CR2") %$% IdRows %>% as.character,
  CR3 = df_clust_rows_l4 %>% subset(ClusterRows == "CR3") %$% IdRows %>% as.character,
  CR4 = df_clust_rows_l4 %>% subset(ClusterRows == "CR4") %$% IdRows %>% as.character,
  CR5 = df_clust_rows_l4 %>% subset(ClusterRows == "CR5") %$% IdRows %>% as.character,
  CR6 = df_clust_rows_l4 %>% subset(ClusterRows == "CR6") %$% IdRows %>% as.character,
  CR7 = df_clust_rows_l4 %>% subset(ClusterRows == "CR7") %$% IdRows %>% as.character,
  CR8 = df_clust_rows_l4 %>% subset(ClusterRows == "CR8") %$% IdRows %>% as.character
)

######select the genes significant genes in L3
toselect_l4 <- mlist_query_l4 %>% unlist %>% as.character %>% unique

######select the cluster 4,6,8 for the different leaves
melted_all_l4 <- Tab_av %>% melt %>%
  dplyr::rename(.data = .,gene_id = Var1,group = Var2)

######select the genes
melted_sub_all_l4 <- melted_all_l4 %>%
  dplyr::filter(.data = .,gene_id %in% toselect_l4) %>% droplevels

######add the cluster information
melted_sub_all_l4$ClusterRows <- match(melted_sub_all_l4$gene_id,
                                       df_clust_rows_l4$IdRows) %>%
  df_clust_rows_l4$ClusterRows[.]

######create the group
melted_sub_all_l4$group <- melted_sub_all_l4$group %>% 
  factor(levels = c("L3_NB","L3_FullSynCom","L3_L3","L3_L4_up","L3_L5",
                    "L4_NB","L4_FullSynCom","L4_L3","L4_L4_up","L4_L5",
                    "L5_NB","L5_FullSynCom","L5_L3","L5_L4_up","L5_L5"))

melted_sub_all_l4$NumLeaf <- "L3"
melted_sub_all_l4$NumLeaf[melted_sub_all_l4$group %>% grep(pattern = "^L4")] <- "L4"
melted_sub_all_l4$NumLeaf[melted_sub_all_l4$group %>% grep(pattern = "^L5")] <- "L5"
melted_sub_all_l4$NumLeaf <- melted_sub_all_l4$NumLeaf %>% factor

melted_sub_all_l4$Treatment <- melted_sub_all_l4$group %>%
  gsub(pattern = "^L[345]_",replacement = "") %>% 
  factor(levels = c("NB","FullSynCom","L3","L4_up","L5"))

######rename the information on the tReatment column
melted_sub_all_l4$Treatment <- melted_sub_all_l4$Treatment %>% as.character
melted_sub_all_l4$Treatment[which(melted_sub_all_l4$Treatment == 'FullSynCom')] <- '+SynCom'
melted_sub_all_l4$Treatment[which(melted_sub_all_l4$Treatment == 'L3')] <- '+SynCom-L3'
melted_sub_all_l4$Treatment[which(melted_sub_all_l4$Treatment == 'L4_up')] <- '+SynCom-L4'
melted_sub_all_l4$Treatment[which(melted_sub_all_l4$Treatment == 'L5')] <- '+SynCom-L5'

######order the Treatment columns
melted_sub_all_l4$Treatment <- melted_sub_all_l4$Treatment %>% 
  factor(levels=c('NB', '+SynCom', '+SynCom-L3', '+SynCom-L4', '+SynCom-L5'))

#######create boxplot
ggplot(data=melted_sub_all_l4, mapping = aes(x=group,y=value))+
  geom_boxplot()+
  geom_point()+
  facet_grid(NumLeaf~ClusterRows)

######select only L3
df_final_l4 <- melted_sub_all_l4 %>% subset(NumLeaf == 'L4') %>% droplevels

######check the expression
ggplot(data=df_final_l4, mapping = aes(x=group,y=value))+
  geom_boxplot()+
  geom_point()+
  facet_grid(NumLeaf~ClusterRows)

#######create a group with group and ClusterRows
df_final_l4$group1 <- paste(df_final_l4$ClusterRows, df_final_l4$group, sep='_')

######transfor the group as factor
df_final_l4$group1 <- df_final_l4$group1 %>% 
  factor()

######statistical assumptions
#1.interdependence of the variables

#2.variance homogeneity
bartlett.test(value ~ group1, data = df_final_l4)
#p-value<0.05. Reject the null hypothesis. The variance are different

#3.normality
qqnorm(df_final_l4$value)
qqline(df_final_l4$value, col='red')

###create the unique names for the for loop
mcluster <- df_final_l4$ClusterRows %>% unique

###Empty dataset
Res_letter_l4 <- NULL

###perform the for loop
for (cluster in mcluster){
  df_clust_l4 <- df_final_l4 %>% subset(ClusterRows == cluster) %>% droplevels
  #kruskal wallis
  m1_l4 <- kruskal.test(value ~ group1,data = df_clust_l4)
  ######
  DT_l4 = dunnTest(value ~ group1,data = df_clust_l4,method="bh")
  PT_l4 = DT_l4$res
  PT_l4
  df_letter_l4 <- cldList(P.adj ~ Comparison,
                          data = PT_l4, remove.zero = F,
                          threshold = 0.05)
  df_letter_l4$ClusterRows <- cluster
  Res_letter_l4 <- rbind(Res_letter_l4, df_letter_l4)
}

######calculate the mean and the confidence interval for each analysis
###create the unique names for the for loop
mgroup1 <- df_final_l4$group1 %>% unique
###Empty dataset
Res_mean_l4 <- NULL
Res_ci_l4 <- NULL

###perform the for loop
for (i in mgroup1){
  df_group1_l4 <- df_final_l4 %>% subset(group1 == i) %>% droplevels
  #calculate means
  df_mean_l4 <- mean(df_group1_l4$value) %>% as.data.frame
  ###calculate ci
  # bootstrapping with 1000 replications 
  results_l4 <- boot(data=df_group1_l4$value, statistic=Bmean, R=1000)
  # get 95% confidence interval 
  results_boot_l4 <- boot.ci(results_l4, type=c("norm", "basic", "perc", "bca"))
  ###create a data frame for ci
  df_ci_l4 <- results_boot_l4$bca %>% as.data.frame
  #add the leaf name on data frame
  df_mean_l4$group1 <- i
  df_ci_l4$group1 <- i
  #combine the results
  Res_mean_l4 <- rbind(Res_mean_l4, df_mean_l4)
  Res_ci_l4 <- rbind(Res_ci_l4, df_ci_l4)
}

######create the dataset formean
colnames(Res_mean_l4)[1]<-'mean'
###select column from 4to6
Res_ci_l4 <- Res_ci_l4[,4:6]
###rename the first column
colnames(Res_ci_l4)[1]<-'lower'
colnames(Res_ci_l4)[2]<-'upper'

######merge the two datasets
res_l4 <- merge(Res_mean_l4, Res_ci_l4, by='group1')

######remove res between group names
res_l4$group1 <- res_l4$group1 %>% gsub(pattern = ' ', replacement = '')

######split group1
res_l4$group1 <- res_l4$group1 %>% as.character
res_l4$ClusterRows <- lapply(strsplit(res_l4$group1, '_'), `[`,1) %>% as.character
res_l4$NumLeaf <- lapply(strsplit(res_l4$group1, '_'), `[`, 2) %>% as.character
res_l4$Treatment <- lapply(strsplit(res_l4$group1, '_'), `[`, 3) %>% as.character

######rename the information on the tReatment column
res_l4$Treatment[which(res_l4$Treatment == 'FullSynCom')] <- '+SynCom'
res_l4$Treatment[which(res_l4$Treatment == 'L3')] <- '+SynCom-L3'
res_l4$Treatment[which(res_l4$Treatment == 'L4')] <- '+SynCom-L4'
res_l4$Treatment[which(res_l4$Treatment == 'L5')] <- '+SynCom-L5'

######change the first column name 
colnames(Res_letter_l4)[1] <- 'group1'

######merge the dataset
res_l4 <- merge(res_l4, Res_letter_l4, by= c('group1', 'ClusterRows'))

######organize the CRs
df_final_l4$ClusterRows <- df_final_l4$ClusterRows %>%
  factor(levels = c('CR1', 'CR2', 'CR3', 'CR4', 'CR5',
                    'CR6', 'CR7', 'CR8'))

res_l4$ClusterRows <- res_l4$ClusterRows %>%
  factor(levels = c('CR1', 'CR2', 'CR3', 'CR4', 'CR5',
                    'CR6', 'CR7', 'CR8'))

######define the color for the dropout
paleta_syncom <- c('#d0d1e6','#a6bddb','#67a9cf',
                   '#3690c0','#02818a','#016450')
names(paleta_syncom) <- c('NB','+SynCom','+SynCom-L1&L2',
                          '+SynCom-L3','+SynCom-L4','+SynCom-L5')

######plot the result
p_l4 <- ggplot(data=df_final_l4, aes(x=Treatment, y=value, color=Treatment))+
  geom_point(shape=21, size=1.5, alpha=1)+
  facet_grid(.~ClusterRows, scales = 'free')+
  geom_pointrange(data=res_l4, aes(y=mean, ymin=lower, ymax=upper,
                                   color=Treatment), 
                  size=1)+
  geom_hline(yintercept=0, linetype="dashed", color = "black", size = 0.5)+
  geom_text(data = res_l4, aes(x = Treatment,y = 2.5,
                               label = Letter),
            inherit.aes = F,size = 6,family ="Arial",color = "black") +
  #scale_color_manual(values = paleta_syncom)+
  xlab(NULL)+ 
  ylab('Standardized expression')+
  scale_color_manual(values = paleta_syncom)+
  clean +
  theme(axis.text.x = element_text(size=12, angle=45, hjust=1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.position = 'none')

p_l4

######save the figure
oh.save.pdf(p = p_l4,
            outname = "figS4B_gene_expression_L4_all_clusters.pdf",
            outdir = "../figures/",width = 30,height = 10)

######select significant genes for l5
Res_sig_l5  <- Res_contrasts %>% 
  subset(padj < 0.1)  %>%
  dplyr::filter(.data = .,Contrast %>% grepl(pattern = "_L5_FullSynCom")) %>%
  droplevels

sig_genes_l5 <- Res_sig_l5$gene_id %>% as.character %>%
  unique

### Map to entrez id 
Tab_sub_l5 <- Tab_av[,colnames(Tab_av) %>% grep(pattern = "^L5")]

Tab_sub_l5 <- match(sig_genes_l5,rownames(Tab_sub_l5)) %>%
  Tab_sub_l5[.,]

######heatmap
res_heatmap_l5 <- chibi.heatmap(Tab = Tab_sub_l5,
                                   dist_method_rows = "euclidean",
                                   hclust_method_rows = "ward.D",
                                   range_fill_heatmap = c(-2,2),
                                   k_rows = 4,panel_spacing = 0.2,
                                   k_cols = 5,size_axis_text_row = 6,
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

d <- res_heatmap_l5$heatmap
d

######select the clusters
df_clust_rows_l5 <- res_heatmap_l5$df_clust_rows

######transform the Tab_sub_l3 in data frame
df_l5 <- Tab_sub_l5 %>% data.frame

######create the gene_id column
df_l5$gene_id <- row.names(df_l5)

######clean the rows
row.names(df_l5) <- NULL

########select the clusters
mlist_query_l5 <- list(
  CR1 = df_clust_rows_l5 %>% subset(ClusterRows == "CR1") %$% IdRows %>% as.character,
  CR2 = df_clust_rows_l5 %>% subset(ClusterRows == "CR2") %$% IdRows %>% as.character,
  CR3 = df_clust_rows_l5 %>% subset(ClusterRows == "CR3") %$% IdRows %>% as.character,
  CR4 = df_clust_rows_l5 %>% subset(ClusterRows == "CR4") %$% IdRows %>% as.character
)

######select the genes significant genes in L3
toselect_l5 <- mlist_query_l5 %>% unlist %>% as.character %>% unique

######select the cluster 4,6,8 for the different leaves
melted_all_l5 <- Tab_av %>% melt %>%
  dplyr::rename(.data = .,gene_id = Var1,group = Var2)

######select the genes
melted_sub_all_l5 <- melted_all_l5 %>%
  dplyr::filter(.data = .,gene_id %in% toselect_l5) %>% droplevels

######add the cluster information
melted_sub_all_l5$ClusterRows <- match(melted_sub_all_l5$gene_id,
                                       df_clust_rows_l5$IdRows) %>%
  df_clust_rows_l5$ClusterRows[.]

######create the group
melted_sub_all_l5$group <- melted_sub_all_l5$group %>% 
  factor(levels = c("L3_NB","L3_FullSynCom","L3_L3","L3_L4_up","L3_L5",
                    "L4_NB","L4_FullSynCom","L4_L3","L4_L4_up","L4_L5",
                    "L5_NB","L5_FullSynCom","L5_L3","L5_L4_up","L5_L5"))

melted_sub_all_l5$NumLeaf <- "L3"
melted_sub_all_l5$NumLeaf[melted_sub_all_l5$group %>% grep(pattern = "^L4")] <- "L4"
melted_sub_all_l5$NumLeaf[melted_sub_all_l5$group %>% grep(pattern = "^L5")] <- "L5"
melted_sub_all_l5$NumLeaf <- melted_sub_all_l5$NumLeaf %>% factor

melted_sub_all_l5$Treatment <- melted_sub_all_l5$group %>%
  gsub(pattern = "^L[345]_",replacement = "") %>% 
  factor(levels = c("NB","FullSynCom","L3","L4_up","L5"))

######rename the information on the tReatment column
melted_sub_all_l5$Treatment <- melted_sub_all_l5$Treatment %>% as.character
melted_sub_all_l5$Treatment[which(melted_sub_all_l5$Treatment == 'FullSynCom')] <- '+SynCom'
melted_sub_all_l5$Treatment[which(melted_sub_all_l5$Treatment == 'L3')] <- '+SynCom-L3'
melted_sub_all_l5$Treatment[which(melted_sub_all_l5$Treatment == 'L4_up')] <- '+SynCom-L4'
melted_sub_all_l5$Treatment[which(melted_sub_all_l5$Treatment == 'L5')] <- '+SynCom-L5'

######order the Treatment columns
melted_sub_all_l5$Treatment <- melted_sub_all_l5$Treatment %>% 
  factor(levels=c('NB', '+SynCom', '+SynCom-L3', '+SynCom-L4', '+SynCom-L5'))

#######create boxplot
ggplot(data=melted_sub_all_l5, mapping = aes(x=group,y=value))+
  geom_boxplot()+
  geom_point()+
  facet_grid(NumLeaf~ClusterRows)

######select only L3
df_final_l5 <- melted_sub_all_l5 %>% subset(NumLeaf == 'L5') %>% droplevels

######check the expression
ggplot(data=df_final_l5, mapping = aes(x=group,y=value))+
  geom_boxplot()+
  geom_point()+
  facet_grid(NumLeaf~ClusterRows)

#######create a group with group and ClusterRows
df_final_l5$group1 <- paste(df_final_l5$ClusterRows, df_final_l5$group, sep='_')

######transfor the group as factor
df_final_l5$group1 <- df_final_l5$group1 %>% 
  factor()

######statistical assumptions
#1.interdependence of the variables

#2.variance homogeneity
bartlett.test(value ~ group1, data = df_final_l5)
#p-value<0.05. Reject the null hypothesis. The variance are different

#3.normality
qqnorm(df_final_l5$value)
qqline(df_final_l5$value, col='red')

###create the unique names for the for loop
mcluster <- df_final_l5$ClusterRows %>% unique

###Empty dataset
Res_letter_l5 <- NULL

###perform the for loop
for (cluster in mcluster){
  df_clust_l5 <- df_final_l5 %>% subset(ClusterRows == cluster) %>% droplevels
  #kruskal wallis
  m1_l5 <- kruskal.test(value ~ group1,data = df_clust_l5)
  ######
  DT_l5 = dunnTest(value ~ group1,data = df_clust_l5,method="bh")
  PT_l5 = DT_l5$res
  PT_l5
  df_letter_l5 <- cldList(P.adj ~ Comparison,
                          data = PT_l5, remove.zero = F,
                          threshold = 0.05)
  df_letter_l5$ClusterRows <- cluster
  Res_letter_l5 <- rbind(Res_letter_l5, df_letter_l5)
}

######calculate the mean and the confidence interval for each analysis
###create the unique names for the for loop
mgroup1 <- df_final_l5$group1 %>% unique
###Empty dataset
Res_mean_l5 <- NULL
Res_ci_l5 <- NULL

###perform the for loop
for (i in mgroup1){
  df_group1_l5 <- df_final_l5 %>% subset(group1 == i) %>% droplevels
  #calculate means
  df_mean_l5 <- mean(df_group1_l5$value) %>% as.data.frame
  ###calculate ci
  # bootstrapping with 1000 replications 
  results_l5 <- boot(data=df_group1_l5$value, statistic=Bmean, R=1000)
  # get 95% confidence interval 
  results_boot_l5 <- boot.ci(results_l5, type=c("norm", "basic", "perc", "bca"))
  ###create a data frame for ci
  df_ci_l5 <- results_boot_l5$bca %>% as.data.frame
  #add the leaf name on data frame
  df_mean_l5$group1 <- i
  df_ci_l5$group1 <- i
  #combine the results
  Res_mean_l5 <- rbind(Res_mean_l5, df_mean_l5)
  Res_ci_l5 <- rbind(Res_ci_l5, df_ci_l5)
}

######create the dataset formean
colnames(Res_mean_l5)[1]<-'mean'
###select column from 4to6
Res_ci_l5 <- Res_ci_l5[,4:6]
###rename the first column
colnames(Res_ci_l5)[1]<-'lower'
colnames(Res_ci_l5)[2]<-'upper'

######merge the two datasets
res_l5 <- merge(Res_mean_l5, Res_ci_l5, by='group1')

######remove res between group names
res_l5$group1 <- res_l5$group1 %>% gsub(pattern = ' ', replacement = '')

######split group1
res_l5$group1 <- res_l5$group1 %>% as.character
res_l5$ClusterRows <- lapply(strsplit(res_l5$group1, '_'), `[`,1) %>% as.character
res_l5$NumLeaf <- lapply(strsplit(res_l5$group1, '_'), `[`, 2) %>% as.character
res_l5$Treatment <- lapply(strsplit(res_l5$group1, '_'), `[`, 3) %>% as.character

######rename the information on the tReatment column
res_l5$Treatment[which(res_l5$Treatment == 'FullSynCom')] <- '+SynCom'
res_l5$Treatment[which(res_l5$Treatment == 'L3')] <- '+SynCom-L3'
res_l5$Treatment[which(res_l5$Treatment == 'L4')] <- '+SynCom-L4'
res_l5$Treatment[which(res_l5$Treatment == 'L5')] <- '+SynCom-L5'

######change the first column name 
colnames(Res_letter_l5)[1] <- 'group1'

######merge the dataset
res_l5 <- merge(res_l5, Res_letter_l5, by= c('group1', 'ClusterRows'))

######organize the CRs
df_final_l5$ClusterRows <- df_final_l5$ClusterRows %>%
  factor(levels = c('CR1', 'CR2', 'CR3', 'CR4'))

res_l5$ClusterRows <- res_l5$ClusterRows %>%
  factor(levels = c('CR1', 'CR2', 'CR3', 'CR4'))

######define the color for the dropout
paleta_syncom <- c('#d0d1e6','#a6bddb','#67a9cf',
                   '#3690c0','#02818a','#016450')
names(paleta_syncom) <- c('NB','+SynCom','+SynCom-L1&L2',
                          '+SynCom-L3','+SynCom-L4','+SynCom-L5')

######plot the result
p_l5 <- ggplot(data=df_final_l5, aes(x=Treatment, y=value, color=Treatment))+
  geom_point(shape=21, size=1.5, alpha=1)+
  facet_grid(.~ClusterRows, scales = 'free')+
  geom_pointrange(data=res_l5, aes(y=mean, ymin=lower, ymax=upper,
                                   color=Treatment), 
                  size=1)+
  geom_hline(yintercept=0, linetype="dashed", color = "black", size = 0.5)+
  geom_text(data = res_l5, aes(x = Treatment,y = 2.5,
                               label = Letter),
            inherit.aes = F,size = 6,family ="Arial",color = "black") +
  #scale_color_manual(values = paleta_syncom)+
  xlab(NULL)+ 
  ylab('Standardized expression')+
  scale_color_manual(values = paleta_syncom)+
  clean +
  theme(axis.text.x = element_text(size=12, angle=45, hjust=1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.position = 'none')

p_l5

######save the figure
oh.save.pdf(p = p_l5,
            outname = "figS4k_gene_expression_L5_all_clusters.pdf",
            outdir = "../figures/",width = 20,height = 10)
