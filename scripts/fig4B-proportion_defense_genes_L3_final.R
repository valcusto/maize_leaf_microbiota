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
df_gene <- read.csv('../cleandata/defense_genes_in_CR468_CR8C.csv')

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

######prepare the dataset
melted_sub$group <- melted_sub$group %>% 
  factor(levels = c("L3_NB","L3_FullSynCom","L3_L3","L3_L4_up","L3_L5"))

######create the NumLeaf column
melted_sub$NumLeaf <- "L3"

######create the Treatment column
melted_sub$Treatment <- melted_sub$group %>%
  gsub(pattern = "^L[345]_",replacement = "") %>% as.character

######merge the melted and df_clustgene
melted_sub <- merge(melted_sub, df_clustgene, by='gene_id')

######save the melted_sum
write.table(x=melted_sub, file='../cleandata/subclusters_L3.csv',
            sep=',', row.names = F, col.names = T)

######rename the information on the tReatment column
melted_sub$Treatment[which(melted_sub$Treatment == 'FullSynCom')] <- '+SynCom'
melted_sub$Treatment[which(melted_sub$Treatment == 'L3')] <- '+SynCom-L3'
melted_sub$Treatment[which(melted_sub$Treatment == 'L4_up')] <- '+SynCom-L4'
melted_sub$Treatment[which(melted_sub$Treatment == 'L5')] <- '+SynCom-L5'

######order the Treatment columns
melted_sub$Treatment <- melted_sub$Treatment %>% 
  factor(levels=c('NB', '+SynCom', '+SynCom-L3','+SynCom-L4','+SynCom-L5'))


######how many genes we have per cluster
cr1 <- melted_sub %>%
  subset(Clustgene == 'CR1') %>% droplevels
total_cr1 <- cr1$gene_id %>% unique %>% length

cr2 <- melted_sub %>%
  subset(Clustgene == 'CR2') %>% droplevels
total_cr2 <- cr2$gene_id %>% unique %>% length

cr3 <- melted_sub %>%
  subset(Clustgene == 'CR3') %>% droplevels
total_cr3 <- cr3$gene_id %>% unique %>% length

cr4 <- melted_sub %>%
  subset(Clustgene == 'CR4') %>% droplevels
total_cr4 <- cr4$gene_id %>% unique %>% length

######find the defense genes 
df_defense <- df_gene %>%
  subset(category == 'defense') %>% droplevels

######find the unique defense genes
unique_defense <- df_defense$IdRows %>% unique

######count the number of defense genes per cluster
melted_defense <- melted_sub[which(melted_sub$gene_id %in% unique_defense),]

#######find the number of defense genes per cluster
cr1_defense <- melted_defense %>%
  subset(Clustgene == 'CR1') %>% droplevels
defense_cr1 <- cr1_defense$gene_id %>% unique %>% length

cr2_defense <- melted_defense %>%
  subset(Clustgene == 'CR2') %>% droplevels
defense_cr2 <- cr2_defense$gene_id %>% unique %>% length

cr3_defense <- melted_defense %>%
  subset(Clustgene == 'CR3') %>% droplevels
defense_cr3 <- cr3_defense$gene_id %>% unique %>% length

cr4_defense <- melted_defense %>%
  subset(Clustgene == 'CR4') %>% droplevels
defense_cr4 <- cr4_defense$gene_id %>% unique %>% length

######create a dataframe
df_final <- data_frame(cluster = c('genome', 'CR1', 'CR2', 'CR3', 'CR4'),
                       total = c(39756, total_cr1, total_cr2, total_cr3, total_cr4),
                       defense = c(4914, defense_cr1, defense_cr2, defense_cr3,
                                   defense_cr4))

######calculate the proportion of defense genes per cluster
df_final$proportion <- df_final$defense/df_final$total

######calculate the hyprogeometric analysis
fcr1 <- phyper(defense_cr1-1, total_cr1, 39756-total_cr1,4914,
               lower.tail = FALSE)
fcr2 <- phyper(defense_cr2-1, total_cr2, 39756-total_cr2,4914,
               lower.tail = FALSE)
fcr3 <- phyper(defense_cr3-1, total_cr3, 39756-total_cr3,4914,
               lower.tail = FALSE)
fcr4 <- phyper(defense_cr4-1, total_cr4, 39756-total_cr4,4914,
               lower.tail = FALSE)

######add the pvalue
df_final$pvalue <- c(1, fcr1, fcr2,fcr3,fcr4)

######adjust the pvalue
df_final$padjust <- p.adjust(df_final$pvalue, method = 'fdr')

######add the significance column
df_final$significance <- 'NS'
df_final$significance[which(df_final$pvalue < 0.05)] <- 'p < 0.05'

######order the factors
df_final$cluster <- df_final$cluster %>%
  factor(levels = c('genome', 'CR1', 'CR2', 'CR3', 'CR4'))

######create the plot
p <- ggplot(df_final, aes(x=cluster, y=proportion, fill = significance))+
  geom_col()+
  geom_hline(yintercept = 0.1236040, linetype = 'dashed', color = 'black',
             size=0.3)+
  ylab('Proportion of defence genes (%)')+
  scale_y_continuous(labels = scales::percent, limits=c(0,0.5))+
  scale_fill_manual(name = 'Significance',labels = c('NS',
                                                     'p < 0.05'),
                    values=c('#999999', '#d8daeb'))+
  clean+
  theme(axis.text.x = element_text(size=18,vjust=1, angle=45, hjust=1),
        axis.text.y = element_text(size=18),
        legend.key.height= unit(0.3, 'cm'),
        legend.key.width = unit(0.6, 'cm'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        legend.position = 'right')

p

######save figures
oh.save.pdf(p = p,outname = "fig4b_defence_genes_proportion_fourclusters_L3_final.pdf",
            outdir = "../figures/",width = 6,height = 5)
