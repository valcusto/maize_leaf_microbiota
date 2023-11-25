######load packages
library(tidyverse)
library(ohchibi)
library(paletteer)
library(scales)
library(dplyr)
library(ggtree)
library(clusterProfiler) #GO analysis
library(biomaRt)

######set seed
set.seed(130816)

######set work directory
setwd("C:/Users/Castrillo Lab/OneDrive - The University of Nottingham/Publication 2019/Paper preparation/manuscripts/scripts")

######helper functions
source('0-Clean_up_plots.R')
source('0_5-heatmap_ohchibi.R')
source('0_6-theme_ohchibi.R')
source('0_7-chibi.heatmap.R')

split_and_repeat <- function(data, col_name) {
  split_values <- strsplit(data[[col_name]], "/")
  rows <- lengths(split_values)
  data <- data[rep(seq_len(nrow(data)), rows), ]
  data[[col_name]] <- unlist(split_values)
  rownames(data) <- NULL
  return(data)
}

### Correlation data
df_cor <- read.table(file = "../rawdata/res_cor_geneexp_length.tsv")

######open the expression dataset
Dat_av  <- readRDS(file = "../cleandata/res_tab_av_rnaseq.RDS")
Tab_av <- Dat_av$Tab_av

##########################Step01:Data preparation
######rename the columns
colnames(df_cor) <- c('IdRows', 'NumLeaf', 'r_pearson', 'pvalue')

######open the gene backgroun
term2gene <- read.csv("../rawdata/Zmays.GO.anno_rice_At.term_gene.txt", header=F, sep="\t")
term2name <- read.csv("../rawdata/Zmays.GO.anno_rice_At.term_name.txt", header=F, sep="\t")

######select L3
df_l3 <- df_cor %>% subset(pvalue < 0.01) %>% subset(NumLeaf == "L3")
up_l3 <- df_l3 %>% subset(r_pearson > 0) %$% IdRows
down_l3 <- df_l3 %>% subset(r_pearson < 0) %$% IdRows

######select the significants for L3
sig_genes_l3 <- df_l3$IdRows %>% as.character %>%
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

######select only the signifcant ones
sig_genes_l3 <- df_clust_rows_l3$IdRows

######how many genes significant
length(sig_genes_l3) #459 genes

######select each cluster for L3
cr1c_l3 <- df_clust_rows_l3 %>%
  subset(ClusterRows == 'CR1C') %>% droplevels
cr2c_l3 <- df_clust_rows_l3 %>%
  subset(ClusterRows == 'CR2C') %>% droplevels
cr3c_l3 <- df_clust_rows_l3 %>%
  subset(ClusterRows == 'CR3C') %>% droplevels
cr4c_l3 <- df_clust_rows_l3 %>%
  subset(ClusterRows == 'CR4C') %>% droplevels
cr5c_l3 <- df_clust_rows_l3 %>%
  subset(ClusterRows == 'CR5C') %>% droplevels
cr6c_l3 <- df_clust_rows_l3 %>%
  subset(ClusterRows == 'CR6C') %>% droplevels
cr7c_l3 <- df_clust_rows_l3 %>%
  subset(ClusterRows == 'CR7C') %>% droplevels
cr8c_l3 <- df_clust_rows_l3 %>%
  subset(ClusterRows == 'CR8C') %>% droplevels
cr9c_l3 <- df_clust_rows_l3 %>%
  subset(ClusterRows == 'CR9C') %>% droplevels

######cr1
gene_cr1c_l3 <- as.factor(cr1c_l3$IdRows)

######GO enrichment analysis using each clusters list
x_cr1c_l3 <- enricher(gene_cr1c_l3, TERM2GENE=term2gene, TERM2NAME=term2name,
                     pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)

######What is the percentage of genes with GO?
df_cr1c_l3 <- split_and_repeat(x_cr1c_l3@result, "geneID")
unique_cr1c_l3 <- df_cr1c_l3$geneID %>% unique
###answer:
cr1c_l3_total <- length(gene_cr1c_l3)
cr1c_l3_go <- length(unique_cr1c_l3)
prop_go_cr1c_l3 <- cr1c_l3_go/cr1c_l3_total #76% GO identified

######cr2c
gene_cr2c_l3 <- as.factor(cr2c_l3$IdRows)

x_cr2c_l3 <- enricher(gene_cr2c_l3, TERM2GENE=term2gene, TERM2NAME=term2name,
                      pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)

######What is the percentage of genes with GO?
df_cr2c_l3 <- split_and_repeat(x_cr2c_l3@result, "geneID")
unique_cr2c_l3 <- df_cr2c_l3$geneID %>% unique
###answer:
cr2c_l3_total <- length(gene_cr2c_l3)
cr2c_l3_go <- length(unique_cr2c_l3)
prop_go_cr2c_l3 <- cr2c_l3_go/cr2c_l3_total #56% GO identified

######cr3c
gene_cr3c_l3 <- as.factor(cr3c_l3$IdRows)

x_cr3c_l3 <- enricher(gene_cr3c_l3, TERM2GENE=term2gene, TERM2NAME=term2name,
                      pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)

######What is the percentage of genes with GO?
df_cr3c_l3 <- split_and_repeat(x_cr3c_l3@result, "geneID")
unique_cr3c_l3 <- df_cr3c_l3$geneID %>% unique
###answer:
cr3c_l3_total <- length(gene_cr3c_l3)
cr3c_l3_go <- length(unique_cr3c_l3)
prop_go_cr3c_l3 <- cr3c_l3_go/cr3c_l3_total #94% GO identified

######cr4c
gene_cr4c_l3 <- as.factor(cr4c_l3$IdRows)

x_cr4c_l3 <- enricher(gene_cr4c_l3, TERM2GENE=term2gene, TERM2NAME=term2name,
                      pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)

######What is the percentage of genes with GO?
df_cr4c_l3 <- split_and_repeat(x_cr4c_l3@result, "geneID")
unique_cr4c_l3 <- df_cr4c_l3$geneID %>% unique
###answer:
cr4c_l3_total <- length(gene_cr4c_l3)
cr4c_l3_go <- length(unique_cr4c_l3)
prop_go_cr4c_l3 <- cr4c_l3_go/cr4c_l3_total #63% GO identified

######cr5c
gene_cr5c_l3 <- as.factor(cr5c_l3$IdRows)

x_cr5c_l3 <- enricher(gene_cr5c_l3, TERM2GENE=term2gene, TERM2NAME=term2name,
                      pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)

######What is the percentage of genes with GO?
df_cr5c_l3 <- split_and_repeat(x_cr5c_l3@result, "geneID")
unique_cr5c_l3 <- df_cr5c_l3$geneID %>% unique
###answer:
cr5c_l3_total <- length(gene_cr5c_l3)
cr5c_l3_go <- length(unique_cr5c_l3)
prop_go_cr5c_l3 <- cr5c_l3_go/cr5c_l3_total #62% GO identified

######cr6c
gene_cr6c_l3 <- as.factor(cr6c_l3$IdRows)

x_cr6c_l3 <- enricher(gene_cr6c_l3, TERM2GENE=term2gene, TERM2NAME=term2name,
                      pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)

######What is the percentage of genes with GO?
df_cr6c_l3 <- split_and_repeat(x_cr6c_l3@result, "geneID")
unique_cr6c_l3 <- df_cr6c_l3$geneID %>% unique
###answer:
cr6c_l3_total <- length(gene_cr6c_l3)
cr6c_l3_go <- length(unique_cr6c_l3)
prop_go_cr6c_l3 <- cr6c_l3_go/cr6c_l3_total #48% GO identified

######cr7c
gene_cr7c_l3 <- as.factor(cr7c_l3$IdRows)

x_cr7c_l3 <- enricher(gene_cr7c_l3, TERM2GENE=term2gene, TERM2NAME=term2name,
                      pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)

######What is the percentage of genes with GO?
df_cr7c_l3 <- split_and_repeat(x_cr7c_l3@result, "geneID")
unique_cr7c_l3 <- df_cr7c_l3$geneID %>% unique
###answer:
cr7c_l3_total <- length(gene_cr7c_l3)
cr7c_l3_go <- length(unique_cr7c_l3)
prop_go_cr7c_l3 <- cr7c_l3_go/cr7c_l3_total #31% GO identified

######cr8c
gene_cr8c_l3 <- as.factor(cr8c_l3$IdRows)

x_cr8c_l3 <- enricher(gene_cr8c_l3, TERM2GENE=term2gene, TERM2NAME=term2name,
                      pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)

######What is the percentage of genes with GO?
df_cr8c_l3 <- split_and_repeat(x_cr8c_l3@result, "geneID")
unique_cr8c_l3 <- df_cr8c_l3$geneID %>% unique
###answer:
cr8c_l3_total <- length(gene_cr8c_l3)
cr8c_l3_go <- length(unique_cr8c_l3)
prop_go_cr8c_l3 <- cr8c_l3_go/cr8c_l3_total #75% GO identified

######cr9c
gene_cr9c_l3 <- as.factor(cr9c_l3$IdRows)

x_cr9c_l3 <- enricher(gene_cr9c_l3, TERM2GENE=term2gene, TERM2NAME=term2name,
                      pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)

######What is the percentage of genes with GO?
df_cr9c_l3 <- split_and_repeat(x_cr9c_l3@result, "geneID")
unique_cr9c_l3 <- df_cr9c_l3$geneID %>% unique
###answer:
cr9c_l3_total <- length(gene_cr9c_l3)
cr9c_l3_go <- length(unique_cr9c_l3)
prop_go_cr9c_l3 <- cr9c_l3_go/cr9c_l3_total #66% GO identified

######merge results
df_results_L3 <- merge_result(list(CR9C = x_cr9c_l3,
                                   CR7C = x_cr7c_l3,
                                   CR4C = x_cr4c_l3,
                                   CR6C = x_cr6c_l3,
                                   CR8C = x_cr8c_l3,
                                   CR3C = x_cr3c_l3,
                                   CR2C = x_cr2c_l3,
                                   CR1C = x_cr1c_l3,
                                   CR5C = x_cr5c_l3))

df_results_l3_final <- merge_result(list(CR1C = x_cr1c_l3,
                                         CR2C = x_cr2c_l3,
                                         CR3C = x_cr3c_l3,
                                         CR4C = x_cr4c_l3,
                                         CR5C = x_cr5c_l3,
                                         CR6C = x_cr6c_l3,
                                         CR7C = x_cr7c_l3,
                                         CR8C = x_cr8c_l3,
                                         CR9C = x_cr9c_l3))

######save the dataframe
df_results_l3_1 <- df_results_l3_final@compareClusterResult

#write.table(df_results_l3_1, "../cleandata/df_results_correlation_GO_Jia_final_L3.txt", quote=F, 
            #row.names=F, col.names=T, sep="\t")

######doplot
p_L3 <- dotplot(df_results_L3, showCategory=5, title = "L3", 
                includeAll=FALSE)+
  scale_color_paletteer_c(name = 'q-value',"viridis::plasma",
                          na.value = "#BFBFBF")

p_L3 <- p_L3+
  ggtitle('L3')+
  clean+
  theme(plot.title = element_text(size=30, face = 'bold', hjust=0.5),
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=12))

p_L3

######save the L3 GO
oh.save.pdf(p = p_L3,
            outname = "figS4G_dotplot_GO_cor_L3_jia_final.pdf",
            outdir = '../figures/',
            width = 20,height = 25)


######select L4
df_l4 <- df_cor %>% subset(pvalue < 0.01) %>% subset(NumLeaf == "L4")
up_l4 <- df_l4 %>% subset(r_pearson > 0) %$% IdRows
down_l4 <- df_l4 %>% subset(r_pearson < 0) %$% IdRows

######select the significants for L3
sig_genes_l4 <- df_l4$IdRows %>% as.character %>%
  unique

### Map to entrez id 
Tab_sub_l4 <- Tab_av[,colnames(Tab_av) %>% grep(pattern = "^L4")]

Tab_sub_l4 <- match(sig_genes_l4,rownames(Tab_sub_l4)) %>%
  Tab_sub_l4[.,]

######rename the Tab_sub_l3 column names
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
df_clust_rows_l4 <- mclust_rows_l4 %>% cutree(k = 6) %>%
  data.frame(IdRows = names(.), ClusterRows = paste0("CR",.,"C"),row.names = NULL)

######remove the first column
df_clust_rows_l4 <- df_clust_rows_l4[,-1]

######select only the signifcant ones
sig_genes_l4 <- df_clust_rows_l4$IdRows

######how many genes significant
length(sig_genes_l4) #122 genes

######select each cluster for L3
cr1c_l4 <- df_clust_rows_l4 %>%
  subset(ClusterRows == 'CR1C') %>% droplevels
cr2c_l4 <- df_clust_rows_l4 %>%
  subset(ClusterRows == 'CR2C') %>% droplevels
cr3c_l4 <- df_clust_rows_l4 %>%
  subset(ClusterRows == 'CR3C') %>% droplevels
cr4c_l4 <- df_clust_rows_l4 %>%
  subset(ClusterRows == 'CR4C') %>% droplevels
cr5c_l4 <- df_clust_rows_l4 %>%
  subset(ClusterRows == 'CR5C') %>% droplevels
cr6c_l4 <- df_clust_rows_l4 %>%
  subset(ClusterRows == 'CR6C') %>% droplevels

######cr1
gene_cr1c_l4 <- as.factor(cr1c_l4$IdRows)

######GO enrichment analysis using each clusters list
x_cr1c_l4 <- enricher(gene_cr1c_l4, TERM2GENE=term2gene, TERM2NAME=term2name,
                      pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)

######What is the percentage of genes with GO?
df_cr1c_l4 <- split_and_repeat(x_cr1c_l4@result, "geneID")
unique_cr1c_l4 <- df_cr1c_l4$geneID %>% unique
###answer:
cr1c_l4_total <- length(gene_cr1c_l4)
cr1c_l4_go <- length(unique_cr1c_l4)
prop_go_cr1c_l4 <- cr1c_l4_go/cr1c_l4_total #53% GO identified

######cr2c
gene_cr2c_l4 <- as.factor(cr2c_l4$IdRows)

x_cr2c_l4 <- enricher(gene_cr2c_l4, TERM2GENE=term2gene, TERM2NAME=term2name,
                      pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)

######What is the percentage of genes with GO?
df_cr2c_l4 <- split_and_repeat(x_cr2c_l4@result, "geneID")
unique_cr2c_l4 <- df_cr2c_l4$geneID %>% unique
###answer:
cr2c_l4_total <- length(gene_cr2c_l4)
cr2c_l4_go <- length(unique_cr2c_l4)
prop_go_cr2c_l4 <- cr2c_l4_go/cr2c_l4_total #52% GO identified

######cr3c
gene_cr3c_l4 <- as.factor(cr3c_l4$IdRows)

x_cr3c_l4 <- enricher(gene_cr3c_l4, TERM2GENE=term2gene, TERM2NAME=term2name,
                      pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)

######What is the percentage of genes with GO?
df_cr3c_l4 <- split_and_repeat(x_cr3c_l4@result, "geneID")
unique_cr3c_l4 <- df_cr3c_l4$geneID %>% unique
###answer:
cr3c_l4_total <- length(gene_cr3c_l4)
cr3c_l4_go <- length(unique_cr3c_l4)
prop_go_cr3c_l4 <- cr3c_l4_go/cr3c_l4_total #87% GO identified

######cr4c
gene_cr4c_l4 <- as.factor(cr4c_l4$IdRows)

x_cr4c_l4 <- enricher(gene_cr4c_l4, TERM2GENE=term2gene, TERM2NAME=term2name,
                      pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)

######What is the percentage of genes with GO?
df_cr4c_l4 <- split_and_repeat(x_cr4c_l4@result, "geneID")
unique_cr4c_l4 <- df_cr4c_l4$geneID %>% unique
###answer:
cr4c_l4_total <- length(gene_cr4c_l4)
cr4c_l4_go <- length(unique_cr4c_l4)
prop_go_cr4c_l4 <- cr4c_l4_go/cr4c_l4_total #70% GO identified

######cr5c
gene_cr5c_l4 <- as.factor(cr5c_l4$IdRows)

x_cr5c_l4 <- enricher(gene_cr5c_l4, TERM2GENE=term2gene, TERM2NAME=term2name,
                      pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)

######What is the percentage of genes with GO?
df_cr5c_l4 <- split_and_repeat(x_cr5c_l4@result, "geneID")
unique_cr5c_l4 <- df_cr5c_l4$geneID %>% unique
###answer:
cr5c_l4_total <- length(gene_cr5c_l4)
cr5c_l4_go <- length(unique_cr5c_l4)
prop_go_cr5c_l4 <- cr5c_l4_go/cr5c_l4_total #47% GO identified

######cr6c
gene_cr6c_l4 <- as.factor(cr6c_l4$IdRows)

x_cr6c_l4 <- enricher(gene_cr6c_l4, TERM2GENE=term2gene, TERM2NAME=term2name,
                      pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)

######What is the percentage of genes with GO?
df_cr6c_l4 <- split_and_repeat(x_cr6c_l4@result, "geneID")
unique_cr6c_l4 <- df_cr6c_l4$geneID %>% unique
###answer:
cr6c_l4_total <- length(gene_cr6c_l4)
cr6c_l4_go <- length(unique_cr6c_l4)
prop_go_cr6c_l4 <- cr6c_l4_go/cr6c_l4_total #100% GO identified

######merge results
df_results_l4 <- merge_result(list(
  CR4C = x_cr4c_l4,
  CR6C = x_cr6c_l4,
  CR3C = x_cr3c_l4,
  CR2C = x_cr2c_l4,
  CR1C = x_cr1c_l4,
  CR5C = x_cr5c_l4))

df_results_l4_final <- merge_result(list(CR1C = x_cr1c_l4,
                                         CR2C = x_cr2c_l4,
                                         CR3C = x_cr3c_l4,
                                         CR4C = x_cr4c_l4,
                                         CR5C = x_cr5c_l4,
                                         CR6C = x_cr6c_l4))

######save the dataframe
df_results_l4_1 <- df_results_l4_final@compareClusterResult

#write.table(df_results_l4_1, "../cleandata/df_results_correlation_GO_Jia_final_L4.txt", quote=F, 
            #row.names=F, col.names=T, sep="\t")

######doplot
p_l4 <- dotplot(df_results_l4, showCategory=5, title = "L4", 
                includeAll=FALSE)+
  scale_color_paletteer_c(name = 'q-value',"viridis::plasma",
                          na.value = "#BFBFBF")

p_l4 <- p_l4+
  ggtitle('L4')+
  clean+
  theme(plot.title = element_text(size=30, face = 'bold', hjust=0.5),
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=12))

p_l4

######save the L4 GO
oh.save.pdf(p = p_l4,
            outname = "figS4G_dotplot_GO_cor_L4_jia_final.pdf",
            outdir = '../figures/',
            width = 20,height = 25)

######select L5
df_l5 <- df_cor %>% subset(pvalue < 0.01) %>% subset(NumLeaf == "L5")
up_l5 <- df_l5 %>% subset(r_pearson > 0) %$% IdRows
down_l5 <- df_l5 %>% subset(r_pearson < 0) %$% IdRows

######select the significants for L3
sig_genes_l5 <- df_l5$IdRows %>% as.character %>%
  unique

### Map to entrez id 
Tab_sub_l5 <- Tab_av[,colnames(Tab_av) %>% grep(pattern = "^L5")]

Tab_sub_l5 <- match(sig_genes_l5,rownames(Tab_sub_l5)) %>%
  Tab_sub_l5[.,]

######rename the Tab_sub_l3 column names
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
df_clust_rows_l5 <- mclust_rows_l5 %>% cutree(k = 7) %>%
  data.frame(IdRows = names(.), ClusterRows = paste0("CR",.,"C"),row.names = NULL)

######remove the first column
df_clust_rows_l5 <- df_clust_rows_l5[,-1]

######select only the signifcant ones
sig_genes_l5 <- df_clust_rows_l5$IdRows

######how many genes significant
length(sig_genes_l5) #2337 genes

######select each cluster for L3
cr1c_l5 <- df_clust_rows_l5 %>%
  subset(ClusterRows == 'CR1C') %>% droplevels
cr2c_l5 <- df_clust_rows_l5 %>%
  subset(ClusterRows == 'CR2C') %>% droplevels
cr3c_l5 <- df_clust_rows_l5 %>%
  subset(ClusterRows == 'CR3C') %>% droplevels
cr4c_l5 <- df_clust_rows_l5 %>%
  subset(ClusterRows == 'CR4C') %>% droplevels
cr5c_l5 <- df_clust_rows_l5 %>%
  subset(ClusterRows == 'CR5C') %>% droplevels
cr6c_l5 <- df_clust_rows_l5 %>%
  subset(ClusterRows == 'CR6C') %>% droplevels
cr7c_l5 <- df_clust_rows_l5 %>%
  subset(ClusterRows == 'CR7C') %>% droplevels

######cr1
gene_cr1c_l5 <- as.factor(cr1c_l5$IdRows)

######GO enrichment analysis using each clusters list
x_cr1c_l5 <- enricher(gene_cr1c_l5, TERM2GENE=term2gene, TERM2NAME=term2name,
                      pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)

######What is the percentage of genes with GO?
df_cr1c_l5 <- split_and_repeat(x_cr1c_l5@result, "geneID")
unique_cr1c_l5 <- df_cr1c_l5$geneID %>% unique
###answer:
cr1c_l5_total <- length(gene_cr1c_l5)
cr1c_l5_go <- length(unique_cr1c_l5)
prop_go_cr1c_l5 <- cr1c_l5_go/cr1c_l5_total #64% GO identified

######cr2c
gene_cr2c_l5 <- as.factor(cr2c_l5$IdRows)

x_cr2c_l5 <- enricher(gene_cr2c_l5, TERM2GENE=term2gene, TERM2NAME=term2name,
                      pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)

######What is the percentage of genes with GO?
df_cr2c_l5 <- split_and_repeat(x_cr2c_l5@result, "geneID")
unique_cr2c_l5 <- df_cr2c_l5$geneID %>% unique
###answer:
cr2c_l5_total <- length(gene_cr2c_l5)
cr2c_l5_go <- length(unique_cr2c_l5)
prop_go_cr2c_l5 <- cr2c_l5_go/cr2c_l5_total #64% GO identified

######cr3c
gene_cr3c_l5 <- as.factor(cr3c_l5$IdRows)

x_cr3c_l5 <- enricher(gene_cr3c_l5, TERM2GENE=term2gene, TERM2NAME=term2name,
                      pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)

######What is the percentage of genes with GO?
df_cr3c_l5 <- split_and_repeat(x_cr3c_l5@result, "geneID")
unique_cr3c_l5 <- df_cr3c_l5$geneID %>% unique
###answer:
cr3c_l5_total <- length(gene_cr3c_l5)
cr3c_l5_go <- length(unique_cr3c_l5)
prop_go_cr3c_l5 <- cr3c_l5_go/cr3c_l5_total #69% GO identified

######cr4c
gene_cr4c_l5 <- as.factor(cr4c_l5$IdRows)

x_cr4c_l5 <- enricher(gene_cr4c_l5, TERM2GENE=term2gene, TERM2NAME=term2name,
                      pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)

######What is the percentage of genes with GO?
df_cr4c_l5 <- split_and_repeat(x_cr4c_l5@result, "geneID")
unique_cr4c_l5 <- df_cr4c_l5$geneID %>% unique
###answer:
cr4c_l5_total <- length(gene_cr4c_l5)
cr4c_l5_go <- length(unique_cr4c_l5)
prop_go_cr4c_l5 <- cr4c_l5_go/cr4c_l5_total #46% GO identified

######cr5c
gene_cr5c_l5 <- as.factor(cr5c_l5$IdRows)

x_cr5c_l5 <- enricher(gene_cr5c_l5, TERM2GENE=term2gene, TERM2NAME=term2name,
                      pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)

######What is the percentage of genes with GO?
df_cr5c_l5 <- split_and_repeat(x_cr5c_l5@result, "geneID")
unique_cr5c_l5 <- df_cr5c_l5$geneID %>% unique
###answer:
cr5c_l5_total <- length(gene_cr5c_l5)
cr5c_l5_go <- length(unique_cr5c_l5)
prop_go_cr5c_l5 <- cr5c_l5_go/cr5c_l5_total #67% GO identified

######cr6c
gene_cr6c_l5 <- as.factor(cr6c_l5$IdRows)

x_cr6c_l5 <- enricher(gene_cr6c_l5, TERM2GENE=term2gene, TERM2NAME=term2name,
                      pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)

######What is the percentage of genes with GO?
df_cr6c_l5 <- split_and_repeat(x_cr6c_l5@result, "geneID")
unique_cr6c_l5 <- df_cr6c_l5$geneID %>% unique
###answer:
cr6c_l5_total <- length(gene_cr6c_l5)
cr6c_l5_go <- length(unique_cr6c_l5)
prop_go_cr6c_l5 <- cr6c_l5_go/cr6c_l5_total #63% GO identified

######cr7c
gene_cr7c_l5 <- as.factor(cr7c_l5$IdRows)

x_cr7c_l5 <- enricher(gene_cr7c_l5, TERM2GENE=term2gene, TERM2NAME=term2name,
                      pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)

######What is the percentage of genes with GO?
df_cr7c_l5 <- split_and_repeat(x_cr7c_l5@result, "geneID")
unique_cr7c_l5 <- df_cr7c_l5$geneID %>% unique
###answer:
cr7c_l5_total <- length(gene_cr7c_l5)
cr7c_l5_go <- length(unique_cr7c_l5)
prop_go_cr7c_l5 <- cr7c_l5_go/cr7c_l5_total #73% GO identified

######merge results
df_results_l5 <- merge_result(list(CR7C = x_cr7c_l5,
                                   CR4C = x_cr4c_l5,
                                   CR6C = x_cr6c_l5,
                                   CR3C = x_cr3c_l5,
                                   CR2C = x_cr2c_l5,
                                   CR1C = x_cr1c_l5,
                                   CR5C = x_cr5c_l5))

df_results_l5_final <- merge_result(list(CR1C = x_cr1c_l5,
                                         CR2C = x_cr2c_l5,
                                         CR3C = x_cr3c_l5,
                                         CR4C = x_cr4c_l5,
                                         CR5C = x_cr5c_l5,
                                         CR6C = x_cr6c_l5,
                                         CR7C = x_cr7c_l5))

######save the dataframe
df_results_l5_1 <- df_results_l5_final@compareClusterResult

#write.table(df_results_l5_1, "../cleandata/df_results_correlation_GO_Jia_final_L5.txt", quote=F, 
            #row.names=F, col.names=T, sep="\t")

######doplot
p_l5 <- dotplot(df_results_l5, showCategory=5, title = "L5", 
                includeAll=FALSE)+
  scale_color_paletteer_c(name = 'q-value',"viridis::plasma",
                          na.value = "#BFBFBF")

p_l5 <- p_l5+
  ggtitle('L5')+
  clean+
  theme(plot.title = element_text(size=30, face = 'bold', hjust=0.5),
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=12))

p_l5

######save the l5 GO
oh.save.pdf(p = p_l5,
            outname = "figS4G_dotplot_GO_cor_L5_jia_final.pdf",
            outdir = '../figures/',
            width = 20,height = 25)


######composition
composition = egg::ggarrange(p_L3, p_l4, p_l5,
                             nrow=1)

######save the plot as pdf
oh.save.pdf(p = composition,
            outname = "figS4G_dotplot_GO_cor_geneexp_length_jia.pdf",
            outdir = '../figures/',
            width = 60,height = 30)
