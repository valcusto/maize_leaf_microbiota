######load packages
library(tidyverse)
library(ohchibi)
library(paletteer)
library(scales)
library(dplyr)
library(ggtree)
library(clusterProfiler) #GO analysis
library(biomaRt)

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

######load data
d1 <- readRDS('../cleandata/res_rnaseq_contrasts_syncom_do.RDS') # DESeq2 results
Dat_av  <- readRDS(file = "../cleandata/res_tab_av_rnaseq.RDS")
Tab_av <- Dat_av$Tab_av

######open the gene
term2gene <- read.csv("../rawdata/Zmays.GO.anno_rice_At.term_gene.txt", header=F, sep="\t")
term2name <- read.csv("../rawdata/Zmays.GO.anno_rice_At.term_name.txt", header=F, sep="\t")

######select only genes for L3
Res_sig  <- d1 %>% 
  subset(padj < 0.1 )  %>%
  dplyr::filter(.data = .,Contrast %>% grepl(pattern = "_L3_FullSynCom")) %>%
  droplevels

######extract only the Res_sig file
sig_genes <- Res_sig$gene_id %>% as.character %>%
  unique

######how many genes significant
length(sig_genes) # 476 genes to work with

######create the clusterRows dataset
Tab_sub <- Tab_av[,colnames(Tab_av) %>% grep(pattern = "^L3")]

Tab_sub <- match(sig_genes,rownames(Tab_sub)) %>%
  Tab_sub[.,]

res_heatmap <- chibi.heatmap(Tab = Tab_sub,
                             dist_method_rows = "euclidean",hclust_method_rows = "ward.D",
                             range_fill_heatmap = c(-2,2),
                             k_rows = 11,panel_spacing = 0.2,
                             k_cols = 5,size_axis_text_row = 6,axis_ticks_row = T)

a <- res_heatmap$heatmap
a

######extract the cluster information
df_clust_rows_l3 <- res_heatmap$df_clust_rows

######select each cluster for L3
cr1_l3 <- df_clust_rows_l3 %>%
  subset(ClusterRows == 'CR1') %>% droplevels
cr2_l3 <- df_clust_rows_l3 %>%
  subset(ClusterRows == 'CR2') %>% droplevels
cr3_l3 <- df_clust_rows_l3 %>%
  subset(ClusterRows == 'CR3') %>% droplevels
cr4_l3 <- df_clust_rows_l3 %>%
  subset(ClusterRows == 'CR4') %>% droplevels
cr5_l3 <- df_clust_rows_l3 %>%
  subset(ClusterRows == 'CR5') %>% droplevels
cr6_l3 <- df_clust_rows_l3 %>%
  subset(ClusterRows == 'CR6') %>% droplevels
cr7_l3 <- df_clust_rows_l3 %>%
  subset(ClusterRows == 'CR7') %>% droplevels
cr8_l3 <- df_clust_rows_l3 %>%
  subset(ClusterRows == 'CR8') %>% droplevels
cr9_l3 <- df_clust_rows_l3 %>%
  subset(ClusterRows == 'CR9') %>% droplevels
cr10_l3 <- df_clust_rows_l3 %>%
  subset(ClusterRows == 'CR10') %>% droplevels
cr11_l3 <- df_clust_rows_l3 %>%
  subset(ClusterRows == 'CR11') %>% droplevels

######cr1
gene_cr1_l3 <- as.factor(cr1_l3$IdRows)

######GO enrichment analysis using each clusters list
x_cr1_l3 <- enricher(gene_cr1_l3, TERM2GENE=term2gene, TERM2NAME=term2name,
              pvalueCutoff = 0.1, pAdjustMethod = "fdr", qvalueCutoff = 0.1)

######What is the percentage of genes with GO?
df_cr1_l3 <- split_and_repeat(x_cr1_l3@result, "geneID")
unique_cr1_l3 <- df_cr1_l3$geneID %>% unique
###answer:
cr1_l3_total <- length(gene_cr1_l3)
cr1_l3_go <- length(unique_cr1_l3)
prop_go_cr1_l3 <- cr1_l3_go/cr1_l3_total #63% GO identified

######cr2
gene_cr2_l3 <- as.factor(cr2_l3$IdRows)

######GO enrichment analysis using each clusters list
x_cr2_l3 <- enricher(gene_cr2_l3, TERM2GENE=term2gene, TERM2NAME=term2name,
                     pvalueCutoff = 0.1, pAdjustMethod = "fdr", 
                     qvalueCutoff = 0.1)

######view the results
View(x_cr2_l3@result)

######What is the percentage of genes with GO?
df_cr2_l3 <- split_and_repeat(x_cr2_l3@result, "geneID")
unique_cr2_l3 <- df_cr2_l3$geneID %>% unique
###answer:
cr2_l3_total <- length(gene_cr2_l3)
cr2_l3_go <- length(unique_cr2_l3)
prop_go_cr2_l3 <- cr2_l3_go/cr2_l3_total #65% GO identified

######cr3
gene_cr3_l3 <- as.factor(cr3_l3$IdRows)

######GO enrichment analysis using each clusters list
x_cr3_l3 <- enricher(gene_cr3_l3, TERM2GENE=term2gene, TERM2NAME=term2name,
                     pvalueCutoff = 0.1, pAdjustMethod = "fdr", 
                     qvalueCutoff = 0.1)

######view the results
View(x_cr3_l3@result)

######What is the percentage of genes with GO?
df_cr3_l3 <- split_and_repeat(x_cr3_l3@result, "geneID")
unique_cr3_l3 <- df_cr3_l3$geneID %>% unique
###answer:
cr3_l3_total <- length(gene_cr3_l3)
cr3_l3_go <- length(unique_cr3_l3)
prop_go_cr3_l3 <- cr3_l3_go/cr3_l3_total #54% GO identified

######cr4
gene_cr4_l3 <- as.factor(cr4_l3$IdRows)

######GO enrichment analysis using each clusters list
x_cr4_l3 <- enricher(gene_cr4_l3, TERM2GENE=term2gene, TERM2NAME=term2name,
                     pvalueCutoff = 0.1, pAdjustMethod = "fdr", 
                     qvalueCutoff = 0.1)

######view the results
View(x_cr4_l3@result)

######What is the percentage of genes with GO?
df_cr4_l3 <- split_and_repeat(x_cr4_l3@result, "geneID")
unique_cr4_l3 <- df_cr4_l3$geneID %>% unique
###answer:
cr4_l3_total <- length(gene_cr4_l3)
cr4_l3_go <- length(unique_cr4_l3)
prop_go_cr4_l3 <- cr4_l3_go/cr4_l3_total #76% GO identified

######cr5
gene_cr5_l3 <- as.factor(cr5_l3$IdRows)

######GO enrichment analysis using each clusters list
x_cr5_l3 <- enricher(gene_cr5_l3, TERM2GENE=term2gene, TERM2NAME=term2name,
                     pvalueCutoff = 0.1, pAdjustMethod = "fdr", 
                     qvalueCutoff = 0.1)

######view the results
View(x_cr5_l3@result)

######What is the percentage of genes with GO?
df_cr5_l3 <- split_and_repeat(x_cr5_l3@result, "geneID")
unique_cr5_l3 <- df_cr5_l3$geneID %>% unique
###answer:
cr5_l3_total <- length(gene_cr5_l3)
cr5_l3_go <- length(unique_cr5_l3)
prop_go_cr5_l3 <- cr5_l3_go/cr5_l3_total #58% GO identified

#######cr6
gene_cr6_l3 <- as.factor(cr6_l3$IdRows)

######GO enrichment analysis using each clusters list
x_cr6_l3 <- enricher(gene_cr6_l3, TERM2GENE=term2gene, TERM2NAME=term2name,
                     pvalueCutoff = 0.1, pAdjustMethod = "fdr", 
                     qvalueCutoff = 0.1)

######view the results
View(x_cr6_l3@result)

######What is the percentage of genes with GO?
df_cr6_l3 <- split_and_repeat(x_cr6_l3@result, "geneID")
unique_cr6_l3 <- df_cr6_l3$geneID %>% unique
###answer:
cr6_l3_total <- length(gene_cr6_l3)
cr6_l3_go <- length(unique_cr6_l3)
prop_go_cr6_l3 <- cr6_l3_go/cr6_l3_total #81% GO identified

#######cr7
gene_cr7_l3 <- as.factor(cr7_l3$IdRows)

######GO enrichment analysis using each clusters list
x_cr7_l3 <- enricher(gene_cr7_l3, TERM2GENE=term2gene, TERM2NAME=term2name,
                     pvalueCutoff = 0.1, pAdjustMethod = "fdr", 
                     qvalueCutoff = 0.1)

######view the results
View(x_cr7_l3@result)

######What is the percentage of genes with GO?
df_cr7_l3 <- split_and_repeat(x_cr7_l3@result, "geneID")
unique_cr7_l3 <- df_cr7_l3$geneID %>% unique
###answer:
cr7_l3_total <- length(gene_cr7_l3)
cr7_l3_go <- length(unique_cr7_l3)
prop_go_cr7_l3 <- cr7_l3_go/cr7_l3_total #77% GO identified

#######cr8
gene_cr8_l3 <- as.factor(cr8_l3$IdRows)

######GO enrichment analysis using each clusters list
x_cr8_l3 <- enricher(gene_cr8_l3, TERM2GENE=term2gene, TERM2NAME=term2name,
                     pvalueCutoff = 0.1, pAdjustMethod = "fdr", 
                     qvalueCutoff = 0.1)

######view the results
View(x_cr8_l3@result)

######What is the percentage of genes with GO?
df_cr8_l3 <- split_and_repeat(x_cr8_l3@result, "geneID")
unique_cr8_l3 <- df_cr8_l3$geneID %>% unique
###answer:
cr8_l3_total <- length(gene_cr8_l3)
cr8_l3_go <- length(unique_cr8_l3)
prop_go_cr8_l3 <- cr8_l3_go/cr8_l3_total #74% GO identified

#######cr9
gene_cr9_l3 <- as.factor(cr9_l3$IdRows)

######GO enrichment analysis using each clusters list
x_cr9_l3 <- enricher(gene_cr9_l3, TERM2GENE=term2gene, TERM2NAME=term2name,
                     pvalueCutoff = 0.1, pAdjustMethod = "fdr", 
                     qvalueCutoff = 0.1)

######view the results
View(x_cr9_l3@result)

######What is the percentage of genes with GO?
df_cr9_l3 <- split_and_repeat(x_cr9_l3@result, "geneID")
unique_cr9_l3 <- df_cr9_l3$geneID %>% unique
###answer:
cr9_l3_total <- length(gene_cr9_l3)
cr9_l3_go <- length(unique_cr9_l3)
prop_go_cr9_l3 <- cr9_l3_go/cr9_l3_total #70% GO identified

#######cr10
gene_cr10_l3 <- as.factor(cr10_l3$IdRows)

######GO enrichment analysis using each clusters list
x_cr10_l3 <- enricher(gene_cr10_l3, TERM2GENE=term2gene, TERM2NAME=term2name,
                     pvalueCutoff = 0.1, pAdjustMethod = "fdr", 
                     qvalueCutoff = 0.1)

######view the results
View(x_cr10_l3@result)

######What is the percentage of genes with GO?
df_cr10_l3 <- split_and_repeat(x_cr10_l3@result, "geneID")
unique_cr10_l3 <- df_cr10_l3$geneID %>% unique
###answer:
cr10_l3_total <- length(gene_cr10_l3)
cr10_l3_go <- length(unique_cr10_l3)
prop_go_cr10_l3 <- cr10_l3_go/cr10_l3_total #77% GO identified

#######cr11
gene_cr11_l3 <- as.factor(cr11_l3$IdRows)

######GO enrichment analysis using each clusters list
x_cr11_l3 <- enricher(gene_cr11_l3, TERM2GENE=term2gene, TERM2NAME=term2name,
                      pvalueCutoff = 0.1, pAdjustMethod = "fdr", 
                      qvalueCutoff = 0.1)

######view the results
View(x_cr11_l3@result)

######What is the percentage of genes with GO?
df_cr11_l3 <- split_and_repeat(x_cr11_l3@result, "geneID")
unique_cr11_l3 <- df_cr11_l3$geneID %>% unique
###answer:
cr11_l3_total <- length(gene_cr11_l3)
cr11_l3_go <- length(unique_cr11_l3)
prop_go_cr11_l3 <- cr11_l3_go/cr11_l3_total #72% GO identified

######merge results
df_results <- merge_result(list(CR9 = x_cr9_l3,
                                CR10=x_cr10_l3,
                                CR7 = x_cr7_l3,
                                CR4 = x_cr4_l3,
                                CR6 = x_cr6_l3,
                                CR8 = x_cr8_l3,
                                CR11=x_cr11_l3,
                                CR3 = x_cr3_l3,
                                CR2 = x_cr2_l3,
                                CR1 = x_cr1_l3,
                                CR5 = x_cr5_l3))

df_results_l3 <- merge_result(list(CR1 = x_cr1_l3,
                                   CR2 = x_cr2_l3,
                                   CR3 = x_cr3_l3,
                                   CR4 = x_cr4_l3,
                                   CR5 = x_cr5_l3,
                                   CR6 = x_cr6_l3,
                                   CR7 = x_cr7_l3,
                                   CR8 = x_cr8_l3,
                                   CR9 = x_cr9_l3,
                                   CR10=x_cr10_l3,
                                   CR11=x_cr11_l3))

######save the dataframe
df_results_l3_1 <- df_results_l3@compareClusterResult

write.table(df_results_l3_1, "../cleandata/df_results_GO_Jia_final_L3.txt", quote=F, 
            row.names=F, col.names=T, sep="\t")

#####dotplot
#How dotplot works?
#dotplot method implemented in clusterProfiler try to make comparison,
#among clusters more informative and reasonable. After extracting,
#e.g. 5 categories for each clusters (showcategory=5), clusterProfiler
#try to collect and overlap of these categories among clusters.
#https://guangchuangyu.github.io/2016/11/showcategory-parameter-for-visualizing-comparecluster-output/
p <- dotplot(df_results,showCategory=5, title = "L3", includeAll=FALSE)+
  scale_color_paletteer_c( "viridis::plasma",na.value = "#BFBFBF")

p <- p+
  ggtitle('L3')+
  #coord_flip() +
  clean+
  theme(plot.title = element_text(size=30, face = 'bold', hjust=0.5),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=12),
        legend.key.height = unit(1, "cm")) 

p

######save the L3 GO
oh.save.pdf(p = p,
            outname = "../figures/figS4A_Heatmap_GO_all_L3_maize_jia_signifcant.pdf",
            width = 20,height = 30)

######select only genes for L3
Res_sig_l4  <- d1 %>% 
  subset(padj < 0.1 )  %>%
  dplyr::filter(.data = .,Contrast %>% grepl(pattern = "_L4_FullSynCom")) %>%
  droplevels

######extract only the Res_sig file
sig_genes_l4 <- Res_sig_l4$gene_id %>% as.character %>%
  unique

######how many genes significant
length(sig_genes_l4) # 123 genes to work with

######create the clusterRows dataset
Tab_sub_l4 <- Tab_av[,colnames(Tab_av) %>% grep(pattern = "^L4")]

Tab_sub_l4 <- match(sig_genes_l4,rownames(Tab_sub_l4)) %>%
  Tab_sub_l4[.,]

res_heatmap_l4 <- chibi.heatmap(Tab = Tab_sub_l4,
                                dist_method_rows = "euclidean",hclust_method_rows = "ward.D",
                                range_fill_heatmap = c(-2,2),
                                k_rows = 8,panel_spacing = 0.2,
                                k_cols = 5,size_axis_text_row = 6,axis_ticks_row = T)

b <- res_heatmap_l4$heatmap
b

######extract the cluster information
df_clust_rows_l4 <- res_heatmap_l4$df_clust_rows

######select each cluster for L4
cr1_l4 <- df_clust_rows_l4 %>%
  subset(ClusterRows == 'CR1') %>% droplevels
cr2_l4 <- df_clust_rows_l4 %>%
  subset(ClusterRows == 'CR2') %>% droplevels
cr3_l4 <- df_clust_rows_l4 %>%
  subset(ClusterRows == 'CR3') %>% droplevels
cr4_l4 <- df_clust_rows_l4 %>%
  subset(ClusterRows == 'CR4') %>% droplevels
cr5_l4 <- df_clust_rows_l4 %>%
  subset(ClusterRows == 'CR5') %>% droplevels
cr6_l4 <- df_clust_rows_l4 %>%
  subset(ClusterRows == 'CR6') %>% droplevels
cr7_l4 <- df_clust_rows_l4 %>%
  subset(ClusterRows == 'CR7') %>% droplevels
cr8_l4 <- df_clust_rows_l4 %>%
  subset(ClusterRows == 'CR8') %>% droplevels

######cr1
gene_cr1_l4 <- as.factor(cr1_l4$IdRows)

######GO enrichment analysis using each clusters list
x_cr1_l4 <- enricher(gene_cr1_l4, TERM2GENE=term2gene, TERM2NAME=term2name,
                     pvalueCutoff = 0.1, pAdjustMethod = "fdr", 
                     qvalueCutoff = 0.1)

######view the results
View(x_cr1_l4@result)

######What is the percentage of genes with GO?
df_cr1_l4 <- split_and_repeat(x_cr1_l4@result, "geneID")
unique_cr1_l4 <- df_cr1_l4$geneID %>% unique
###answer:
cr1_l4_total <- length(gene_cr1_l4)
cr1_l4_go <- length(unique_cr1_l4)
prop_go_cr1_l4 <- cr1_l4_go/cr1_l4_total #70% GO identified

######cr2
gene_cr2_l4 <- as.factor(cr2_l4$IdRows)

######GO enrichment analysis using each clusters list
x_cr2_l4 <- enricher(gene_cr2_l4, TERM2GENE=term2gene, TERM2NAME=term2name,
                     pvalueCutoff = 0.1, pAdjustMethod = "fdr", 
                     qvalueCutoff = 0.1)

######view the results
View(x_cr2_l4@result)

######What is the percentage of genes with GO?
df_cr2_l4 <- split_and_repeat(x_cr2_l4@result, "geneID")
unique_cr2_l4 <- df_cr2_l4$geneID %>% unique
###answer:
cr2_l4_total <- length(gene_cr2_l4)
cr2_l4_go <- length(unique_cr2_l4)
prop_go_cr2_l4 <- cr2_l4_go/cr2_l4_total #66% GO identified

######cr3
gene_cr3_l4 <- as.factor(cr3_l4$IdRows)

######GO enrichment analysis using each clusters list
x_cr3_l4 <- enricher(gene_cr3_l4, TERM2GENE=term2gene, TERM2NAME=term2name,
                     pvalueCutoff = 0.1, pAdjustMethod = "fdr", 
                     qvalueCutoff = 0.1)

######view the results
View(x_cr3_l4@result)

######What is the percentage of genes with GO?
df_cr3_l4 <- split_and_repeat(x_cr3_l4@result, "geneID")
unique_cr3_l4 <- df_cr3_l4$geneID %>% unique
###answer:
cr3_l4_total <- length(gene_cr3_l4)
cr3_l4_go <- length(unique_cr3_l4)
prop_go_cr3_l4 <- cr3_l4_go/cr3_l4_total #40% GO identified

######cr4
gene_cr4_l4 <- as.factor(cr4_l4$IdRows)

######GO enrichment analysis using each clusters list
x_cr4_l4 <- enricher(gene_cr4_l4, TERM2GENE=term2gene, TERM2NAME=term2name,
                     pvalueCutoff = 0.1, pAdjustMethod = "fdr", qvalueCutoff = 0.1)

######view the results
View(x_cr4_l4@result)

######What is the percentage of genes with GO?
df_cr4_l4 <- split_and_repeat(x_cr4_l4@result, "geneID")
unique_cr4_l4 <- df_cr4_l4$geneID %>% unique
###answer:
cr4_l4_total <- length(gene_cr4_l4)
cr4_l4_go <- length(unique_cr4_l4)
prop_go_cr4_l4 <- cr4_l4_go/cr4_l4_total #36% GO identified

######cr5
gene_cr5_l4 <- as.factor(cr5_l4$IdRows)

######GO enrichment analysis using each clusters list
x_cr5_l4 <- enricher(gene_cr5_l4, TERM2GENE=term2gene, TERM2NAME=term2name,
                     pvalueCutoff = 0.1, pAdjustMethod = "fdr", 
                     qvalueCutoff = 0.1)

######view the results
View(x_cr5_l4@result)

######What is the percentage of genes with GO?
df_cr5_l4 <- split_and_repeat(x_cr5_l4@result, "geneID")
unique_cr5_l4 <- df_cr5_l4$geneID %>% unique
###answer:
cr5_l4_total <- length(gene_cr5_l4)
cr5_l4_go <- length(unique_cr5_l4)
prop_go_cr5_l4 <- cr5_l4_go/cr5_l4_total #57% GO identified

#######cr6
gene_cr6_l4 <- as.factor(cr6_l4$IdRows)

######GO enrichment analysis using each clusters list
x_cr6_l4 <- enricher(gene_cr6_l4, TERM2GENE=term2gene, TERM2NAME=term2name,
                     pvalueCutoff = 0.1, pAdjustMethod = "fdr", qvalueCutoff = 0.1)

######view the results
View(x_cr6_l4@result)

######What is the percentage of genes with GO?
df_cr6_l4 <- split_and_repeat(x_cr6_l4@result, "geneID")
unique_cr6_l4 <- df_cr6_l4$geneID %>% unique
###answer:
cr6_l4_total <- length(gene_cr6_l4)
cr6_l4_go <- length(unique_cr6_l4)
prop_go_cr6_l4 <- cr6_l4_go/cr6_l4_total #58% GO identified

#######cr7
gene_cr7_l4 <- as.factor(cr7_l4$IdRows)

######GO enrichment analysis using each clusters list
x_cr7_l4 <- enricher(gene_cr7_l4, TERM2GENE=term2gene, TERM2NAME=term2name,
                     pvalueCutoff = 0.1, pAdjustMethod = "fdr", qvalueCutoff = 0.1)

######view the results
View(x_cr7_l4@result)

######What is the percentage of genes with GO?
df_cr7_l4 <- split_and_repeat(x_cr7_l4@result, "geneID")
unique_cr7_l4 <- df_cr7_l4$geneID %>% unique
###answer:
cr7_l4_total <- length(gene_cr7_l4)
cr7_l4_go <- length(unique_cr7_l4)
prop_go_cr7_l4 <- cr7_l4_go/cr7_l4_total #60% GO identified

#######cr8
gene_cr8_l4 <- as.factor(cr8_l4$IdRows)

######GO enrichment analysis using each clusters list
x_cr8_l4 <- enricher(gene_cr8_l4, TERM2GENE=term2gene, TERM2NAME=term2name,
                     pvalueCutoff = 0.1, pAdjustMethod = "fdr",
                     qvalueCutoff = 0.1)

######view the results
View(x_cr8_l4@result)

######What is the percentage of genes with GO?
df_cr8_l4 <- split_and_repeat(x_cr8_l4@result, "geneID")
unique_cr8_l4 <- df_cr8_l4$geneID %>% unique
###answer:
cr8_l4_total <- length(gene_cr8_l4)
cr8_l4_go <- length(unique_cr8_l4)
prop_go_cr8_l4 <- cr8_l4_go/cr8_l4_total #42% GO identified

######merge results
df_results_l4_final <- merge_result(list(CR7 = x_cr7_l4,
                                         CR4 = x_cr4_l4,
                                         CR6 = x_cr6_l4,
                                         CR8 = x_cr8_l4,
                                         CR3 = x_cr3_l4,
                                         CR2 = x_cr2_l4,
                                         CR1 = x_cr1_l4,
                                         CR5 = x_cr5_l4))

df_results_l4 <- merge_result(list(CR1 = x_cr1_l4,
                                   CR2 = x_cr2_l4,
                                   CR3 = x_cr3_l4,
                                   CR4 = x_cr4_l4,
                                   CR5 = x_cr5_l4,
                                   CR6 = x_cr6_l4,
                                   CR7 = x_cr7_l4,
                                   CR8 = x_cr8_l4))

######save the dataframe
df_results_l4_1 <- df_results_l4@compareClusterResult

write.table(df_results_l4_1, "../cleandata/df_results_GO_Jia_final_L4.txt", quote=F, 
            row.names=F, col.names=T, sep="\t")

#####dotplot
#How dotplot works?
#dotplot method implemented in clusterProfiler try to make comparison,
#among clusters more informative and reasonable. After extracting,
#e.g. 5 categories for each clusters (showcategory=5), clusterProfiler
#try to collect and overlap of these categories among clusters.
#https://guangchuangyu.github.io/2016/11/showcategory-parameter-for-visualizing-comparecluster-output/
p_l4 <- dotplot(df_results_l4_final,showCategory=5, title = "L4", 
                includeAll=FALSE)+
  scale_color_paletteer_c( "viridis::plasma",na.value = "#BFBFBF")

p_l4 <- p_l4+
  ggtitle('L4')+
  clean+
  theme(plot.title = element_text(size=30, face = 'bold', hjust=0.5),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=12),
        legend.key.height = unit(1, "cm")) 
p_l4

######save the l4 GO
oh.save.pdf(p = p_l4,
            outname = "../figures/figS4A_Heatmap_GO_all_L4_maize_jia.pdf",
            width = 20,height = 30)

######select only genes for L5
Res_sig_l5  <- d1 %>% 
  subset(padj < 0.1 )  %>%
  dplyr::filter(.data = .,Contrast %>% grepl(pattern = "_L5_FullSynCom")) %>%
  droplevels

######extract only the Res_sig file
sig_genes_l5 <- Res_sig_l5$gene_id %>% as.character %>%
  unique

######how many genes significant
length(sig_genes_l5) # 57 genes to work with

######create the clusterRows dataset
Tab_sub_l5 <- Tab_av[,colnames(Tab_av) %>% grep(pattern = "^L5")]

Tab_sub_l5 <- match(sig_genes_l5,rownames(Tab_sub_l5)) %>%
  Tab_sub_l5[.,]

res_heatmap_l5 <- chibi.heatmap(Tab = Tab_sub_l5,
                                dist_method_rows = "euclidean",hclust_method_rows = "ward.D",
                                range_fill_heatmap = c(-2,2),
                                k_rows = 4,panel_spacing = 0.2,
                                k_cols = 5,size_axis_text_row = 6,axis_ticks_row = T)

c <- res_heatmap_l5$heatmap
c

######extract the cluster information
df_clust_rows_l5 <- res_heatmap_l5$df_clust_rows

######select each cluster for l5
cr1_l5 <- df_clust_rows_l5 %>%
  subset(ClusterRows == 'CR1') %>% droplevels
cr2_l5 <- df_clust_rows_l5 %>%
  subset(ClusterRows == 'CR2') %>% droplevels
cr3_l5 <- df_clust_rows_l5 %>%
  subset(ClusterRows == 'CR3') %>% droplevels
cr4_l5 <- df_clust_rows_l5 %>%
  subset(ClusterRows == 'CR4') %>% droplevels

######cr1
gene_cr1_l5 <- as.factor(cr1_l5$IdRows)

######GO enrichment analysis using each clusters list
x_cr1_l5 <- enricher(gene_cr1_l5, TERM2GENE=term2gene, TERM2NAME=term2name,
                     pvalueCutoff = 0.1, pAdjustMethod = "fdr", 
                     qvalueCutoff = 0.1)

######view the results
View(x_cr1_l5@result)

######What is the percentage of genes with GO?
df_cr1_l5 <- split_and_repeat(x_cr1_l5@result, "geneID")
unique_cr1_l5 <- df_cr1_l5$geneID %>% unique
###answer:
cr1_l5_total <- length(gene_cr1_l5)
cr1_l5_go <- length(unique_cr1_l5)
prop_go_cr1_l5 <- cr1_l5_go/cr1_l5_total #33% GO identified

######cr2
gene_cr2_l5 <- as.factor(cr2_l5$IdRows)

######GO enrichment analysis using each clusters list
x_cr2_l5 <- enricher(gene_cr2_l5, TERM2GENE=term2gene, TERM2NAME=term2name,
                     pvalueCutoff = 0.1, pAdjustMethod = "fdr", 
                     qvalueCutoff = 0.1)

######view the results
View(x_cr2_l5@result)

######What is the percentage of genes with GO?
df_cr2_l5 <- split_and_repeat(x_cr2_l5@result, "geneID")
unique_cr2_l5 <- df_cr2_l5$geneID %>% unique
###answer:
cr2_l5_total <- length(gene_cr2_l5)
cr2_l5_go <- length(unique_cr2_l5)
prop_go_cr2_l5 <- cr2_l5_go/cr2_l5_total #25% GO identified

######cr3
gene_cr3_l5 <- as.factor(cr3_l5$IdRows)

######GO enrichment analysis using each clusters list
x_cr3_l5 <- enricher(gene_cr3_l5, TERM2GENE=term2gene, TERM2NAME=term2name,
                     pvalueCutoff = 0.1, pAdjustMethod = "fdr", 
                     qvalueCutoff = 0.1)

######view the results
View(x_cr3_l5@result)

######What is the percentage of genes with GO?
df_cr3_l5 <- split_and_repeat(x_cr3_l5@result, "geneID")
unique_cr3_l5 <- df_cr3_l5$geneID %>% unique
###answer:
cr3_l5_total <- length(gene_cr3_l5)
cr3_l5_go <- length(unique_cr3_l5)
prop_go_cr3_l5 <- cr3_l5_go/cr3_l5_total #83% GO identified

######cr4
gene_cr4_l5 <- as.factor(cr4_l5$IdRows)

######GO enrichment analysis using each clusters list
x_cr4_l5 <- enricher(gene_cr4_l5, TERM2GENE=term2gene, TERM2NAME=term2name,
                     pvalueCutoff = 0.1, pAdjustMethod = "fdr", qvalueCutoff = 0.1)

######view the results
View(x_cr4_l5@result)

######What is the percentage of genes with GO?
df_cr4_l5 <- split_and_repeat(x_cr4_l5@result, "geneID")
unique_cr4_l5 <- df_cr4_l5$geneID %>% unique
###answer:
cr4_l5_total <- length(gene_cr4_l5)
cr4_l5_go <- length(unique_cr4_l5)
prop_go_cr4_l5 <- cr4_l5_go/cr4_l5_total #76% GO identified

######merge results
df_results_l5_final <- merge_result(list(CR4 = x_cr4_l5,
                                         CR3 = x_cr3_l5,
                                         CR2 = x_cr2_l5,
                                         CR1 = x_cr1_l5))

df_results_l5 <- merge_result(list(CR1 = x_cr1_l5,
                                   CR2 = x_cr2_l5,
                                   CR3 = x_cr3_l5,
                                   CR4 = x_cr4_l5))

######save the dataframe
df_results_l5_1 <- df_results_l5@compareClusterResult

write.table(df_results_l5_1, "../cleandata/df_results_GO_Jia_final_L5.txt", quote=F, 
            row.names=F, col.names=T, sep="\t")

#####dotplot
#How dotplot works?
#dotplot method implemented in clusterProfiler try to make comparison,
#among clusters more informative and reasonable. After extracting,
#e.g. 5 categories for each clusters (showcategory=5), clusterProfiler
#try to collect and overlap of these categories among clusters.
#https://guangchuangyu.github.io/2016/11/showcategory-parameter-for-visualizing-comparecluster-output/
p_l5 <- dotplot(df_results_l5_final,showCategory=5, title = "l5", 
                includeAll=FALSE)+
  scale_color_paletteer_c( "viridis::plasma",na.value = "#BFBFBF")

p_l5 <- p_l5+
  ggtitle('L5')+
  clean+
  theme(plot.title = element_text(size=30, face = 'bold', hjust=0.5),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=12),
        legend.key.height = unit(1, "cm")) 
p_l5

######save the l5 GO
oh.save.pdf(p = p_l5,
            outname = "../figures/figS4A_Heatmap_GO_all_L5_maize_jia.pdf",
            width = 20,height = 30)


######composition
composition = egg::ggarrange(p, p_l4, p_l5,
                             nrow=1)

######save the plot as pdf
oh.save.pdf(p = composition,
            outname = "figS4C_Heatmap_GO_all_clusters_maize_jia_final.pdf",
            width = 60,height = 30)
