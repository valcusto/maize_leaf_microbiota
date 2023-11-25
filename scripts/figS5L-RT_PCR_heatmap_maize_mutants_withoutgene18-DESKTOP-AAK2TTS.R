######load package
library(ohchibi)
library(paletteer)
library(scales)
library(ggh4x)
library(emmeans)

######set seed
set.seed(12456)

######set work directory
setwd("C:/Users/Castrillo Lab/OneDrive - The University of Nottingham/Publication 2019/Paper preparation/manuscripts/scripts")

######helper functions
source("0-Clean_up_plots.R")
source('0_11-theme_ohchibi_2.R')

######load data
df <- read.csv('../rawdata/rt_pcr_final_for_fold_change.csv')

######correct the names for the treatment, genes and mutants
df$Sample_ID <- df$Sample_ID %>% as.character
df$Treatment[which(df$Treatment == 'Full_SynCom')] <- '+SynCom'
df$gene[which(df$gene == 'gene2')] <- 'Zm00001eb022810'
df$gene[which(df$gene == 'gene3')] <- 'Zm00001eb381280'
df$gene[which(df$gene == 'gene5')] <- 'Zm00001eb290350'
df$gene[which(df$gene == 'gene18')] <- 'Zm0000eb398010'
df$gene[which(df$gene == 'gene24')] <- 'Zm0000eb089630'
df$gene[which(df$gene == 'gene36')] <- 'Zm00001eb149250'

df$Sample_ID[which(df$Sample_ID == '2')] <- 'zm00001eb022810'
df$Sample_ID[which(df$Sample_ID == '5')] <- 'zm00001eb290350'
df$Sample_ID[which(df$Sample_ID == '18')] <- 'zm0000eb398010'

######remove the zm0000eb398010
df_1 <- df %>%
  subset(Sample_ID != 'zm0000eb398010')

######Average the Contrast
Tab <- acast(data = df, 
             formula = gene ~ Sample_ID+Treatment+NumLeaf, 
             value.var = "expression",fun.aggregate = mean)

#######Cluster and order the family patterns
mclust <- hclust(d = (1-cor(Tab %>%t)) %>% as.dist,method = "ward.D2")
order <- mclust$order %>% mclust$labels[.]

######Visualise the dendogram
mclust %>% plot

###### Melt the Tab structures with the structure used to plot
melted <- Tab %>% melt
colnames(melted) <- c("gene","group","mean")

#####Divide the dendogram in cluster
df_samples <- mclust %>% dendextend::cutree(tree=., k=2) %>%
  data.frame(gene=names(.), Clustgene=paste0("Cl", .),
             row.names = NULL) 

#remove first column
df_samples <- df_samples[,-1]

######merge the melted and df_samples
melted_sum <- merge(melted, df_samples, by='gene')

#####separate the mutants from the treatment
melted_sum$group <- melted_sum$group %>% as.character
melted_sum$Sample_ID <- sapply(strsplit(melted_sum$group, '_'), `[`, 1)
melted_sum$Treatment <- sapply(strsplit(melted_sum$group, '_'), `[`, 2)
melted_sum$NumLeaf <- sapply(strsplit(melted_sum$group, '_'), `[`, 3)

######organize the factor
melted_sum$Sample_ID <- melted_sum$Sample_ID %>%
  factor(levels = c("W22",
                    'zm00001eb022810',
                    'zm00001eb290350'))

melted_sum$Treatment <- melted_sum$Treatment %>%
  factor(levels = c('NB', '+SynCom'))

melted_sum$NumLeaf <- melted_sum$NumLeaf %>%
  factor(levels = c('L3', 'L4', 'L5'))

melted_sum$gene <- melted_sum$gene %>%
  factor(levels = c("Zm0000eb398010","Zm0000eb089630","Zm00001eb381280",
                    "Zm00001eb290350","Zm00001eb149250","Zm00001eb022810"))

###create the unique names for the for loop
mleaf <- df$NumLeaf %>% unique
mgene <- df$gene %>% unique
mmutants <- df$Sample_ID %>% unique

###Empty dataset
Res_em <- NULL
Res_pval <- NULL

###perform the for loop for ttest
for (i in mgene) {
  df_gene <- df %>% subset(gene == i) %>% droplevels
  for (m in mmutants) {
    df_mutants <- df_gene %>% subset(Sample_ID == m) %>% droplevels
    for (leaf in mleaf){
      df_leaf <- df_mutants %>% subset(NumLeaf == leaf) %>% droplevels
      #linear model
      m1 <- lm(formula = expression ~ Treatment,data = df_leaf)
      ######
      m1_em <- emmeans(m1,pairwise ~ Treatment,
                       adjust = "none")
      df_em <- m1_em$emmeans %>% as.data.frame
      df_pval <- m1_em$contrasts %>% as.data.frame
      df_pval$gene<- i
      df_pval$Sample_ID <- m
      df_pval$NumLeaf<- leaf
      df_em$gene  <- i
      df_em$Sample_ID <- m
      df_em$NumLeaf  <- leaf
      #Extract pvalue
      Res_em <- rbind(Res_em,df_em)
      Res_pval <- rbind(Res_pval,df_pval)
  }
   }
}

#####add significance on the column
pthres = 0.1
Res_pval$Significance <- "NS"
Res_pval$Significance[which(Res_pval$p.value < pthres)] <- 'p<0.1'

######select the column 7 to 10
pval <- Res_pval[,c(7:10)]

#####merge Res_em and pval
Res <- merge(Res_em, pval, by=c('gene', 'Sample_ID', 'NumLeaf'))

######combine the Res with df_2
melted_sum <- merge(melted_sum, Res, 
                    by=c('Treatment', 'gene', 'Sample_ID', 'NumLeaf'))

######heatmap
p_heatmap_cl1 <- ggplot(data = melted_sum, mapping = aes(x=Treatment,y=gene))+
  geom_line(aes(color = Significance),size = 3) +
  geom_raster(aes(fill=mean))+
  geom_tile(aes(color = Significance),fill = '#00000000', 
            size = 1.5,width = 0.99,height = 0, show.legend = F)+
  facet_nested(.~Sample_ID+NumLeaf, scales = "free",space ="free")+
  scale_fill_gradient2(name = 'Standardized\nexpression',low = "#021567",
                       mid = "white",high = "#d95f0e",
                       midpoint = 0,limits = c(0,3),oob = squish)+
  scale_color_manual(name='Significance\nbetween\nNB and\nSynCom',
                     values = c("#00000000","#4d4d4d"),
                     na.value = "#00000000")+
  #geom_text(aes(label = Significance),size = 8) + 
  clean +
  theme_ohchibi( size_axis_title.x = 22,
                 size_axis_title.y = 22,
                 legend_proportion_size = 2,
                 size_title_text = 30,
                 size_legend_text = 20,
                 size_panel_border = 1.5,
                 size_lines_panel = 0,
                 size_axis_text.y = 22,
                 size_axis_text.x = 22) +
  theme(panel.spacing = unit(0.1, 'lines'),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

p_heatmap_cl1

######save figures
oh.save.pdf(p = p_heatmap_cl1,
            outname = "../figures/figS5L_maize_mutants_RT-PCR_heatmap_withoutgene18.pdf",
            outdir = "../figures/",width = 20,height = 10)
