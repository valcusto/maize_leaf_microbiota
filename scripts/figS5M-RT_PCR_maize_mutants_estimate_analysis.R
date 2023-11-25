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
#Define palette of colors for genes
paleta_gene <- c("#01665e","#5ab4ac",
                 "#3288bd","#e08214",
                 "#d8b365","#8c510a")
names(paleta_gene) <- c("Zm00001eb022810", "Zm00001eb381280",
                        "Zm00001eb290350", "Zm00001eb398010",
                        "Zm00001eb089630","Zm00001eb149250")

######load data
df <- read.csv('../rawdata/rt_pcr_final_for_fold_change.csv')

######correct the names for the treatment, genes and mutants
df$Sample_ID <- df$Sample_ID %>% as.character
df$Treatment[which(df$Treatment == 'Full_SynCom')] <- '+SynCom'
df$gene[which(df$gene == 'gene2')] <- 'Zm00001eb022810'
df$gene[which(df$gene == 'gene3')] <- 'Zm00001eb381280'
df$gene[which(df$gene == 'gene5')] <- 'Zm00001eb290350'
df$gene[which(df$gene == 'gene18')] <- 'Zm00001eb398010'
df$gene[which(df$gene == 'gene24')] <- 'Zm00001eb089630'
df$gene[which(df$gene == 'gene36')] <- 'Zm00001eb149250'

df$Sample_ID[which(df$Sample_ID == '2')] <- 'zm00001eb022810'
df$Sample_ID[which(df$Sample_ID == '5')] <- 'zm00001eb290350'
df$Sample_ID[which(df$Sample_ID == '18')] <- 'zm0000eb398010'

######select the W22
df_w22 <- df %>%
  subset(Sample_ID == 'W22') %>% droplevels

###create the unique names for the for loop
mleaf <- df_w22$NumLeaf %>% unique
mgene <- df_w22$gene %>% unique

###Empty dataset
Res_em <- NULL
Res_pval <- NULL

###perform the for loop for ttest
for (i in mgene) {
  df_gene <- df_w22 %>% subset(gene == i) %>% droplevels
  for (leaf in mleaf){
    df_leaf <- df_gene %>% subset(NumLeaf == leaf) %>% droplevels
      #linear model
    m1 <- lm(formula = expression ~ Treatment,data = df_leaf)
      ######
    m1_em <- emmeans(m1,pairwise ~ Treatment,
                     adjust = "none")
    df_em <- m1_em$emmeans %>% as.data.frame
    df_pval <- m1_em$contrasts %>% as.data.frame
    df_pval$gene<- i
    df_pval$NumLeaf<- leaf
    df_em$gene  <- i
    df_em$NumLeaf  <- leaf
    #Extract pvalue
    Res_em <- rbind(Res_em,df_em)
    Res_pval <- rbind(Res_pval,df_pval)
  }
}

#####add significance on the column
pthres = 0.1
Res_pval$Significance <- "NS"
Res_pval$Significance[which(Res_pval$p.value < pthres)] <- 'p<0.1'

######select the column 7 to 10
pval <- Res_pval[,c(7:9)]

#####merge Res_em and pval
Res <- merge(Res_em, pval, by=c('gene', 'NumLeaf'))

######add W22 column
Res$Sample_ID <- 'W22'

#####order the Treatment factor
Res$Treatment <- Res$Treatment %>%
  factor(levels = c('NB', '+SynCom'))

######organize the factor
Res$gene <- Res$gene %>%
  factor(levels = c("Zm00001eb022810", "Zm00001eb381280",
                    "Zm00001eb290350", "Zm00001eb398010",
                    "Zm00001eb089630","Zm00001eb149250"))

######Step10:Plot the informations
p1 <- ggplot(Res, aes(x=Treatment, y=emmean)) + 
  geom_line(aes(group = gene,color = Significance),size = 1) +
  geom_point(shape = 21,aes(fill = gene),size = 10) +
  facet_nested(.~Sample_ID+NumLeaf, scales = "free",space ="free")+
  theme_ohchibi(size_panel_border = 2) +
  scale_fill_manual(name = 'Genes',values = paleta_gene) +
  scale_color_manual(values = c("#D9D9D9","black")) +
  ylim(0.0,2.0)+
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
  theme(legend.position = "right") +
  theme(
    panel.background = element_rect(fill = "white"),
    strip.text.x = element_text(family = "Helvetica",face = "bold",size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=30),
    axis.text.x = element_text(size=20),
    axis.text.y = element_text(size=20)) 

p1

######save the W22
oh.save.pdf(p = p1,
            outname = "../figures/figS5M-W22_RT_PCR_expression.pdf",
            width = 15,height = 10)

######select the 022810
df_022810 <- df %>%
  subset(Sample_ID == 'zm00001eb022810') %>% droplevels

###create the unique names for the for loop
mleaf <- df_022810$NumLeaf %>% unique
mgene <- df_022810$gene %>% unique

###Empty dataset
Res_em <- NULL
Res_pval <- NULL

###perform the for loop for ttest
for (i in mgene) {
  df_gene <- df_022810 %>% subset(gene == i) %>% droplevels
  for (leaf in mleaf){
    df_leaf <- df_gene %>% subset(NumLeaf == leaf) %>% droplevels
    #linear model
    m1 <- lm(formula = expression ~ Treatment,data = df_leaf)
    ######
    m1_em <- emmeans(m1,pairwise ~ Treatment,
                     adjust = "none")
    df_em <- m1_em$emmeans %>% as.data.frame
    df_pval <- m1_em$contrasts %>% as.data.frame
    df_pval$gene<- i
    df_pval$NumLeaf<- leaf
    df_em$gene  <- i
    df_em$NumLeaf  <- leaf
    #Extract pvalue
    Res_em <- rbind(Res_em,df_em)
    Res_pval <- rbind(Res_pval,df_pval)
  }
}

#####add significance on the column
pthres = 0.1
Res_pval$Significance <- "NS"
Res_pval$Significance[which(Res_pval$p.value < pthres)] <- 'p<0.1'

######select the column 7 to 10
pval <- Res_pval[,c(7:9)]

#####merge Res_em and pval
Res <- merge(Res_em, pval, by=c('gene', 'NumLeaf'))

######add zm00001eb022810 column
Res$Sample_ID <- 'zm00001eb022810'

#####order the Treatment factor
Res$Treatment <- Res$Treatment %>%
  factor(levels = c('NB', '+SynCom'))

######organize the factor
Res$gene <- Res$gene %>%
  factor(levels = c("Zm00001eb022810", "Zm00001eb381280",
                    "Zm00001eb290350", "Zm00001eb398010",
                    "Zm00001eb089630","Zm00001eb149250"))

######Step10:Plot the informations
p2 <- ggplot(Res, aes(x=Treatment, y=emmean)) + 
  geom_line(aes(group = gene,color = Significance),size = 1) +
  geom_point(shape = 21,aes(fill = gene),size = 10) +
  facet_nested(.~Sample_ID+NumLeaf, scales = "free",space ="free")+
  theme_ohchibi(size_panel_border = 2) +
  scale_fill_manual(name = 'Genes',values = paleta_gene) +
  scale_color_manual(values = c("#D9D9D9","black")) +
  #ylim(0.0,2.0)+
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
  theme(legend.position = "right") +
  theme(
    panel.background = element_rect(fill = "white"),
    strip.text.x = element_text(family = "Helvetica",face = "bold",size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=30),
    axis.text.x = element_text(size=20),
    axis.text.y = element_text(size=20)) 

p2

######save the W22
oh.save.pdf(p = p2,
            outname = "../figures/figS5M-zm00001eb022810_RT_PCR_expression.pdf",
            width = 15,height = 10)
######select the 290350
df_290350 <- df %>%
  subset(Sample_ID == 'zm00001eb290350') %>% droplevels

###create the unique names for the for loop
mleaf <- df_290350$NumLeaf %>% unique
mgene <- df_290350$gene %>% unique

###Empty dataset
Res_em <- NULL
Res_pval <- NULL

###perform the for loop for ttest
for (i in mgene) {
  df_gene <- df_290350 %>% subset(gene == i) %>% droplevels
  for (leaf in mleaf){
    df_leaf <- df_gene %>% subset(NumLeaf == leaf) %>% droplevels
    #linear model
    m1 <- lm(formula = expression ~ Treatment,data = df_leaf)
    ######
    m1_em <- emmeans(m1,pairwise ~ Treatment,
                     adjust = "none")
    df_em <- m1_em$emmeans %>% as.data.frame
    df_pval <- m1_em$contrasts %>% as.data.frame
    df_pval$gene<- i
    df_pval$NumLeaf<- leaf
    df_em$gene  <- i
    df_em$NumLeaf  <- leaf
    #Extract pvalue
    Res_em <- rbind(Res_em,df_em)
    Res_pval <- rbind(Res_pval,df_pval)
  }
}

#####add significance on the column
pthres = 0.1
Res_pval$Significance <- "NS"
Res_pval$Significance[which(Res_pval$p.value < pthres)] <- 'p<0.1'

######select the column 7 to 10
pval <- Res_pval[,c(7:9)]

#####merge Res_em and pval
Res <- merge(Res_em, pval, by=c('gene', 'NumLeaf'))

######add zm00001eb022810 column
Res$Sample_ID <- 'zm00001eb290350'

#####order the Treatment factor
Res$Treatment <- Res$Treatment %>%
  factor(levels = c('NB', '+SynCom'))

######organize the factor
Res$gene <- Res$gene %>%
  factor(levels = c("Zm00001eb022810", "Zm00001eb381280",
                    "Zm00001eb290350", "Zm00001eb398010",
                    "Zm00001eb089630","Zm00001eb149250"))

######Step10:Plot the informations
p2 <- ggplot(Res, aes(x=Treatment, y=emmean)) + 
  geom_line(aes(group = gene,color = Significance),size = 1) +
  geom_point(shape = 21,aes(fill = gene),size = 10) +
  facet_nested(.~Sample_ID+NumLeaf, scales = "free",space ="free")+
  theme_ohchibi(size_panel_border = 2) +
  scale_fill_manual(name = 'Genes',values = paleta_gene) +
  scale_color_manual(values = c("#D9D9D9","black")) +
  #ylim(0.0,2.0)+
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
  theme(legend.position = "right") +
  theme(
    panel.background = element_rect(fill = "white"),
    strip.text.x = element_text(family = "Helvetica",face = "bold",size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=30),
    axis.text.x = element_text(size=20),
    axis.text.y = element_text(size=20)) 

p2

######save the W22
oh.save.pdf(p = p2,
            outname = "../figures/figS5M-zm00001eb290350_RT_PCR_expression.pdf",
            width = 15,height = 10)
