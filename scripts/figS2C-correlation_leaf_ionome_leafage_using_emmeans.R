######Load packages
library(ohchibi)
library(dplyr)
library(egg)
library(ggrepel)

######set seed
set.seed(12356)

######set work directory
setwd("C:/Users/Castrillo Lab/OneDrive - The University of Nottingham/Publication 2019/Paper preparation/manuscripts/scripts")

######helper functions
source('0-Clean_up_plots.R')

##########################emmeans ions Sao Domingos 50%FC
emmeans_SD <- read.csv("../cleandata/emmeans_ions_SD.csv", header = T, sep = ",",
                       comment.char = "",quote = "")

emmeans_SD$num_leaf <- emmeans_SD$num_leaf %>% 
  gsub(pattern = "L",replacement = "") %>% as.numeric

#Calculate the correlation coeffcient
cor_ions <- emmeans_SD$Ion %>% as.character %>% unique
Res_cor <- NULL
for(ion in cor_ions){
  cor_sub <- emmeans_SD %>%
    subset(subset = Ions == ion) %>% droplevels
  Res <- cor(cor_sub$emmean, cor_sub$num_leaf)
  df <- Res %>% as.data.frame
  df$Ions <- ion
  Res_cor <- rbind(Res_cor, df)
}

colnames(Res_cor) <- c("r","Ion")

#####Loop the p-value of each ion
Res_cor_p <- NULL
for(ion in cor_ions){
  cor_sub <- emmeans_SD %>%
    subset(subset = Ions == ion) %>% droplevels
  Res_p <- cor.test(cor_sub$emmean, cor_sub$num_leaf)$p.value
  df_p <- Res_p %>% as.data.frame
  df_p$Ions <- ion
  Res_cor_p <- rbind(Res_cor_p, df_p)
}

colnames(Res_cor_p) <- c("p.value","Ion")
Res_cor_p$p.adj <- Res_cor_p$p.value %>% p.adjust(method = "fdr")

##Merge
sig <- merge(Res_cor, Res_cor_p, by = "Ion")
sig$Significance <- rep("NoSignificant")
sig$Significance[which(sig$p.adj < 0.05)] <- "Significant"
sig$xp <- rep(5)
sig$yp <- rep (2)
sig$r <- format(round(sig$r, 2), nsmall = 2) %>% as.numeric

######order the coefficient values
sig <- sig[order(sig$r),]
order <- sig$Ion

######oder sig Ions
sig$Ion <- sig$Ion %>% factor(levels = order)

######plot the results
p_cor_sd <- ggplot(data=sig, aes(x=Ion,y=r, fill=Significance))+
  geom_bar(stat='identity')+
  xlab(NULL)+
  ylab('Correlation coefficient vs NumLeaf')+
  ggtitle('São Domingos')+
  scale_fill_manual(labels = c('No significant', 'Significance correlation (q < 0.05)'),
                    values=c('#999999', '#FF7F00')) + 
  clean +
  theme(legend.position = c(0.25, 0.9),
        legend.key.height= unit(0.3, 'cm'),
        legend.key.width = unit(0.6, 'cm'),
        legend.box.just = "right",
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=20),
        plot.title = element_text(hjust = 0.5, size=20, face = 'bold'))

p_cor_sd

######save the plot as pdf
oh.save.pdf(p = p_cor_sd,outname = "figS2B-correlation_ions_leaf.pdf",
            outdir = "../figures/",width = 15,height = 10)

##########################emmeans ions Tarrafal 50%FC
emmeans_t <- read.csv("../cleandata/emmeans_ions_tarrafal.csv", header = T, sep = ",",
                       comment.char = "",quote = "")

emmeans_t$num_leaf <- emmeans_t$num_leaf %>% 
  gsub(pattern = "L",replacement = "") %>% as.numeric

#Calculate the correlation coeffcient
cor_ions <- emmeans_t$Ions %>% as.character %>% unique
Res_cor <- NULL
for(ion in cor_ions){
  cor_sub <- emmeans_t %>%
    subset(subset = Ions == ion) %>% droplevels
  Res <- cor(cor_sub$emmean, cor_sub$num_leaf)
  df <- Res %>% as.data.frame
  df$Ion <- ion
  Res_cor <- rbind(Res_cor, df)
}

colnames(Res_cor) <- c("r","Ion")

#####Loop the p-value of each ion
Res_cor_p <- NULL
for(ion in cor_ions){
  cor_sub <- emmeans_t %>%
    subset(subset = Ions == ion) %>% droplevels
  Res_p <- cor.test(cor_sub$emmean, cor_sub$num_leaf)$p.value
  df_p <- Res_p %>% as.data.frame
  df_p$Ion <- ion
  Res_cor_p <- rbind(Res_cor_p, df_p)
}

colnames(Res_cor_p) <- c("p.value","Ion")
Res_cor_p$p.adj <- Res_cor_p$p.value %>% p.adjust(method = "fdr")

##Merge
sig <- merge(Res_cor, Res_cor_p, by = "Ion")
sig$Significance <- rep("NoSignificant")
sig$Significance[which(sig$p.adj < 0.05)] <- "Significant"
sig$xp <- rep(5)
sig$yp <- rep (2)
sig$r <- format(round(sig$r, 2), nsmall = 2) %>% as.numeric

######order the coefficient values
sig <- sig[order(sig$r),]
order <- sig$Ion

######oder sig Ions
sig$Ion <- sig$Ion %>% factor(levels = order)

######plot the results
p_cor_t <- ggplot(data=sig, aes(x=Ion,y=r, fill=Significance))+
  geom_bar(stat='identity')+
  xlab(NULL)+
  ylab('Correlation coefficient vs NumLeaf')+
  ggtitle('Tarrafal')+
  scale_fill_manual(labels = c('No significant', 'Significance correlation (q < 0.05)'),
                    values=c('#999999', '#A6CEE3')) + 
  clean +
  theme(legend.position = c(0.25, 0.9),
        legend.key.height= unit(0.3, 'cm'),
        legend.key.width = unit(0.6, 'cm'),
        legend.box.just = "right",
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=20),
        plot.title = element_text(hjust = 0.5, size=20, face = 'bold'))

p_cor_t

######save the plot as pdf
oh.save.pdf(p = p_cor_t,outname = "figS2B-correlation_ions_leaf_tarrafal.pdf",
            outdir = "../figures/",width = 15,height = 10)

##########################emmeans ions Sutton Bonington 50%FC
emmeans_sb <- read.csv("../cleandata/emmeans_ions_suttonbonington.csv", header = T, sep = ",",
                      comment.char = "",quote = "")

emmeans_sb$num_leaf <- emmeans_sb$num_leaf %>% 
  gsub(pattern = "L",replacement = "") %>% as.numeric

#Calculate the correlation coeffcient
cor_ions <- emmeans_sb$Ions %>% as.character %>% unique
Res_cor <- NULL
for(ion in cor_ions){
  cor_sub <- emmeans_sb %>%
    subset(subset = Ions == ion) %>% droplevels
  Res <- cor(cor_sub$emmean, cor_sub$num_leaf)
  df <- Res %>% as.data.frame
  df$Ion <- ion
  Res_cor <- rbind(Res_cor, df)
}

colnames(Res_cor) <- c("r","Ion")

#####Loop the p-value of each ion
Res_cor_p <- NULL
for(ion in cor_ions){
  cor_sub <- emmeans_sb %>%
    subset(subset = Ions == ion) %>% droplevels
  Res_p <- cor.test(cor_sub$emmean, cor_sub$num_leaf)$p.value
  df_p <- Res_p %>% as.data.frame
  df_p$Ion <- ion
  Res_cor_p <- rbind(Res_cor_p, df_p)
}

colnames(Res_cor_p) <- c("p.value","Ion")
Res_cor_p$p.adj <- Res_cor_p$p.value %>% p.adjust(method = "fdr")

##Merge
sig <- merge(Res_cor, Res_cor_p, by = "Ion")
sig$Significance <- rep("NoSignificant")
sig$Significance[which(sig$p.adj < 0.05)] <- "Significant"
sig$xp <- rep(5)
sig$yp <- rep (2)
sig$r <- format(round(sig$r, 2), nsmall = 2) %>% as.numeric

######order the coefficient values
sig <- sig[order(sig$r),]
order <- sig$Ion

######oder sig Ions
sig$Ion <- sig$Ion %>% factor(levels = order)

######plot the results
p_cor_sb <- ggplot(data=sig, aes(x=Ion,y=r, fill=Significance))+
  geom_bar(stat='identity')+
  xlab(NULL)+
  ylab('Correlation coefficient vs NumLeaf')+
  ggtitle('Sutton Bonington')+
  scale_fill_manual(labels = c('No significant','Significance correlation (q < 0.05)'),
                    values=c('#999999','#1F78B4')) + 
  clean +
  theme(legend.position = c(0.25, 0.9),
        legend.key.height= unit(0.3, 'cm'),
        legend.key.width = unit(0.6, 'cm'),
        legend.box.just = "right",
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=20),
        plot.title = element_text(hjust = 0.5, size=20, face = 'bold'))

p_cor_sb

######save the plot as pdf
oh.save.pdf(p = p_cor_sb,outname = "figS2B-correlation_ions_leaf_suttonbonington.pdf",
            outdir = "../figures/",width = 15,height = 10)
