######load packages
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

######upload data
df_length <-  read.csv("../cleandata/leaf_metadata_length.csv", 
                header = T, sep = ",", stringsAsFactors = T,
                quote = "",comment.char = "")

df_biomass <-  read.csv("../cleandata/leaf_biomass_corrected_uid.csv", header = T, sep = ",",
                quote = "",comment.char = "")

######correct the name in the first column
colnames(df_length)[1] <- 'UiD'

######select the 50%FC
df_sub_length <- df_length %>%
  subset((treatment == '50%FC') & (num_leaf != 'Stem')) %>% droplevels

df_sub_length <- df_sub_length %>%
  subset(leaf_lenght_cm != 'NA') %>% droplevels

df_sub_biomass <- df_biomass %>%
  subset((treatment == '50%FC') & (num_leaf != 'Stem')) %>% droplevels

df_sub_biomass <- df_sub_biomass %>%
  subset(dw_mg != 'NA') %>% droplevels

######change the na,e of sao domingos
df_sub_length$soil <- df_sub_length$soil %>%
  gsub(pattern = 'Sao Domingos', replacement = 'Sao Domingos')

df_sub_biomass$soil <- df_sub_biomass$soil %>%
  gsub(pattern = 'Sao Domingos', replacement = 'Sao Domingos')

######order the soil names
df_sub_length$soil <- df_sub_length$soil %>%
  factor(levels = c('Sao Domingos', 'Tarrafal','Sutton Bonington'))

df_sub_biomass$soil <- df_sub_biomass$soil %>%
  factor(levels = c('Sao Domingos', 'Tarrafal','Sutton Bonington'))

######combine the genotype, num_leaf and leaf
df_sub_biomass$group <- paste(df_sub_biomass$genotype,
                              df_sub_biomass$soil,
                              df_sub_biomass$num_leaf, sep='_')

df_sub_length$group <- paste(df_sub_length$genotype,
                             df_sub_length$soil,
                             df_sub_length$num_leaf, sep = '_')

######calculate the average and sum for each dataset
df_final_biomass <- aggregate(dw_mg ~ group,df_sub_biomass, sum)

df_final_lenght <- aggregate(leaf_lenght_cm ~ group,df_sub_length, mean)

#######look for the intersect in num_pot
usables_samples <- intersect(df_final_biomass$group, df_final_lenght$group)

######select the samples in the intersect
df_final_biomass <- match(usables_samples, df_final_biomass$group) %>%
  df_final_biomass[.,] %>% droplevels
df_final_lenght <- match(usables_samples, df_final_lenght$group) %>%
  df_final_lenght[.,] %>% droplevels

######combine the two datasets
df_final <- cbind(df_final_biomass,df_final_lenght)

######check if the correlation is linear
plot(dw_mg~leaf_lenght_cm, data = df_final)

######calculate the correlation
cortest <- cor.test(df_final$dw_mg, df_final$leaf_lenght_cm, 
                       use = "pairwise.complete.obs", method='pearson')
cor <- cortest$estimate %>% data.frame
colnames(cor) <- 'cor'
cor$pvalue <- cortest$p.value

######adjust the pvalue
cor$padjust <- cor$pvalue %>%
  p.adjust(p=., method='fdr')

######prepare the dataset for plotting
###add the significance plot
pthres <- 0.05
cor$significance <- rep('NS', nrow(cor))
cor$significance[which(cor$padj < pthres)] <- "q-value<0.05"
cor$significance <- cor$significance %>% factor

######select four unique column
df_final <- df_final[,c(1,2,4)]

######correlation dataset
mtext_t <- 'r=0.83 
p-value=3.23e-18'

######correlation plot
p_cor <- ggplot(data = df_final,aes(dw_mg,leaf_lenght_cm)) +
  geom_smooth(method = "lm",size = 2,color = "red",se = F) +
  geom_point() + 
  theme_ohchibi(size_panel_border = 2) +
  xlab(label = "Leaf dry weight (mg)") +
  ylab(label = "Leaf length (cm)") +
  clean +
  annotate(geom = "text",x = 900,y = 15,label = mtext_t,size = 8)
  

p_cor

######save figures
oh.save.pdf(p = p_cor,outname = "fig1D_correlation_leaf_length_leaf_dry_wgt.pdf",
            outdir = "../figures/",width = 20,height = 15)
