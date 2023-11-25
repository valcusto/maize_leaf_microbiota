######load packages
library(ohchibi)
library(paletteer)
library(scales)

######set seed
set.seed(130816)

######set work directory
setwd("C:/Users/Castrillo Lab/OneDrive - The University of Nottingham/Publication 2019/Paper preparation/manuscripts/scripts")

######helper functions
paleta_marker <- c('#67a9cf','#3690c0','#02818a','#016450', 'grey')
names(paleta_marker) <- c('L1&L2','L3','L4','L5', 'Unselected')
source('0-Clean_up_plots.R')
source('0_4-subset_Dataset.R')

######Read Dat object to compute the RA dataframe
Dat_all <- readRDS(file = "../cleandata/dat_syncomleaves_complete.end.RDS")
Dat_rab <- Dat_all$RelativeAbundance %>%
  subset.Dataset(Fraction == "Shoot" &  Bacteria == "+SynCom",drop = T,clean = T)
melted_rab <- Dat_rab$Tab %>% melt
colnames(melted_rab) <- c("Id","DADA2_Header","RA")
melted_rab <- merge(melted_rab,Dat_rab$Map,by = "DADA2_Header")

######Read the results of the model 
Res_all <- readRDS(file = "../cleandata/res_glm_syncomleaves_complete_5models.end.RDS")

######prevalence
df_prev <- Res_all$df_prev_asv

######add the total number of samples
df_prev$total <- 64

######prevalence
df_prev$prevalence <- df_prev$Freq/df_prev$total

######what is the prevalence for 10 samples
thresh <- 10/64

#Explore the model 
Res <- Res_all$GLMLeaves$ResNutrientSeparated
df_anova <- Res$Res_anova 
df_anova <- merge(df_anova,df_prev, by = "Id",all = TRUE)
df_anova$p.valuelog10 <- -log10(df_anova$p.value)
df_anova$SignificanceP <- NA
df_anova$SignificanceP[which(df_anova$p.value < 0.05)] <- "p < 0.05"
df_em <- Res$Res_em

#######  1/4 Hoagland #############
condition <- "1/4 Hoagland"

df_em_sub <- df_em %>%
  subset(Nutrient == condition) %>% droplevels
df_anova_sub <- (df_anova %>%subset(Nutrient == condition) %>% droplevels)[,-3]
melted_rab_sub <- melted_rab %>% 
  subset(Nutrient == condition) %>% droplevels
df_ra_sub <- aggregate(RA~Id,melted_rab_sub,mean)

######transform Id in character
df_em_sub$Id <- df_em_sub$Id %>% as.character

######create a column with prevalence > 15%
df_em_sub$leaf_marker <- 'Unselected'
df_em_sub$leaf_marker[which(df_em_sub$Id == 'ASV147')] <- 'L1&L2'
df_em_sub$leaf_marker[which(df_em_sub$Id == 'ASV49')] <- 'L1&L2'
df_em_sub$leaf_marker[which(df_em_sub$Id == 'ASV224')] <- 'L1&L2'

df_em_sub$leaf_marker[which(df_em_sub$Id == 'ASV71')] <- 'L3'
df_em_sub$leaf_marker[which(df_em_sub$Id == 'ASV155')] <- 'L3'
df_em_sub$leaf_marker[which(df_em_sub$Id == 'ASV83')] <- 'L3'

df_em_sub$leaf_marker[which(df_em_sub$Id == 'ASV52')] <- 'L4'
df_em_sub$leaf_marker[which(df_em_sub$Id == 'ASV45')] <- 'L4'
df_em_sub$leaf_marker[which(df_em_sub$Id == 'ASV9')] <- 'L4'
df_em_sub$leaf_marker[which(df_em_sub$Id == 'ASV34')] <- 'L4'
df_em_sub$leaf_marker[which(df_em_sub$Id == 'ASV2')] <- 'L4'
df_em_sub$leaf_marker[which(df_em_sub$Id == 'ASV7')] <- 'L4'
df_em_sub$leaf_marker[which(df_em_sub$Id == 'ASV62')] <- 'L4'

df_em_sub$leaf_marker[which(df_em_sub$Id == 'ASV10')] <- 'L5'
df_em_sub$leaf_marker[which(df_em_sub$Id == 'ASV43')] <- 'L5'
df_em_sub$leaf_marker[which(df_em_sub$Id == 'ASV31')] <- 'L5'

#######order the ASVs
Tab <- acast(data = df_em_sub,formula = NumLeaf~Id,value.var = "emmean") %>% 
  scale
mclust <- hclust(d = dist(Tab %>% t),method = "ward.D")
order_asv <- mclust$order %>% mclust$labels[.]

plot(mclust)

df_em_sub$Id <- df_em_sub$Id %>% factor(levels = order_asv)

melted_tab <- melt(Tab)
colnames(melted_tab) <- c("NumLeaf","Id","zscore")
df_em_sub <- merge(df_em_sub,melted_tab, by = c("NumLeaf","Id"),all = TRUE)

df_em_sub <- merge(df_em_sub,df_anova_sub,
                   by = "Id",all.x = TRUE)

df_em_sub <- merge(df_em_sub, df_ra_sub, by = "Id",all.x = TRUE)

######select the values to put in red
selected_values <- df_em_sub %>%
  subset(leaf_marker != 'Unselected') %>% droplevels

######select each leaf and check the zscore
select_l1 <- selected_values %>%
  subset(NumLeaf == 'L1' & leaf_marker == 'L1&L2') %>% droplevels
select_l2 <- selected_values %>%
  subset(NumLeaf == 'L2' & leaf_marker == 'L1&L2') %>% droplevels
select_l3 <- selected_values %>%
  subset(NumLeaf == 'L3' & leaf_marker == 'L3') %>% droplevels
select_l4 <- selected_values %>%
  subset(NumLeaf == 'L4' & leaf_marker == 'L4') %>% droplevels
select_l5 <- selected_values %>%
  subset(NumLeaf == 'L5' & leaf_marker == 'L5') %>% droplevels

selected_values <- selected_values$Id %>% unique %>% as.character

######plot the heatmap
p_low <- ggplot(data = df_em_sub,aes(NumLeaf,Id)) +
  geom_raster(aes(fill = zscore)) +
  facet_grid(.~Nutrient,space = "free",scales = "free") +
  scale_fill_paletteer_c(name = 'Standardized\nRelative Abundance' ,"pals::kovesi.diverging_bwr_55_98_c37") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  clean +
  theme_ohchibi() +
  theme(
    axis.text.y = element_text(size = 7,
                               color = c('black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'red', 'red', 'black', 'red',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black','red',
                                         'red','black', 'black', 'black', 
                                         'black',
                                         'black', 'black', 'black', 'red',
                                         'red',
                                         'black', 'black', 'black',
                                         'black','red', 'red', 'red', 'red',
                                         'black','black',
                                         'black','black',
                                         'black','black','black', 'red',
                                         'black','black','black','black',
                                         'black','black','black', 'red',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black',
                                         'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'red',
                                         'red', 'black', 'red', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black') %>% rev),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_blank(),
    legend.position = "bottom") 

p_low

p1 <- ggplot(data = df_em_sub[,c("Id","prevalence","leaf_marker")] %>% 
               unique,aes(Id,prevalence)) +
  geom_bar(stat = "identity",aes(fill = leaf_marker)) +
  geom_hline(yintercept = 0.15,linetype = "dashed",size = 0.5 , color = "#ef6548") +
  coord_flip() +
  scale_fill_paletteer_c("pals::kovesi.rainbow_bgyrm_35_85_c71",
                         limits = c(0,10), oob = squish) +
  scale_fill_manual(name = 'Leaf marker',values = paleta_marker) +
  scale_y_continuous(breaks = seq(0,1, 0.05),expand = c(0,0), 
                     labels = scales::percent) +
  #theme_ohchibi() + 
  scale_x_discrete(expand = c(0,0)) +
  clean +
  theme_ohchibi()+
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 7, 
                               color = c('black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'red', 'red', 'black', 'red',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black','red',
                                         'red','black', 'black', 'black', 
                                         'black',
                                         'black', 'black', 'black', 'red',
                                         'red',
                                         'black', 'black', 'black',
                                         'black','red', 'red', 'red', 'red',
                                         'black','black',
                                         'black','black',
                                         'black','black','black', 'red',
                                         'black','black','black','black',
                                         'black','black','black', 'red',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black',
                                         'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'red',
                                         'red', 'black', 'red', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black') %>% rev),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 9)
  ) +
  ylab(label = "Prevalence ASV")

p1

p1_1 <- ggplot(data = df_em_sub[,c("Id","RA", "leaf_marker")] %>% unique,aes(Id,RA)) +
  geom_bar(stat = "identity",aes(fill = leaf_marker)) +
  geom_hline(yintercept = 0.001,linetype = "dashed",size = 0.5 , color = "#ef6548") +
  coord_flip()  +
  scale_y_continuous(breaks = seq(0,0.20,0.01),expand = c(0,0)) +
  theme_ohchibi() + 
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(name = 'Leaf marker',values = paleta_marker) +
  #theme_ohchibi() + 
  scale_x_discrete(expand = c(0,0)) +
  clean +
  theme_ohchibi()+
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 7, 
                               color = c('black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'red', 'red', 'black', 'red',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black','red',
                                         'red','black', 'black', 'black', 
                                         'black',
                                         'black', 'black', 'black', 'red',
                                         'red',
                                         'black', 'black', 'black',
                                         'black','red', 'red', 'red', 'red',
                                         'black','black',
                                         'black','black',
                                         'black','black','black', 'red',
                                         'black','black','black','black',
                                         'black','black','black', 'red',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black',
                                         'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'red',
                                         'red', 'black', 'red', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black',
                                         'black', 'black', 'black', 'black') %>% rev),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 9))+
  ylab(label = "Cumulative RA in Shoot")

p1_1

######composition
composition <- egg::ggarrange(p_low, p1,p1_1, nrow = 1)

######save the plot
oh.save.pdf(p = composition,
            outname = "figS3L_leaf_marker_low_nutrient.pdf",
            outdir = "../figures/",width = 16,height = 12)

#####save the png file
png(filename = "../figure/figS3L_leaf_marker_low_nutrient.png", units='cm',
    width=600, height=350, res = 200)
composition
dev.off()

ggsave(plot = composition,'../figure/composition.png',
       units='cm',
       width=40, height=30)

### Plot correlation as heatmap
df <- Res_all$LMShootLength
df <- with(df,order(-r)) %>%
  df[.,]
order_id <- df$Id
df$Id <- df$Id %>% factor(levels = order_id)
df$p.valuelog10 <- -log10(df$p.value )
df$SignificanceP <- NA
df$SignificanceP[which(df$p.value < 0.05)] <- "p < 0.05"

######select the prevalence data
df_prev_leaf <- Res_all$df_prev_asv

######select the ASVs in df
prev_asv_df <- df$Id %>% unique %>% as.character

######select only the ASVs in df
df_prev_leaf1 <- df_prev_leaf[which(df_prev_leaf$Id %in% prev_asv_df),]

######add the total number of samples
df_prev_leaf1$total <- 64

######calculate the prevalence in %
df_prev_leaf1$prevalence <- df_prev_leaf1$Freq/df_prev_leaf1$total

######merge the two datasets
df <- merge(df, df_prev_leaf1, by = 'Id')

######add the growth column
df$growth <- 'Unselected'
df$growth[which(df$Id == 'ASV30')] <- 'growth_down'
df$growth[which(df$Id == 'ASV83')] <- 'growth_down'
df$growth[which(df$Id == 'ASV13')] <- 'growth_down'
df$growth[which(df$Id == 'ASV49')] <- 'growth_down'
df$growth[which(df$Id == 'ASV224')] <- 'growth_down'
df$growth[which(df$Id == 'ASV70')] <- 'growth_down'

df$growth[which(df$Id == 'ASV45')] <- 'growth_up'
df$growth[which(df$Id == 'ASV142')] <- 'growth_up'
df$growth[which(df$Id == 'ASV22')] <- 'growth_up'
df$growth[which(df$Id == 'ASV52')] <- 'growth_up'
df$growth[which(df$Id == 'ASV2')] <- 'growth_up'
df$growth[which(df$Id == 'ASV34')] <- 'growth_up'

######define the color for the growth
paleta_pheno <- c('#d8b365', '#f6e8c3', 'grey')
names(paleta_pheno) <- c('growth_down', 'growth_up', 'Unselected')

p8 <- ggplot(data = df,aes(Id,r)) +
  geom_bar(stat = "identity",aes(fill =growth,color =growth )) +
  geom_hline(yintercept = -0.3,linetype = "dashed",size = 0.5 , color = "#ef6548") +
  geom_hline(yintercept = 0.3,linetype = "dashed",size = 0.5 , color = "#ef6548") +
  coord_flip() +
  scale_fill_paletteer_c("pals::kovesi.rainbow_bgyrm_35_85_c71",
                         limits = c(0,10), oob = squish) +
  scale_fill_manual(name = 'Leaf marker',values = paleta_pheno) +
  scale_color_manual(name = 'Leaf marker',values = paleta_pheno) +
  clean +
  theme_ohchibi() + 
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(breaks = seq(-0.8,0.8,0.1),expand = c(0,0)) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 7,
                               color = c('black','black','black','black','black',
                                         'red','black','red','black','red',
                                         'black','red','black','black','black',
                                         'black','black','black','black','black',
                                         'black','red','red','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','red','black','black','black',
                                         'black','black','red','red','red',
                                         'red','red')%>% rev),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size=15),
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 9)
  ) +
  ylab(label = "Pearson r against leaf length")


p8

p9 <- ggplot(data = df[,c("Id","prevalence","growth")] %>% 
               unique,aes(Id,prevalence)) +
  geom_bar(stat = "identity",aes(fill = growth)) +
  geom_hline(yintercept = 0.15,linetype = "dashed",size = 0.5 , color = "#ef6548") +
  coord_flip() +
  scale_fill_paletteer_c("pals::kovesi.rainbow_bgyrm_35_85_c71",
                         limits = c(0,10), oob = squish) +
  scale_fill_manual(name = 'Leaf marker',values = paleta_pheno) +
  scale_y_continuous(breaks = seq(0,1, 0.05),expand = c(0,0), 
                     labels = scales::percent) +
  #theme_ohchibi() + 
  scale_x_discrete(expand = c(0,0)) +
  clean +
  theme_ohchibi()+
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 7, 
                               color = c('black','black','black','black','black',
                                         'red','black','red','black','red',
                                         'black','red','black','black','black',
                                         'black','black','black','black','black',
                                         'black','red','red','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','red','black','black','black',
                                         'black','black','red','red','red',
                                         'red','red') %>% rev),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 9)
  ) +
  ylab(label = "Prevalence ASV")

p9

######composition
composition <- egg::ggarrange(p9, p8, nrow = 1)

######save the dataset
oh.save.pdf(p = composition,
            outname = "figS3I_leaf_marker_phenotype_leaf_length.pdf",
            outdir = "../figures/",width = 10,height = 12)

##### Nutrient model ##########
Res <- Res_all$GLMNutrient$ResLeavesSeparated
df_em_sub<- Res$Res_em
df_pairs <- Res$Res_pairs

Tab <- acast(data = df_em_sub,formula = NumLeaf+Nutrient~Id,value.var = "emmean",
             fill = 0) %>% 
  scale
mclust <- hclust(d = dist(Tab %>% t),method = "ward.D")
order_asv <- mclust$order %>% mclust$labels[.]

df_em_sub$Id <- df_em_sub$Id %>% factor(levels = order_asv)

melted_tab <- melt(Tab)
colnames(melted_tab) <- c("UId","Id","zscore")

df_em_sub$UId <- paste0(df_em_sub$NumLeaf,"_",df_em_sub$Nutrient)
df_em_sub <- merge(df_em_sub,melted_tab, by = c("UId","Id"),all.x  = TRUE)

#Apend pvalue
df_em_sub <- merge(df_em_sub,df_pairs[,c("Id","NumLeaf","p.value")], by = c("Id","NumLeaf"),all.x = TRUE)
df_em_sub$SignificanceP <- 'NS'
df_em_sub$SignificanceP[which(df_em_sub$p.value < 0.05)] <- "p < 0.05"
df_em_sub$SignificanceP <- df_em_sub$SignificanceP %>%
  factor(levels = c('p < 0.05', 'NS'))

p7 <- ggplot(data = df_em_sub,aes(Nutrient,Id)) +
  geom_raster(aes(fill = zscore)) +
  geom_tile(aes(color = SignificanceP),fill = "#00000000",width = 0.9,height = 0.9,size = 0.5) +
  facet_grid(.~NumLeaf,space = "free",scales = "free") +
  scale_fill_paletteer_c(name = 'Standardized\nRelative abundance',"pals::kovesi.diverging_bwr_55_98_c37") +
  scale_color_manual(name = 'Significance',values = c("black", "#00000000")) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_ohchibi() +
  theme(axis.title.y = element_blank(),
    axis.text.y = element_text(size = 7,
                               color = c('red','black','red','red','black',
                                         'black','red','black','red','black',
                                         'red','black','red','red','red',
                                         'black','black','black','black','black',
                                         'red','red','red','red','red',
                                         'black','red','red','black','black',
                                         'red','black','black','black','black',
                                         'black','black','black','black','black',
                                         'red','red','red','red','black',
                                         'red','black','black','black','black',
                                         'black','black','red','black','red',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'red','black','black','black','black',
                                         'black','black','red','red','black',
                                         'black','red','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black', 'black', 'black', 'black', 'black',
                                         'black', 'black') %>% rev),
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 8),
    legend.position = "bottom"
  )

p7

### Plot prevalence
df_prev$Id <- df_prev$Id %>% factor(levels = order_asv)

######add the nutrient_marker column
df_prev$nutrient_marker <- 'Unselected'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV25')] <- 'NL1&NL2'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV70')] <- 'NL1&NL2'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV71')] <- 'NL1&NL2'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV31')] <- 'NL1&NL2'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV57')] <- 'NL1&NL2'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV13')] <- 'NL1&NL2'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV83')] <- 'NL1&NL2'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV6')] <- 'NL1&NL2'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV142')] <- 'NL1&NL2'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV232')] <- 'NL1&NL2'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV9')] <- 'NL1&NL2'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV38')] <- 'NL1&NL2'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV330')] <- 'NL1&NL2'

df_prev$nutrient_marker[which(df_prev$Id == 'ASV233')] <- 'NL3'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV67')] <- 'NL3'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV155')] <- 'NL3'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV136')] <- 'NL3'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV10')] <- 'NL3'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV63')] <- 'NL3'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV19')] <- 'NL3'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV30')] <- 'NL3'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV2')] <- 'NL3'

df_prev$nutrient_marker[which(df_prev$Id == 'ASV19')] <- 'NL4'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV52')] <- 'NL4'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV574')] <- 'NL4'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV14')] <- 'NL4'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV43')] <- 'NL4'

df_prev$nutrient_marker[which(df_prev$Id == 'ASV7')] <- 'NL5'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV51')] <- 'NL5'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV10')] <- 'NL5'
df_prev$nutrient_marker[which(df_prev$Id == 'ASV63')] <- 'NL5'

######define the color
paleta_nutrient <- c('#67a9cf','#3690c0','#02818a','#016450', 'grey')
names(paleta_nutrient) <- c('NL1&NL2','NL3','NL4','NL5', 'Unselected')

######organize the Leaf marker
df_prev$nutrient_marker <- df_prev$nutrient_marker %>%
  factor(levels = c('NL1&NL2','NL3','NL4','NL5', 'Unselected'))

p3 <- ggplot(data =df_prev,aes(Id,prevalence)) +
  geom_bar(stat = "identity",aes(fill = nutrient_marker)) +
  #geom_hline(yintercept = 0.15,linetype = "dashed",size = 0.5 , color = "#ef6548") +
  coord_flip() +
  scale_y_continuous(breaks = seq(0,1, 0.05),expand = c(0,0), 
                     labels = scales::percent) +
  scale_fill_manual(name = 'Leaf marker',values = paleta_nutrient) +
  theme_ohchibi() + 
  scale_x_discrete(expand = c(0,0)) +
  theme(axis.title.y = element_blank(),
    legend.position = "right",
    axis.text.y = element_text(size = 7,
                               color = c('red','black','red','red','black',
                                         'black','red','black','red','black',
                                         'red','black','red','red','red',
                                         'black','black','black','black','black',
                                         'red','red','red','red','red',
                                         'black','red','red','black','black',
                                         'red','black','black','black','black',
                                         'black','black','black','black','black',
                                         'red','red','red','red','black',
                                         'red','black','black','black','black',
                                         'black','black','red','black','red',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'red','black','black','black','black',
                                         'black','black','red','red','black',
                                         'black','red','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black','black','black','black','black',
                                         'black', 'black', 'black', 'black', 'black',
                                         'black', 'black') %>% rev),
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 9))  +
  ylab(label = "Prevalence ASV")

p3

######composition
composition_nutrient <- egg::ggarrange(p7, p3, nrow = 1)

######save the plot
oh.save.pdf(p = p7,
            outname = "figS3I_leaf_nutrient_marker_no_prevalence.pdf",
            outdir = "../figures/",width = 8,height = 14)
