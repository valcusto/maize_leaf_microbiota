######Load packages
library(ohchibi)
library(multcomp)
library(emmeans)
library(glue)
library(ggtext)

######Set seed
set.seed(130816)

######Set work directory
setwd("C:/Users/Castrillo Lab/OneDrive - The University of Nottingham/Publication 2019/Microbiome/Mutants_screening/scripts")
source("0-Clean_up_plots.R")
######create color for the dataset
paleta_syncom <- c('#80cdc1','#018571')
names(paleta_syncom) <- c('NB','+SynCom')

######Load dataset
df_400 <- read.csv('../rawdata/arabidopsis_mutants_1_400_final_full_nb.csv')

#####Data preparation
###Add growth media dilution column
df_400$dilution <- "1/400"

###Prepare columns with the treatment information and sample_id
df_400$sample_id <- sapply(strsplit(as.character(df_400$Label), "_"),`[`, 1)

df_400$treatment <- sapply(strsplit((df_400$Label), '_'), `[`, 2)

###Rename the factors in the column
df_400$treatment[which(df_400$treatment == "Full")] <- "+SynCom"
df_400$treatment[which(df_400$treatment == "full")] <- "+SynCom"
df_400$treatment[which(df_400$treatment == "Full ")] <- "+SynCom"
df_400$treatment[which(df_400$treatment == "MgCl2 ")] <- "NB"
df_400$treatment[which(df_400$treatment == "MgCl2")] <- "NB"

df_400$sample_id[which(df_400$sample_id == "Col")] <- "Col-0"

###Calculate the area in mm2
df_400$area_mm2 <- df_400$Area_cm2*1000

######Explore the data
ggplot(data = df_400,aes(x = tair,y =area_mm2)) +
  geom_boxplot(aes(color=treatment)) + 
  geom_point(aes(color=treatment), 
             position=position_dodge(width = 0.8)) + 
  facet_grid(.~dilution)+
  theme_ohchibi()

#####Explore the data only for col-0
df_col <- df_400 %>%
  subset(sample_id == "Col-0") %>% droplevels

ggplot(data = df_col,aes(x = sample_id,y =area_mm2)) +
  geom_boxplot(aes(color=treatment)) + 
  geom_point(aes(color=treatment), 
             position=position_dodge(width = 0.75)) + 
  facet_grid(.~dilution)+
  theme_ohchibi()

######Reorder the tair
df_400$tair <- df_400$tair %>%
  factor(levels=c("Col-0_+SynCom","Col-0_NB", "AT5G19200","AT2G45910","AT5G10770",
                  "AT3G58490","AT2G30250","AT5G38560","AT1G66160",
                  "AT5G24030","AT4G11650","AT4G20110",
                  "AT4G27970","AT4G17950","AT4G18030","AT2G45830",
                  "AT3G09030","AT3G48980","AT5G11230","AT1G48880",
                  "AT3G61270","AT4G19080",
                  "AT3G05050","AT2G17610","AT3G02840","AT1G11360",
                  "AT5G07100","AT5G54430","AT2G02070","AT5G22380",
                  "AT5G25400","AT5G07120",
                  "AT1G63420","AT1G62020","AT4G24970",
                  "AT5G23850","AT4G32390","AT5G58440","AT2G34940",
                  "AT3G07480"))

######Calculate emmeans
###Marginal means are basically means extracted 
#from a statistical model, and represent average of 
#response variable (here, area) 
#for each level of predictor variable (here, sample_id).
###Create a unique id to add in the model
df_full <- df_400 %>%
  subset((treatment == "+SynCom") | (sample_id == 'Col-0')) %>% droplevels

######remove the sapce on the label
df_full$Label <- df_full$Label %>%
  gsub(pattern=' ', replacement='')

######group the label
df_full$tair <- df_full$tair %>%
  factor(levels = c("Col-0_+SynCom","Col-0_NB", "AT5G19200","AT2G45910","AT5G10770",
                    "AT3G58490","AT2G30250","AT5G38560","AT1G66160",
                    "AT5G24030","AT4G11650","AT4G20110",
                    "AT4G27970","AT4G17950","AT4G18030","AT2G45830",
                    "AT3G09030","AT3G48980","AT5G11230","AT1G48880",
                    "AT3G61270","AT4G19080",
                    "AT3G05050","AT2G17610","AT3G02840","AT1G11360",
                    "AT5G07100","AT5G54430","AT2G02070","AT5G22380",
                    "AT5G25400","AT5G07120",
                    "AT1G63420","AT1G62020","AT4G24970",
                    "AT5G23850","AT4G32390","AT5G58440","AT2G34940",
                    "AT3G07480"))


m1 <- lm(formula = area_mm2 ~ tair,
         data = df_full)

summary(m1)

msum <- glht(m1, linfct=mcp(tair="Dunnett"), alternative = "greater") %>% summary 
df_em <- emmeans(object = m1,specs = "tair") %>% as.data.frame
df_em$p.value <-  c(1,msum$test$pvalues)

######Reorder the sample_id
df_em$stair <- df_em$tair %>%
  factor(levels=c("Col-0_+SynCom","Col-0_NB", "AT5G19200","AT2G45910","AT5G10770",
                  "AT3G58490","AT2G30250","AT5G38560","AT1G66160",
                  "AT5G24030","AT4G11650","AT4G20110",
                  "AT4G27970","AT4G17950","AT4G18030","AT2G45830",
                  "AT3G09030","AT3G48980","AT5G11230","AT1G48880",
                  "AT3G61270","AT4G19080",
                  "AT3G05050","AT2G17610","AT3G02840","AT1G11360",
                  "AT5G07100","AT5G54430","AT2G02070","AT5G22380",
                  "AT5G25400","AT5G07120",
                  "AT1G63420","AT1G62020","AT4G24970",
                  "AT5G23850","AT4G32390","AT5G58440","AT2G34940",
                  "AT3G07480"))

######Prepare the significnce letter
alpha <- 0.1
df_em$Significance <- "NS"
df_em$Significance[which(df_em$p.value<alpha)] <- "q<0.1"
df_em$Letters <-""
df_em$Letters[which(df_em$Significance == "q<0.1")] <- "*"

######add treatment column at df_em
df_em$treatment <- '+SynCom'
df_em$treatment[which(df_em$tair == 'Col-0_NB')] <- 'NB'

######Prepare a dataset with the values only for full SynCom Col0
df_em_full_col <- df_em %>%
  subset(tair == "Col-0_+SynCom") %>% droplevels

df_em_nb_col <- df_em %>%
  subset(tair == "Col-0_NB") %>% droplevels

df_em_col <- rbind(df_em_full_col, df_em_nb_col)

######transform treatment as factor
df_full$treatment <- df_full$treatment %>%
  factor(levels = c('NB', '+SynCom'))

######order the tair
df_full$tair <- df_full$tair %>%
  factor(levels=c("Col-0_NB","Col-0_+SynCom", "AT5G19200","AT2G45910","AT5G10770",
                  "AT3G58490","AT2G30250","AT5G38560","AT1G66160",
                  "AT5G24030","AT4G11650","AT4G20110",
                  "AT4G27970","AT4G17950","AT4G18030","AT2G45830",
                  "AT3G09030","AT3G48980","AT5G11230","AT1G48880",
                  "AT3G61270","AT4G19080",
                  "AT3G05050","AT2G17610","AT3G02840","AT1G11360",
                  "AT5G07100","AT5G54430","AT2G02070","AT5G22380",
                  "AT5G25400","AT5G07120",
                  "AT1G63420","AT1G62020","AT4G24970",
                  "AT5G23850","AT4G32390","AT5G58440","AT2G34940",
                  "AT3G07480"))

df_em$tair <- df_em$tair %>%
  factor(levels=c("Col-0_NB","Col-0_+SynCom","AT5G19200","AT2G45910","AT5G10770",
                  "AT3G58490","AT2G30250","AT5G38560","AT1G66160",
                  "AT5G24030","AT4G11650","AT4G20110",
                  "AT4G27970","AT4G17950","AT4G18030","AT2G45830",
                  "AT3G09030","AT3G48980","AT5G11230","AT1G48880",
                  "AT3G61270","AT4G19080",
                  "AT3G05050","AT2G17610","AT3G02840","AT1G11360",
                  "AT5G07100","AT5G54430","AT2G02070","AT5G22380",
                  "AT5G25400","AT5G07120",
                  "AT1G63420","AT1G62020","AT4G24970",
                  "AT5G23850","AT4G32390","AT5G58440","AT2G34940",
                  "AT3G07480"))

paleta_iqr <- c('#80cdc1','#018571')
names(paleta_iqr) <- c('Col-0_NB','Col-0_+SynCom')

######plot the result
p_400 <- ggplot(data = df_full,aes(x = tair,y =area_mm2, color=treatment)) +
  geom_point(shape = 21,size = 1, aes(group=treatment), alpha=0.5)+
  geom_text(data = df_em, aes(x = tair,y = 190,label = Letters),
            inherit.aes = F,size = 6,family ="Arial",color = "black") +
  geom_pointrange(
    data = df_em,
    aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
    size = 0.5) +
  scale_color_manual(name = 'Treatment', values = paleta_syncom)+
  geom_rect(data = df_em_full_col,
            mapping = aes(ymin = lower.CL,ymax = upper.CL,xmin = 0,
                          xmax = Inf, fill = tair),
            inherit.aes = F,alpha = 0.2,color = NA,
            alpha =0.3 )+
  scale_fill_manual(name='IQR',
                    values = paleta_iqr,  
                    guide = guide_legend(override.aes = list(alpha = 0.2)))+
  ggtitle("1/400 MS")+
  xlab(element_blank())+
  ylab(label = "Rosette area (mm2)") + clean +
  theme(axis.text.x = element_text(size = 15,angle=45, hjust=1,
                                   color = c('black', 'black', 
                                             'black', 'black',
                                             'red', 'black',
                                             'red', 'black',
                                             'black', 'black',
                                             'black', 'black',
                                             'black', 'black',
                                             'black', 'black',
                                             'black', 'black',
                                             'black', 'black',
                                             'black', 'black',
                                             'black', 'black',
                                             'black', 'black',
                                             'black', 'black',
                                             'black', 'black',
                                             'black', 'red',
                                             'black', 'black',
                                             'black', 'black',
                                             'black', 'black',
                                             'black', 'black')),
        plot.title = element_text(size=30, hjust=0.5))

p_400

######save the data
oh.save.pdf(p = p_400,
            outname = "fig5E_arabidopsis_screening_1_400_without_NB_shade.pdf",
            outdir = "../figures/",width = 20,height = 10)
