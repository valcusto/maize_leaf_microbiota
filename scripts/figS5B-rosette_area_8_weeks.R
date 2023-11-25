######load library
library(ohchibi)
library(tidyverse)
library(emmeans)
library(multcomp)
library(car)
library(FSA)
library(rcompanion)
library(boot)

######set seed
set.seed(12356)

######set work directory
setwd("C:/Users/Castrillo Lab/OneDrive - The University of Nottingham/Publication 2019/Paper preparation/manuscripts/scripts")

######helper functions
source('0-Clean_up_plots.R')
######create color for the dataset
paleta_group <- c('#dfc27d','#a6611a','#80cdc1','#018571')
names(paleta_group) <- c('NB_High', '+SynCom_High', 'NB_Low', '+SynCom_Low')

######upload data
df <-  read.csv("../rawdata/arabidopsis_area_GUS_syncom_final.csv",
                header = T, sep = ",",
                quote = "",comment.char = "")

######change the nutrient name
df$nutrient[which(df$nutrient == '1/1000')] <- 'Low'
df$nutrient[which(df$nutrient == '1/2')] <- 'High'

######combine the treatment and nutrient
df$group <- paste(df$treatment, df$nutrient, sep='_')

######transform the factor
df$group <- df$group %>% factor(levels = c('NB_Low', '+SynCom_Low',
                                           'NB_High', '+SynCom_High'))

##Data exploration
ggplot(data = df,mapping = aes(x = group,y =Area)) +
  geom_boxplot() + 
  geom_point()

######check if the data are balanced
df %>% group_by(group, treatment) %>% 
  summarise(no_rows=length(group)) %>% as.data.frame

######statistical assumptions
#1.independence of the variables

#2.variance homogeneity
###levene test
leveneTest(Area ~ group, data = df)
#3.check the normality of the data
qqnorm(df$Area)
qqline(df$Area, col='red')

######model
m1 <- aov(data = df, formula = Area ~ group)
summary(m1)

msum <- glht(m1, linfct = mcp(group='Tukey')) %>% summary
df_emmeans <- emmeans(m1,specs = "group")  %>% as.data.frame

######set the letters
cld <- cld(msum)
cld <- cld$mcletters$Letters %>% as.data.frame
cld$group <- rownames(cld)
rownames(cld) <- NULL
colnames(cld)[1] <- 'letters'

######combine emmeans and cld
df_emmeans <- merge(df_emmeans, cld, by='group')

######title
mtext <- paste0("Rosette area (cm",round(m1$statistic,2))

######create the plot
p <- ggplot(data=df, aes(x=group, y=Area,
                                 color=group))+
  geom_point(shape=21, size=1.5, alpha=1)+
  geom_pointrange(data=df_emmeans, aes(y=emmean, ymin=lower.CL, 
                                       ymax=upper.CL, 
                               color=group), 
                  size=1)+
  geom_text(data = df_emmeans, aes(x = group,y = 11,
                           label = letters),
            inherit.aes = F,size = 6,family ="Arial",color = "black") +
  ylab(expression(paste("Rosette area (cm"^"2", ")")))+
  scale_color_manual(values = paleta_group)+
  clean +
  theme(axis.text.x = element_text(size=12,vjust=1, angle = 45, hjust=1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, face='bold'),
        axis.title.x = element_blank(),
        legend.position = 'none')

p

######save the plot
oh.save.pdf(p = p,
            outname = "fig5B-arabidopsis_rosette_8_weeks.pdf",
            outdir = "../figures/",width = 10,height = 10)
