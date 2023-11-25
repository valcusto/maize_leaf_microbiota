######load packages
library(ohchibi)
library(rcompanion)
library(emmeans)
library(multcomp)
library(car)
library(FSA)
library(boot)
library(tidyverse)
library(gridExtra)

######set seed
set.seed(12356)

######set work directory
setwd("C:/Users/Castrillo Lab/OneDrive - The University of Nottingham/Publication 2019/Paper preparation/manuscripts/scripts")

######helper functions
source('0-Clean_up_plots.R')
# function to obtain the mean
Bmean <- function(data, indices) {
  d <- data[indices] # allows boot to select sample 
  return(mean(d))
}

######upload data
df <- read.csv('../cleandata/leaf_lenght_community_context_1.csv')

#####create the total column
df$total <- df$length_cm+df$stem

###rename the full hoagland and 1/4 hoagland
df$Nutrient[which(df$Nutrient == 'Full Hoagland')] <- 'High'
df$Nutrient[which(df$Nutrient == '1/4 Hoagland')] <- 'Low'

######select only low nutrient
df_low <- df %>% subset(Nutrient == 'Low') %>% droplevels

#######change the SynCom names
df_low$Treatment <- df_low$Treatment %>% as.character
df_low$Treatment[which(df_low$Treatment == 'No Bacteria')] <- 'NB'
df_low$Treatment[which(df_low$Treatment == 'Full syncom')] <- '+SynCom'
df_low$Treatment[which((df_low$Treatment == 'SynCom_1') | (df_low$Treatment == 'SynCom_5'))] <- "+L1&L2"
df_low$Treatment[which((df_low$Treatment == 'SynCom_2') | (df_low$Treatment == 'SynCom_6'))] <- "+L3"
df_low$Treatment[which((df_low$Treatment == 'SynCom_3') | (df_low$Treatment == 'SynCom_7'))] <- "+L4"
df_low$Treatment[which((df_low$Treatment == 'SynCom_4') | (df_low$Treatment == 'SynCom_8'))] <- "+L5"

######sum the leaf lenght
df_sum <- aggregate(total ~ Treatment+Sample_ID,
                    df_low, sum)

######check the data
ggplot(data = df_sum, aes(x=Treatment, y=total, color = Treatment))+
  geom_boxplot()+
  geom_point()

###1.interdependence of the variables
###2.variance homogeneity
leveneTest(total ~ Treatment, data = df_sum)
#results: p-value>0.05. Do not Reject the null hypothesis.
#All the groups have similar variance

###3.normality
qqnorm(df_sum$total)
qqline(df_sum$total, col='red')
hist(df_sum$total)

######transform Treatment as a factor
df_sum$Treatment <- df_sum$Treatment %>% factor

######remove the outliers
df_sum_new <- df_sum[-c(20,29,32,48,52,71),]

######apply anova
m1 <- kruskal.test(total ~ Treatment,
                      data = df_sum_new)

m1$p.value

DT = dunnTest(total ~ Treatment,data=df_sum_new,method="bh")

PT = DT$res

PT

df_letter <- cldList(P.unadj ~ Comparison,
                        data = PT,
                        threshold = 0.1)

######calculate the mean and the confidence interval for each analysis
###create the unique names for the for loop
mtreat <- df_sum$Treatment %>% unique
###Empty dataset
Res_mean <- NULL
Res_ci <- NULL

###perform the for loop
for (treat in mtreat) {
    df_treat <- df_sum %>% subset(Treatment == treat) %>%droplevels
    #calculate means
    df_mean <- mean(df_treat$total) %>% as.data.frame
    ###calculate ci
    # bootstrapping with 1000 replications 
    results <- boot(data=df_treat$total, statistic=Bmean, R=1000)
    # get 95% confidence interval 
    results_boot <- boot.ci(results, conf = 0.90, 
                            type=c("norm", "basic", "perc", "bca"))
    ###create a data frame for ci
    df_ci <- results_boot$bca %>% as.data.frame
    df_mean$Treatment <- treat
    df_ci$Treatment <- treat
    #combine the results
    Res_mean <- rbind(Res_mean, df_mean)
    Res_ci <- rbind(Res_ci, df_ci) 
}

######create the dataset formean
colnames(Res_mean)[1]<-'mean'
###select column from 4to7
Res_ci <- Res_ci[,4:6]
###rename the first column
colnames(Res_ci)[1]<-'lower'
colnames(Res_ci)[2]<-'upper'

######merge the two datasets
res <- merge(Res_mean, Res_ci, by='Treatment')

######slipt the group names in df_letter
colnames(df_letter)[1] <- 'Treatment'

######merge the two objects
df_em <- merge(res, df_letter, by = 'Treatment')

######create color for the dataset
paleta_syncom <- c('#d0d1e6','#a6bddb','#67a9cf',
                   '#3690c0','#02818a','#016450')
names(paleta_syncom) <- c('NB','+SynCom','+L1&L2',
                          '+L3','+L4','+L5')

#######organize the order 
df_sum$Treatment <- df_sum$Treatment %>%
  factor(levels = c('NB','+SynCom','+L1&L2',
                    '+L3','+L4','+L5'))

df_em$Treatment <- df_em$Treatment %>%
  factor(levels = c('NB','+SynCom','+L1&L2',
                    '+L3','+L4','+L5'))

######plot the results
p_leaf <- ggplot(data=df_sum, aes(x=Treatment, y=total,
                                 color=Treatment))+
  geom_point(shape=21, size=1.5, alpha=1)+
  geom_pointrange(data=df_em, aes(y=mean,
                                  ymin=lower, ymax=upper, color= Treatment),
                  size=1)+
  geom_text(data = df_em, aes(x = Treatment,y = 150,
                           label = Letter),
            inherit.aes = F,size = 10,family ="Arial",color = "black") + 
  ylab('Leaf length (cm)')+
  #ggtitle('Small bacterial communities with stem')+
  scale_color_manual(values = paleta_syncom)+
  clean +
  theme(axis.text.x = element_text(size=20,vjust=1, angle = 45, hjust=1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, face='bold'),
        axis.title.x = element_blank(),
        title = element_text(size=20, hjust=0.5),
        legend.position = 'none')

p_leaf

######save the figures
oh.save.pdf(p = p_leaf,
            outname = "figS3N_leaf_lenght_community_context_with_stem_kruskalwallis.pdf",
            outdir = "../figures/",width = 12,height = 8)

ggsave(plot = p_leaf,'../figure/figS3N_leaf_lenght_community_context_with_stem_kruskalwallis.png',
       units='cm',
       width=40, height=20)
