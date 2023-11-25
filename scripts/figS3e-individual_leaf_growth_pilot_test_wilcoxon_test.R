######load packages
library(ggplot2)
library(scales)
library(reshape2)
library(ohchibi)
library(gridExtra)
library(dplyr)
library(rcompanion)
library(emmeans)
library(multcomp)
library(car)
library(FSA)
library(boot)
library(tidyverse)

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
df <-  read.csv("../cleandata/leaf_length_stem_pilot_test_final_test.csv", header = T, sep = ",",
                quote = "",comment.char = "")

######transform the factor
df$treatment <- df$treatment %>% factor(levels = c('Nonutrient','quarter',
                                                   'half','Full'))
df$num_leaf <- df$num_leaf %>% factor

######create color treatment
paleta_new <- c('#c7eae5','#80cdc1','#f6e8c3','#dfc27d')
names(paleta_new) <- c('No\nnutrient', '1/4\nHoagland', 
                         '1/2\nHoagland', 'Full\nHoagland')

######paleta leaf
paleta_leaf <- c('#d01c8b','#f1b6da','#7b3294',
                 '#c2a5cf','#a6dba0','#008837')
names(paleta_leaf) <- c('L1', 'L2', 'L3', 'L4', 'L5', 'L6')

######create a new column group with the num_leaf and treatment
df$group <- paste(df$treatment, df$num_leaf, sep='_')

######statistical assumptions
#1.independence of the variables

#2.variance homogeneity
#perform levene test
leveneTest(length_stem ~ group, data = df)
#p<0.05. Reject H0.
#Not similar variance
#3.check the normality of the data
qqnorm(df$length_stem)
qqline(df$length_stem, col='red')

##Data exploration
ggplot(data = df,mapping = aes(x = num_leaf,y = length_stem)) +
  geom_boxplot() + 
  geom_point()+
  facet_grid(.~treatment)


#######apply the kruskal wallis for each leaves
l1 <- df %>% subset(num_leaf == 'L1') %>% droplevels

######calculate the median and the IQR for L1
group_by(l1, treatment) %>%
  summarise(
    count = n(),
    median = median(length_stem, na.rm = TRUE),
    IQR = IQR(length_stem, na.rm = TRUE)
  )

######define the groups
nn <- l1 %>% subset(treatment == 'Nonutrient') %>% droplevels
quarter <- l1 %>% subset(treatment == 'quarter') %>% droplevels
half <- l1 %>% subset(treatment == 'half') %>% droplevels
full <- l1 %>% subset(treatment == 'full') %>% droplevels

######apply kruskal-wallis test
m1_l1 <- kruskal.test(length_stem ~ treatment,
                      data = l1)

DT_l1 = dunnTest(length_stem ~ treatment,data=l1,method="bh")

PT_l1 = DT_l1$res

PT_l1

df_letter_l1 <- cldList(P.adj ~ Comparison,
                        data = PT_l1,
                        threshold = 0.05)

df_letter_l1$num_leaf <- 'L1'

######L2
l2 <- df %>% subset(num_leaf == 'L2') %>% droplevels

######apply kruskal-wallis test
m1_l2 <- kruskal.test(length_stem ~ treatment,
                      data = l2)

DT_l2 = dunnTest(length_stem ~ treatment,data=l2,method="bh")

PT_l2 = DT_l2$res

PT_l2

df_letter_l2 <- cldList(P.adj ~ Comparison,
                        data = PT_l2,
                        threshold = 0.05)

df_letter_l2$num_leaf <- 'L2'

######L3
l3 <- df %>% subset(num_leaf == 'L3') %>% droplevels

######apply kruskal-wallis test
m1_l3 <- kruskal.test(length_stem ~ treatment,
                      data = l3)

DT_l3 = dunnTest(length_stem ~ treatment,data=l3,method="bh")

PT_l3 = DT_l3$res

PT_l3

df_letter_l3 <- cldList(P.adj ~ Comparison,
                        data = PT_l3,
                        threshold = 0.1)

df_letter_l3$num_leaf <- 'L3'

######L4
l4 <- df %>% subset(num_leaf == 'L4') %>% droplevels


nn <- l4 %>% subset(treatment == 'Nonutrient') %>% droplevels
quarter <- l4 %>% subset(treatment == 'quarter') %>% droplevels
half <- l4 %>% subset(treatment == 'half') %>% droplevels
full <- l4 %>% subset(treatment == 'Full') %>% droplevels

res <- wilcox.test(full$length_stem, quarter$length_stem)
res

m1_l4 <- kruskal.test(length_stem ~ treatment,
                      data = l4)

DT_l4 = dunnTest(length_stem ~ treatment,data=l4,method="bh")

PT_l4 = DT_l4$res

PT_l4

df_letter_l4 <- cldList(P.adj ~ Comparison,
                        data = PT_l4,
                        threshold = 0.1)

df_letter_l4$num_leaf <- 'L4'

######L5
l5 <- df %>% subset(num_leaf == 'L5') %>% droplevels

######apply kruskal-wallis test
m1_l5 <- kruskal.test(length_stem ~ treatment,
                      data = l5)

DT_l5 = dunnTest(length_stem ~ treatment,data=l5,method="bh")

PT_l5 = DT_l5$res

PT_l5

df_letter_l5 <- cldList(P.adj ~ Comparison,
                        data = PT_l5,
                        threshold = 0.05)

df_letter_l5$num_leaf <- 'L5'

######combine the dataset
df_letter <- rbind(df_letter_l1, df_letter_l2, df_letter_l3,
                   df_letter_l4, df_letter_l5)

######create the for loop
######calculate the mean and the confidence interval for each analysis
###create the unique names for the for loop
mtreat <- df$group %>% unique
###Empty dataset
Res_mean <- NULL
Res_ci <- NULL

###perform the for loop
for (treat in mtreat){
  df_treat <- df %>% subset(group == treat) %>% droplevels
  #calculate means
  df_mean <- mean(df_treat$length_stem) %>% as.data.frame
  ###calculate ci
  # bootstrapping with 1000 replications 
  results <- boot(data=df_treat$length_stem, statistic=Bmean, R=1000)
  # get 95% confidence interval 
  results_boot <- boot.ci(results, type=c("norm", "basic", "perc", "bca"))
  ###create a data frame for ci
  df_ci <- results_boot$bca %>% as.data.frame
  #add the leaf name on data frame
  df_mean$group <- treat
  df_ci$group <- treat
  #combine the results
  Res_mean <- rbind(Res_mean, df_mean)
  Res_ci <- rbind(Res_ci, df_ci)
}

######create the dataset formean
colnames(Res_mean)[1]<-'mean'
###select column from 4to6
Res_ci <- Res_ci[,4:6]
###rename the first column
colnames(Res_ci)[1]<-'lower'
colnames(Res_ci)[2]<-'upper'

######merge the two datasets
res <- merge(Res_mean, Res_ci, by='group')

######remove res between group names
res$group <- res$group %>% gsub(pattern = ' ', replacement = '')

######slipt the group names in df_letter
colnames(df_letter)[1] <- 'treatment'
df_letter$treatment <- df_letter$treatment %>% gsub(pattern = ' ', replacement = '')
######prepare a new column with the location information
res$treatment <- sapply(strsplit(as.character(res$group),
                                  '_'),`[`, 1)
res$num_leaf <- sapply(strsplit(as.character(res$group),
                                      '_'),`[`, 2)

######merge res and letter
df_letter <- merge(res, df_letter, by=c('treatment', 'num_leaf'))

######plot the results
p <- ggplot(data =df,aes(treatment,length_stem, color = treatment)) +
  geom_point(shape=21, linewidth=4, alpha=10)+
  geom_pointrange(data = df_letter,
                  aes(y = mean, ymin = lower, ymax = upper),
                  size = 1.5)+
  facet_grid(.~num_leaf)

p

######edit the treatment name
df$new <- ''
df$new[which(df$treatment=='Nonutrient')] <- 'No\nnutrient'
df$new[which(df$treatment=='quarter')] <- '1/4\nHoagland'
df$new[which(df$treatment=='half')] <- '1/2\nHoagland'
df$new[which(df$treatment=='Full')] <- 'Full\nHoagland'

df_letter$new <- ''
df_letter$new[which(df_letter$treatment=='Nonutrient')] <- 'No\nnutrient'
df_letter$new[which(df_letter$treatment=='quarter')] <- '1/4\nHoagland'
df_letter$new[which(df_letter$treatment=='half')] <- '1/2\nHoagland'
df_letter$new[which(df_letter$treatment=='Full')] <- 'Full\nHoagland'

######order the new position
df$new <- df$new %>% factor(levels = c('No\nnutrient', '1/4\nHoagland',
                                               '1/2\nHoagland', 'Full\nHoagland'))

######plot the result
p_leaf <- ggplot(data=df, aes(x=new, y=length_stem,color=new))+
  geom_point(shape=21, size=4, alpha=30)+
  geom_pointrange(data=df_letter, aes(y=mean, ymin=lower, ymax=upper), 
                  size=1.5)+
  geom_text(data = df_letter, aes(x = new,y = 60,
                           label = Letter),
            inherit.aes = F,size = 6,family ="Arial",color = "black") +
  scale_color_manual(values = paleta_new)+
  xlab(NULL)+ 
  ylab('Leaf length (cm)')+
  ylim(0,60)+
  facet_grid(.~num_leaf)+
  clean +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 30),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size=30, face = 'bold'))

p_leaf

######save the figures
oh.save.pdf(p = p_leaf,
            outname = "figS3C-individual_leaf_lenght_pilot_final.pdf",
            outdir = "../figures/",width = 30,height = 10)
