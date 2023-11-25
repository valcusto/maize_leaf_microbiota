######load packages
library(ohchibi)
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
df <- read.csv('../rawdata/figure3_leaf_length.csv')

#####################data preparation
###rename the first column
colnames(df)[1] <- 'pot'

###rename the treatment column
df$Treatment[which(df$Treatment == 'Full syncom')] <- '+SynCom'
df$Treatment[which(df$Treatment == 'No Bacteria')] <- 'NB'

###order the Treatment column
df$Treatment <- df$Treatment %>% factor(levels = c('NB', '+SynCom'))

###rename the full hoagland and 1/4 hoagland
df$Nutrient[which(df$Nutrient == 'Full Hoagland')] <- 'High'
df$Nutrient[which(df$Nutrient == '1/4 Hoagland')] <- 'Low'

###order the nutrient 
df$Nutrient <- df$Nutrient %>% factor(levels = c('Low','High'))

######explore the data
ggplot(data = df, aes(x=Nutrient, y=leaf_lenght_cm, color = Treatment))+
  geom_boxplot()+
  geom_point()+
  facet_grid(.~NumLeaf)

######remove the outliers
df <- df[-c(42,116,199,209,254,374),]

######statistical assumptions
###create a group with the treatment and nutrient
df$group <- paste(df$Treatment, df$Nutrient, sep='_')

###transform group as a factor
df$group <- df$group %>% factor (levels = c('NB_Low','+SynCom_Low',
                                            'NB_High','+SynCom_High'))

######find influential points on the dataset
model <- lm(leaf_lenght_cm ~ group, data = df)

#find Cook's distance for each observation in the dataset
cooksD <- cooks.distance(model)

# Plot Cook's Distance with a horizontal line at 4/n to see which observations
#exceed this thresdhold
n <- nrow(df)
plot(cooksD, main = "Cooks Distance for Influential Obs")
abline(h = 4/n, lty = 2, col = "steelblue") # add cutoff line

#identify influential points
influential_obs <- as.numeric(names(cooksD)[(cooksD > (4/n))])

#define new data frame with influential points removed
df_removed <- df[-influential_obs, ]

###1.interdependence of the variables
###2.variance homogeneity
leveneTest(leaf_lenght_cm ~ group, data = df_removed)
#results: p-value<0.05. Reject the null hypothesis.
#All the groups ahave similar variance

###3.normality
qqnorm(df_removed$leaf_lenght_cm)
qqline(df_removed$leaf_lenght_cm, col='red')
hist(df_removed$leaf_lenght_cm)

######remove L6 for the analysis
df <- df_removed %>% subset(NumLeaf != 'L6') %>% droplevels

######remove NA from the dataset
df <- df %>% subset(leaf_lenght_cm != 'NA')

######perform kruskal walis for each leaves
l1 <- df %>% subset(NumLeaf == 'L1') %>% droplevels

######apply kruskal-wallis test
m1_l1 <- kruskal.test(leaf_lenght_cm ~ group,
                      data = l1)

m1_l1$p.value

DT_l1 = dunnTest(leaf_lenght_cm ~ group,data=l1,method="bh")

PT_l1 = DT_l1$res

PT_l1

df_letter_l1 <- cldList(P.adj ~ Comparison,
                        data = PT_l1,
                        threshold = 0.1)

df_letter_l1$NumLeaf <- 'L1'

###########subset leaf 2
l2 <- df %>% subset(NumLeaf == 'L2') %>% droplevels

######apply kruskal-wallis test
m1_l2 <- kruskal.test(leaf_lenght_cm ~ group,
                      data = l2)

DT_l2 = dunnTest(leaf_lenght_cm ~ group,data=l2,method="bh")

PT_l2 = DT_l2$res

PT_l2

df_letter_l2 <- cldList(P.adj ~ Comparison,
                        data = PT_l2,
                        threshold = 0.1)

df_letter_l2$NumLeaf <- 'L2'

###########subset leaf 3
l3 <- df %>% subset(NumLeaf == 'L3') %>% droplevels

######apply kruskal-wallis test
m1_l3 <- kruskal.test(leaf_lenght_cm ~ group,
                      data = l3)

DT_l3 = dunnTest(leaf_lenght_cm ~ group,data=l3,method="bh")

PT_l3 = DT_l3$res

PT_l3

df_letter_l3 <- cldList(P.adj ~ Comparison,
                        data = PT_l3,
                        threshold = 0.1)

df_letter_l3$NumLeaf <- 'L3'

######change the MonoLetter order
df_letter_l3$MonoLetter[which(df_letter_l3$Letter == 'c')] <- 'a'
df_letter_l3$MonoLetter[which(df_letter_l3$Letter == 'a')] <- 'c'

###########subset leaf 4
l4 <- df %>% subset(NumLeaf == 'L4') %>% droplevels

######apply kruskal-wallis test
m1_l4 <- kruskal.test(leaf_lenght_cm ~ group,
                      data = l4)

DT_l4 = dunnTest(leaf_lenght_cm ~ group,data=l4,method="bh")

PT_l4 = DT_l4$res

PT_l4

df_letter_l4 <- cldList(P.adj ~ Comparison,
                        data = PT_l4,
                        threshold = 0.1)

df_letter_l4$NumLeaf <- 'L4'

######change the MonoLetter order
df_letter_l4$MonoLetter[which(df_letter_l4$Letter == 'c')] <- 'a'
df_letter_l4$MonoLetter[which(df_letter_l4$Letter == 'a')] <- 'c'

###########subset leaf 5
l5 <- df %>% subset(NumLeaf == 'L5') %>% droplevels

######apply kruskal-wallis test
m1_l5 <- kruskal.test(leaf_lenght_cm ~ group,
                      data = l5)

DT_l5 = dunnTest(leaf_lenght_cm ~ group,data=l5,method="bh")

PT_l5 = DT_l5$res

PT_l5

df_letter_l5 <- cldList(P.adj ~ Comparison,
                        data = PT_l5,
                        threshold = 0.1)

df_letter_l5$NumLeaf <- 'L5'

######change the MonoLetter order
df_letter_l5$MonoLetter[which(df_letter_l5$Letter == 'a')] <- 'b'
df_letter_l5$MonoLetter[which(df_letter_l5$Letter == 'b')] <- 'a'

######combine the data sets
df_letter <- rbind(df_letter_l1, df_letter_l2, df_letter_l3,
                   df_letter_l4, df_letter_l5)

######calculate the mean and the confidence interval for each analysis
###create the unique names for the for loop
mtreat <- df$group %>% unique
mleaf <- df$NumLeaf %>% unique
###Empty dataset
Res_mean <- NULL
Res_ci <- NULL

###perform the for loop
for (leaf in mleaf){
  df_leaf <- df %>% subset(NumLeaf == leaf) %>% droplevels
  for (treat in mtreat) {
    df_treat <- df_leaf %>% subset(group == treat) %>%droplevels
    #calculate means
    df_mean <- mean(df_treat$leaf_lenght_cm) %>% as.data.frame
    ###calculate ci
    # bootstrapping with 1000 replications 
    results <- boot(data=df_treat$leaf_lenght_cm, statistic=Bmean, R=1000)
    # get 95% confidence interval 
    results_boot <- boot.ci(results, conf = 0.90, type=c("norm", "basic", "perc", "bca"))
    ###create a data frame for ci
    df_ci <- results_boot$bca %>% as.data.frame
    #add the leaf name on data frame
    df_mean$NumLeaf <- leaf
    df_mean$group <- treat
    df_ci$NumLeaf <- leaf
    df_ci$group <- treat
    #combine the results
    Res_mean <- rbind(Res_mean, df_mean)
    Res_ci <- rbind(Res_ci, df_ci) 
  }
}

######create the dataset formean
colnames(Res_mean)[1]<-'mean'
###select column from 4to7
Res_ci <- Res_ci[,4:7]
###rename the first column
colnames(Res_ci)[1]<-'lower'
colnames(Res_ci)[2]<-'upper'

######merge the two datasets
res <- merge(Res_mean, Res_ci, by=c('NumLeaf','group'))

######slipt the group names in df_letter
colnames(df_letter)[1] <- 'group'

######prepare a new column with the location information
df_letter$Treatment <- sapply(strsplit(as.character(df_letter$group),
                                  '_'),`[`, 1)
df_letter$Nutrient <- sapply(strsplit(as.character(df_letter$group),
                                     '_'),`[`, 2)

######merge res and letter
df_res <- merge(res, df_letter, by=c('NumLeaf','group'))

######save the df_res for the syncom selection
mlist <- list(
  df_res = df_res
)
saveRDS(object = mlist,file = "../cleandata/res_syncomleaves_complete_length.RDS")

######regroup the order
df_res$group <- df_res$group %>% factor(levels = c('NB_Low', '+SynCom_Low',
                                                   'NB_High', '+SynCom_High'))

df$group <- df$group %>% factor(levels = c('NB_Low', '+SynCom_Low',
                                           'NB_High', '+SynCom_High'))

######create color for the dataset
paleta_group <- c('#dfc27d','#a6611a','#80cdc1','#018571')
names(paleta_group) <- c('NB_High', '+SynCom_High', 'NB_Low', '+SynCom_Low')

######letters
df_letter$MonoLetter <- df_letter$MonoLetter %>% 
  gsub(pattern = ' ', replacement = '')

######plot the results
p_leaf <- ggplot(data=df, aes(x=group, y=leaf_lenght_cm,
                                 color=group))+
  geom_point(shape=21, size=1.5, alpha=1)+
  facet_grid(.~NumLeaf)+
  geom_pointrange(data=df_res, aes(y=mean, ymin=lower, ymax=upper, 
                               color=group), 
                  size=1)+
  geom_text(data = df_res, aes(x = group,y = 40,
                           label = MonoLetter),
            inherit.aes = F,size = 5,family ="Arial",color = "black") + 
  ylab('Leaf length (cm)')+
  scale_color_manual(values = paleta_group)+
  clean +
  theme(axis.text.x = element_text(size=12,vjust=1, angle = 45, hjust=1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, face='bold'),
        axis.title.x = element_blank(),
        legend.position = 'none')

p_leaf

######save the figures
oh.save.pdf(p = p_leaf,
            outname = "fig3A_syncom_single_leaf.pdf",
            outdir = "../figures/",width = 12,height = 8)
