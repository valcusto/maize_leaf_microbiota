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
df <- read.csv('../rawdata/leaf_lenght_august_dropout.csv')

###rename the full hoagland and 1/4 hoagland
df$nutrient[which(df$nutrient == 'Full Hoagland')] <- 'High'
df$nutrient[which(df$nutrient == '1/4 Hoagland')] <- 'Low'

#######organize the order 
df$SynCom <- df$SynCom %>%
  factor(levels = c('NoBacteria','Fullsyncom','L1&L2',
                    'L3','L4','L5'))

######select only low nutrient
df_low <- df %>% subset(nutrient == 'Low' &
                        NumLeaf != 'L6') %>% droplevels

df_low_sub <- df_low %>% subset(UiD != 130 &
                                UiD != 131 & UiD != 150 & UiD != 151 &
                                UiD != 205 & UiD != 206 & UiD != 453 &
                                UiD != 483 & UiD != 487)


######find influential points on the dataset
model <- lm(length_cm ~ SynCom, data = df_low_sub)

#find Cook's distance for each observation in the dataset
cooksD <- cooks.distance(model)

# Plot Cook's Distance with a horizontal line at 4/n to see which observations
#exceed this thresdhold
n <- nrow(df_low_sub)
plot(cooksD, main = "Cooks Distance for Influential Obs")
abline(h = 4/n, lty = 2, col = "steelblue") # add cutoff line

#identify influential points
influential_obs <- as.numeric(names(cooksD)[(cooksD > (4/n))])

#define new data frame with influential points removed
df_low_removed <- df_low_sub[-influential_obs, ]

######explore the data
ggplot(data = df_low_removed, aes(x=SynCom, y=length_cm, color = SynCom))+
  geom_boxplot()+
  geom_point()+
  facet_grid(.~NumLeaf)

###1.interdependence of the variables
###2.variance homogeneity
leveneTest(length_cm ~ SynCom, data = df_low_removed)
#results: p-value>0.05. Do not reject the null hypothesis.
#All the groups ahave similar variance

###3.normality
qqnorm(df_low_removed$length_cm)
qqline(df_low_removed$length_cm, col='red')
hist(df_low_removed$length_cm)

######remove NA from the dataset
df_low <- df_low_removed %>% subset(length_cm != 'NA')

######perform kruskal walis for each leaves
l1 <- df_low %>% subset(NumLeaf == 'L1') %>% droplevels

######apply kruskal-wallis test
m1_l1 <- kruskal.test(length_cm ~ SynCom,
                      data = l1)

m1_l1$p.value

DT_l1 = dunnTest(length_cm ~ SynCom,data=l1,method="bh")

PT_l1 = DT_l1$res

PT_l1

df_letter_l1 <- cldList(P.adj ~ Comparison,
                        data = PT_l1,
                        threshold = 0.1)

df_letter_l1$NumLeaf <- 'L1'

######perform kruskal walis for each leaves
l2 <- df_low %>% subset(NumLeaf == 'L2') %>% droplevels

######apply kruskal-wallis test
m1_l2 <- kruskal.test(length_cm ~ SynCom,
                      data = l2)

m1_l2$p.value

DT_l2 = dunnTest(length_cm ~ SynCom,data=l2,method="bh")

PT_l2 = DT_l2$res

PT_l2

df_letter_l2 <- cldList(P.adj ~ Comparison,
                        data = PT_l2,
                        threshold = 0.1)

df_letter_l2$NumLeaf <- 'L2'

###########subset leaf 3
l3 <- df_low %>% subset(NumLeaf == 'L3') %>% droplevels

######apply kruskal-wallis test
m1_l3 <- kruskal.test(length_cm ~ SynCom,
                      data = l3)

DT_l3 = dunnTest(length_cm ~ SynCom,data=l3,method="bh")

PT_l3 = DT_l3$res

PT_l3

df_letter_l3 <- cldList(P.adj ~ Comparison,
                        data = PT_l3,
                        threshold = 0.1)

df_letter_l3$NumLeaf <- 'L3'

######change the MonoLetter order
df_letter_l3$MonoLetter[which(df_letter_l3$Letter == 'bc')] <- 'ab'
df_letter_l3$MonoLetter[which(df_letter_l3$Letter == 'a')] <- 'c'
df_letter_l3$MonoLetter[which(df_letter_l3$Letter == 'ab')] <- 'bc'
df_letter_l3$MonoLetter[which(df_letter_l3$Letter == 'c')] <- 'b'

###########subset leaf 4
l4 <- df_low %>% subset(NumLeaf == 'L4') %>% droplevels

######apply kruskal-wallis test
m1_l4 <- kruskal.test(length_cm ~ SynCom,
                      data = l4)

DT_l4 = dunnTest(length_cm ~ SynCom,data=l4,method="bh")

PT_l4 = DT_l4$res

PT_l4

df_letter_l4 <- cldList(P.adj ~ Comparison,
                        data = PT_l4,
                        threshold = 0.1)

df_letter_l4$NumLeaf <- 'L4'

######change the MonoLetter order
df_letter_l4$MonoLetter[which(df_letter_l4$Letter == 'c')] <- 'a'
df_letter_l4$MonoLetter[which(df_letter_l4$Letter == 'ab')] <- 'bc'
df_letter_l4$MonoLetter[which(df_letter_l4$Letter == 'a')] <- 'c'
df_letter_l4$MonoLetter[which(df_letter_l4$Letter == 'bc')] <- 'ac'

###########subset leaf 5
l5 <- df_low %>% subset(NumLeaf == 'L5') %>% droplevels

######apply kruskal-wallis test
m1_l5 <- kruskal.test(length_cm ~ SynCom,
                      data = l5)

DT_l5 = dunnTest(length_cm ~ SynCom,data=l5,method="bh")

PT_l5 = DT_l5$res

PT_l5

df_letter_l5 <- cldList(P.adj ~ Comparison,
                        data = PT_l5,
                        threshold = 0.1)

df_letter_l5$NumLeaf <- 'L5'

######change the MonoLetter order
#df_letter_l5$MonoLetter[which(df_letter_l5$Letter == 'a')] <- 'b'
#df_letter_l5$MonoLetter[which(df_letter_l5$Letter == 'b')] <- 'a'

######combine the data sets
df_letter <- rbind(df_letter_l1, df_letter_l2, df_letter_l3,
                   df_letter_l4, df_letter_l5)

######calculate the mean and the confidence interval for each analysis
###create the unique names for the for loop
mtreat <- df_low$SynCom %>% unique
mleaf <- df_low$NumLeaf %>% unique
###Empty dataset
Res_mean <- NULL
Res_ci <- NULL

###perform the for loop
for (leaf in mleaf){
  df_leaf <- df %>% subset(NumLeaf == leaf) %>% droplevels
  for (treat in mtreat) {
    df_treat <- df_leaf %>% subset(SynCom == treat) %>%droplevels
    #calculate means
    df_mean <- mean(df_treat$length_cm) %>% as.data.frame
    ###calculate ci
    # bootstrapping with 1000 replications 
    results <- boot(data=df_treat$length_cm, statistic=Bmean, R=1000)
    # get 95% confidence interval 
    results_boot <- boot.ci(results, conf = 0.90, 
                            type=c("norm", "basic", "perc", "bca"))
    ###create a data frame for ci
    df_ci <- results_boot$bca %>% as.data.frame
    #add the leaf name on data frame
    df_mean$NumLeaf <- leaf
    df_mean$SynCom <- treat
    df_ci$NumLeaf <- leaf
    df_ci$SynCom <- treat
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
res <- merge(Res_mean, Res_ci, by=c('NumLeaf','SynCom'))

######slipt the group names in df_letter
colnames(df_letter)[1] <- 'SynCom'

######merge res and letter
df_res <- merge(res, df_letter, by=c('NumLeaf','SynCom'))

#######change the SynCom names
df_low$SynCom <- df_low$SynCom %>% as.character
df_low$SynCom[which(df_low$SynCom == 'NoBacteria')] <- 'NB'
df_low$SynCom[which(df_low$SynCom == 'Fullsyncom')] <- '+SynCom'
df_low$SynCom[which(df_low$SynCom == 'L1&L2')] <- "+SynCom-L1&L2"
df_low$SynCom[which(df_low$SynCom == 'L3')] <- "+SynCom-L3"
df_low$SynCom[which(df_low$SynCom == 'L4')] <- "+SynCom-L4"
df_low$SynCom[which(df_low$SynCom == 'L5')] <- "+SynCom-L5"

df_res$SynCom <- df_res$SynCom %>% as.character
df_res$SynCom[which(df_res$SynCom == 'NoBacteria')] <- 'NB'
df_res$SynCom[which(df_res$SynCom == 'Fullsyncom')] <- '+SynCom'
df_res$SynCom[which(df_res$SynCom == 'L1&L2')] <- "+SynCom-L1&L2"
df_res$SynCom[which(df_res$SynCom == 'L3')] <- "+SynCom-L3"
df_res$SynCom[which(df_res$SynCom == 'L4')] <- "+SynCom-L4"
df_res$SynCom[which(df_res$SynCom == 'L5')] <- "+SynCom-L5"

######create color for the dataset
paleta_syncom <- c('#d0d1e6','#a6bddb','#67a9cf',
                   '#3690c0','#02818a','#016450')
names(paleta_syncom) <- c('NB','+SynCom','+SynCom-L1&L2',
                          '+SynCom-L3','+SynCom-L4','+SynCom-L5')

#######organize the order 
df_low$SynCom <- df_low$SynCom %>%
  factor(levels = c('NB','+SynCom','+SynCom-L1&L2',
                    '+SynCom-L3','+SynCom-L4','+SynCom-L5'))

df_res$SynCom <- df_res$SynCom %>%
  factor(levels = c('NB','+SynCom','+SynCom-L1&L2',
                    '+SynCom-L3','+SynCom-L4','+SynCom-L5'))

######plot the results
p_leaf <- ggplot(data=df_low, aes(x=SynCom, y=length_cm,
                                 color=SynCom))+
  geom_point(shape=21, size=1.5, alpha=1)+
  facet_grid(.~NumLeaf)+
  geom_pointrange(data=df_res, aes(y=mean, ymin=lower, ymax=upper, 
                               color=SynCom), 
                  size=1)+
  geom_text(data = df_res, aes(x = SynCom,y = 40,
                           label = MonoLetter),
            inherit.aes = F,size = 5,family ="Arial",color = "black") + 
  ylab('Leaf length (cm)')+
  scale_color_manual(values = paleta_syncom)+
  clean +
  theme(axis.text.x = element_text(size=12,vjust=1, angle = 45, hjust=1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, face='bold'),
        axis.title.x = element_blank(),
        legend.position = 'none')

p_leaf

######save the figures
oh.save.pdf(p = p_leaf,
            outname = "fig3C_leaf_lenght_dropout_final.pdf",
            outdir = "../figures/",width = 12,height = 8)


ggsave(plot = p_leaf,'../figure/fig3D_1_leaf_lenght_dropout.png',
       units='cm',
       width=40, height=20)
