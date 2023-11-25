######load library
library(ohchibi)
library(FSA)
library(rcompanion)
library(boot)
library(emmeans)
library(multcomp)
library(car)
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
dat <- readRDS(file = "../cleandata/dat_syncomleaves_complete.end.RDS")

######paleta fractions
paleta_fra <- c('#008837', '#fec44f', '#d95f0e')
names(paleta_fra) <- c('Shoot', 'Root', 'Substrate')

######paleta leaf
paleta_leaf <- c('#d01c8b','#f1b6da','#7b3294',
                 '#c2a5cf','#a6dba0','#008837', '#fec44f', '#d95f0e')
names(paleta_leaf) <- c('L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'Root', 'Substrate')

######paleta final
paleta_final <- c('#d01c8b','#f1b6da','#7b3294',
                 '#c2a5cf','#a6dba0','#008837', '#fec44f', '#d95f0e', '#41ae76')
names(paleta_final) <- c('L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'Root', 
                         'Substrate', 'Total\nLeaf')

#Rarefied
Dat_rar <- dat$Rarefied

######################## Alpha diversity
#Calculate Shannon and Richness
Dat_rar$Map$Shannon <- vegan::diversity(x = Dat_rar$Tab, 
                                        index = "shannon", MARGIN = 2)
#Richness is number of ASVs per sample
Dat_rar$Map$Richness <- colSums(Dat_rar$Tab > 0)

######select the metadata
Shan <- Dat_rar$Map

######select only syncom
Shan <- Shan %>%
  subset(Bacteria != 'NB') %>% droplevels

###### Analyse samples distribution for NumLeaf
ggplot(data = Shan,aes(x = NumLeaf,y =Shannon, color=NumLeaf)) +
  geom_boxplot() + geom_point()+
  facet_grid(.~Nutrient)+
  scale_color_manual(values = paleta_leaf)

###### Analyse samples distribution for fractions
ggplot(data = Shan,aes(x = Fraction,y =Shannon, color=Fraction)) +
  geom_boxplot() + geom_point()+
  facet_grid(.~Nutrient)+
  scale_color_manual(values = paleta_fra)

######select only shoot dataset
total <- Shan %>%
  subset(Fraction == 'Shoot') %>% droplevels

######create a column where you add total leaf
Shan$final <- Shan$NumLeaf

######add the column to the total
total$final <- 'Total\nLeaf'

######merge the two datasets
Shan <- rbind(Shan, total)

######plot the dataset
ggplot(data = Shan,aes(x = final,y =Shannon, color=final)) +
  geom_boxplot() + geom_point()+
  facet_grid(.~Nutrient)+
  scale_color_manual(values = paleta_final)

######statistical assumptions
#1.interdependence of the variables

#2.variance homogeneity
leveneTest(Shannon ~ final, data = Shan)
#p<0.05. Reject H0. At least one group have different variance.

#3.normality
qqnorm(Shan$Shannon)
qqline(Shan$Shannon, col='red')

######create a group_leaf column
Shan$group_leaf <- paste(Shan$Nutrient, Shan$final, sep='_')

######transform in factor
Shan$group_leaf <- Shan$group_leaf %>% factor

######select the low hoagland
shan_low <- Shan %>%
  subset(Nutrient == '1/4 Hoagland') %>% droplevels

######Kruskal wallis on low nutrient
m1_low <- kruskal.test(formula = Shannon ~ group_leaf,
                   data = shan_low)

m1_low$p.value

DT_low = dunnTest(Shannon ~ group_leaf,data = shan_low,method="bh")

PT_low = DT_low$res

PT_low

df_letter_low <- cldList(P.unadj ~ Comparison,
                        data = PT_low,
                        threshold = 0.1)

df_letter_low$Nutrient <- '1/4 Hoagland'

######select the high hoagland
shan_full <- Shan %>%
  subset(Nutrient != '1/4 Hoagland') %>% droplevels

######Kruskal wallis on full nutrient
m1_full <- kruskal.test(formula = Shannon ~ group_leaf,
                       data = shan_full)

m1_full$p.value

DT_full = dunnTest(Shannon ~ group_leaf,data = shan_full,method="bh")

PT_full = DT_full$res

PT_full

df_letter_full <- cldList(P.unadj ~ Comparison,
                         data = PT_full,
                         threshold = 0.1)

df_letter_full$Nutrient <- 'Full Hoagland'

######bind the rows
df_letter <- rbind(df_letter_low, df_letter_full)

######Prepare the treatment and sample_id column
df_letter$final <- sapply(strsplit(as.character(df_letter$Group), '_'), 
                        `[`, 2)

######calculate the mean and the confidence interval for each analysis
###create the unique names for the for loop
mleaf <- Shan$group_leaf %>% unique
###Empty dataset
Res_mean <- NULL
Res_ci <- NULL

###perform the for loop
for (leaf in mleaf){
  df_leaf <- Shan %>% subset(group_leaf == leaf) %>% droplevels
  #calculate means
  df_mean <- mean(df_leaf$Shannon) %>% as.data.frame
  ###calculate ci
  # bootstrapping with 1000 replications 
  results <- boot(data=df_leaf$Shannon, statistic=Bmean, R=1000)
  # get 95% confidence interval 
  results_boot <- boot.ci(results, type=c("norm", "basic", "perc", "bca"))
  ###create a data frame for ci
  df_ci <- results_boot$bca %>% as.data.frame
  #add the leaf name on data frame
  df_mean$group_leaf <- leaf
  df_ci$group_leaf <- leaf
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
res_single <- merge(Res_mean, Res_ci, by='group_leaf')

######remove res between group names
res_single$group_leaf <- res_single$group_leaf %>% 
  gsub(pattern = ' ', replacement = '')

######slipt the group names in df_letter
colnames(df_letter)[1] <- 'group_leaf'
df_letter$group_leaf <- df_letter$group_leaf %>% 
  gsub(pattern = ' ', replacement = '')

######prepare a new column with the location information
df_letter$Nutrient <- sapply(strsplit(as.character(df_letter$group_leaf),
                                  '_'),`[`, 1)

######merge res and letter
df_leaf <- merge(res_single, df_letter, by='group_leaf')
df_leaf$final <- df_leaf$final %>% 
  gsub(pattern = ' ', replacement = '')

######correct the soil names
df_leaf$Nutrient <- df_leaf$Nutrient %>%
  gsub(pattern = '1/4Hoagland', replacement = '1/4 Hoagland')

df_leaf$Nutrient <- df_leaf$Nutrient %>%
  gsub(pattern = 'FullHoagland', replacement = 'Full Hoagland')

######NumLeaf as factor
Shan$final <- Shan$final %>% 
  factor(levels = c('Substrate', 'Root','Total\nLeaf','L1', 
                    'L2', 'L3', 'L4', 'L5', 'L6'))

df_leaf$final <- df_leaf$final %>% 
  factor(levels = c('Substrate', 'Root','Total\nLeaf',
                    'L1', 'L2', 'L3', 'L4', 'L5', 'L6'))

######plot the result
p_leaf <- ggplot(data=Shan, aes(x=final, y=Shannon,
                                color=final))+
  geom_point(shape=21, size=1.5, alpha=1)+
  facet_grid(.~Nutrient)+
  geom_pointrange(data=df_leaf, aes(y=mean, ymin=lower, ymax=upper, 
                               color=final), 
                  size=1)+
  geom_text(data = df_leaf, aes(x = final,y = 3,
                           label = Letter),
            inherit.aes = F,size = 8,family ="Arial",color = "black") +
  scale_color_manual(values = paleta_final)+
  xlab(NULL)+ 
  ylab('Shannon index')+
  clean +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=20),
        plot.title = element_text(hjust = 0.5, size=20),
        legend.position = 'none')

p_leaf

######save figures
oh.save.pdf(p = p_leaf,
            outname = "figS3H-alpha_diversity_syncom_experiment.pdf",
            outdir = "../figures/",width = 20,height = 10)
