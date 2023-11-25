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
df <-  read.csv("../cleandata/leaf_length_stem_pilot_test_final.csv", header = T, sep = ",",
                quote = "",comment.char = "")

######transform the factor
df$treatment <- df$treatment %>% factor
df$num_leaf <- df$num_leaf %>% factor

######create color treatment
paleta_treat <- c('#c7eae5','#80cdc1','#f6e8c3','#dfc27d')
names(paleta_treat) <- c('Nonutrient', 'quarter', 
                         'half', 'Full')

paleta_new <- c('#c7eae5','#80cdc1','#f6e8c3','#dfc27d')
names(paleta_new) <- c('No\nnutrient', '1/4\nHoagland', 
                         '1/2\nHoagland', 'Full\nHoagland')

######paleta leaf
paleta_leaf <- c('#d01c8b','#f1b6da','#7b3294',
                 '#c2a5cf','#a6dba0','#008837')
names(paleta_leaf) <- c('L1', 'L2', 'L3', 'L4', 'L5', 'L6')


######average per treatment
df_sum <- aggregate(length_stem~UiD+treatment,df,sum)

#######find influential points before run ANOVA
###Step01:fit the linear regression model to the dataset with outliers
model <- lm(length_stem ~ treatment, data = df_sum)
###Step02:find Cook's distance for each observation in the dataset
cooksD <- cooks.distance(model)
###Step03:Plot Cook's Distance with a horizontal line at 4/n to 
#see which observations exceed this thresdhold
n <- nrow(df_sum)
plot(cooksD, main = "Cooks Distance for Influential Obs")
abline(h = 4/n, lty = 2, col = "steelblue") # add cutoff line
###Step04:identify influential points
influential_obs <- as.numeric(names(cooksD)[(cooksD > (4/n))])
#define new data frame with influential points removed
df_removed <- df[-influential_obs, ]

######statistical assumptions
#1.independence of the variables

#2.variance homogeneity
#perform levene test
leveneTest(length_stem ~ treatment, data = df_sum)
#p>0.05. Do not Reject H0.
#Not similar variance
#3.check the normality of the data
qqnorm(df_sum$length_stem)
qqline(df_sum$length_stem, col='red')

##Data exploration
ggplot(data = df_sum,mapping = aes(x = treatment,y = length_stem)) +
  geom_boxplot() + 
  geom_point()

######organize the levels of df_sum
df_sum$treatment <- df_sum$treatment %>%
  factor(levels = c('Nonutrient','quarter','half','Full'))

#######Fit the ANOVA model (one-way model)
m1 <- aov(data = df_sum, 
          formula = length_stem ~ treatment)
summary(m1)
#have a significant group effect

######Elaborate the letter data frame
msum <- glht(m1, linfct = mcp(treatment='Tukey')) %>% summary
df_treat <- emmeans(m1,specs = "treatment")  %>% as.data.frame
cld_treat <- cld(msum)

cld_treat <- cld_treat$mcletters$Letters %>% as.data.frame
cld_treat$treatment <- rownames(cld_treat)
rownames(cld_treat) <- NULL
colnames(cld_treat)[1] <- 'letters'

df_treat <- merge(df_treat, cld_treat, by='treatment')

######edit the treatment name
df_sum$new <- ''
df_sum$new[which(df_sum$treatment=='Nonutrient')] <- 'No\nnutrient'
df_sum$new[which(df_sum$treatment=='quarter')] <- '1/4\nHoagland'
df_sum$new[which(df_sum$treatment=='half')] <- '1/2\nHoagland'
df_sum$new[which(df_sum$treatment=='Full')] <- 'Full\nHoagland'

df_treat$new <- ''
df_treat$new[which(df_treat$treatment=='Nonutrient')] <- 'No\nnutrient'
df_treat$new[which(df_treat$treatment=='quarter')] <- '1/4\nHoagland'
df_treat$new[which(df_treat$treatment=='half')] <- '1/2\nHoagland'
df_treat$new[which(df_treat$treatment=='Full')] <- 'Full\nHoagland'

######order the new position
df_sum$new <- df_sum$new %>% factor(levels = c('No\nnutrient', '1/4\nHoagland',
                                               '1/2\nHoagland', 'Full\nHoagland'))

######plot the result
p_leaf <- ggplot(data=df_sum, aes(x=new, y=length_stem,color=new))+
  geom_point(shape=21, size=4, alpha=5)+
  geom_pointrange(data=df_treat, aes(y=emmean, ymin=lower.CL, ymax=upper.CL), 
                  size=1.5)+
  geom_text(data = df_treat, aes(x = new,y = 200,
                           label = letters),
            inherit.aes = F,size = 6,family ="Arial",color = "black") +
  scale_color_manual(values = paleta_new)+
  xlab(NULL)+ 
  ylab('Leaf length (cm)')+
  ylim(0,200)+
  clean +
  theme(strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 30),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size=30, face = 'bold'))

p_leaf

######save the figures
oh.save.pdf(p = p_leaf,
            outname = "figS3C_leaf_lenght_pilot_final_ANOVA.pdf",
            outdir = "../figures/",width = 10,height = 10)
