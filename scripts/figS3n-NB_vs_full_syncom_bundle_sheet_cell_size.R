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
# function to obtain the mean
Bmean <- function(data, indices) {
  d <- data[indices] # allows boot to select sample 
  return(mean(d))
}
source('0-Clean_up_plots.R')

######upload data
df <-  read.csv("../rawdata/cell_size_final_final_final.csv", header = T, sep = ",",
                quote = "",comment.char = "")

######transform the factor
df$treatment <- df$treatment %>% factor
df$marker <- df$marker %>% factor
df$nutrient <- df$nutrient %>% factor
df$NumLeaf <- df$NumLeaf %>% factor

######coorect the names
df$NumLeaf <- df$NumLeaf %>% gsub(pattern = ' ', replacement = '')

######select only low hoagland
df_low <- df %>% subset(nutrient == '1/4 Hoagland') %>% droplevels

######average per treatment
df_mean <- aggregate(cell_size_um~UiD+sample_id+treatment+marker+
                        box+NumLeaf,
                     df_low,mean)

######create the dataframe
df_nb <- df_mean %>% subset(treatment == 'No Bacteria')
df_fs <- df_mean %>% subset(treatment == 'Full syncom')

######re-start the index
rownames(df_nb) <- NULL
rownames(df_fs) <- NULL

######change the name for the last column
colnames(df_nb)[7] <- 'nb'
colnames(df_fs)[7] <- 'fs'

######count the number of samples per L1
#Count the number of interview per municipalities
df_nb$NumLeaf <- df_nb$NumLeaf %>% factor
df_fs$NumLeaf <- df_fs$NumLeaf %>% factor

as.data.frame(df_nb %>% xtabs(formula = ~ NumLeaf))
as.data.frame(df_fs %>% xtabs(formula = ~ NumLeaf))

######average per treatment
df_mean_nb <- aggregate(nb~sample_id+treatment+NumLeaf,
                     df_nb,mean)
df_mean_fs <- aggregate(fs~sample_id+treatment+NumLeaf,
                        df_fs,mean)

######remove the outliers from the df_mean datasets
df_mean_nb_t <- aggregate(nb~treatment+NumLeaf,
                          df_nb,mean)
df_mean_fs_t <- aggregate(fs~treatment+NumLeaf,
                          df_fs,mean)

###error
df_error_nb_t <- aggregate(nb~treatment+NumLeaf,
                          df_nb,sd)
df_error_fs_t <- aggregate(fs~treatment+NumLeaf,
                          df_fs,sd)

######merge the mean and the error
df_mean_final_nb <- merge(df_mean_nb_t, df_error_nb_t, by=c('treatment',
                                                         'NumLeaf'))
df_mean_final_fs <- merge(df_mean_fs_t, df_error_fs_t, by=c('treatment',
                                                            'NumLeaf'))

######calculate the 1.5x
df_mean_final_nb$lower <- df_mean_final_nb$nb.x-(1.5*df_mean_final_nb$nb.y)
df_mean_final_nb$upper <- df_mean_final_nb$nb.x+(1.5*df_mean_final_nb$nb.y)

df_mean_final_fs$lower <- df_mean_final_fs$fs.x-(1.5*df_mean_final_fs$fs.y)
df_mean_final_fs$upper <- df_mean_final_fs$fs.x+(1.5*df_mean_final_fs$fs.y)

######remove the plants with lower or higher than sd
df_nb <- df_nb[-9,]

######count the number of 
as.data.frame(df_nb %>% xtabs(formula = ~ NumLeaf))
as.data.frame(df_fs %>% xtabs(formula = ~ NumLeaf))

######re-start the index
rownames(df_nb) <- NULL
rownames(df_fs) <- NULL

######remove one samples for L2 nb and 1 samples for L4 and L5 for fs
#######compare NB vs fs
df_nb <- df_nb[-21,]
df_fs_nb <- df_fs[-c(1,8),]

as.data.frame(df_nb %>% xtabs(formula = ~ NumLeaf))
as.data.frame(df_fs_nb %>% xtabs(formula = ~ NumLeaf))

#######select only the NumLeaf
df_final_nbfs <- cbind(df_nb, df_fs_nb)

######select the column
df_final_nbfs <- df_final_nbfs[,c(6,7,14)]

######calculate the ratio
df_final_nbfs$ratio <- df_final_nbfs$nb/df_final_nbfs$fs
df_final_nbfs$contrast <- 'No Bacteria vs Full SynCom'
df_final_nbfs <- df_final_nbfs[,c(1,4,5)]

######combine the datasets
df_final <- rbind(df_final_nbfs)

######statistical assumptions
#1.independence of the variables

#2.variance homogeneity
#perform levene test
leveneTest(ratio ~ NumLeaf, data = df_final)
#H0:variance among each group are equal.
#p>0.05.Do not reject H0. group have similar variance
#3.check the normality of the data
qqnorm(df_final$ratio)
qqline(df_final$ratio, col='red')

##Data exploration
ggplot(data = df_final,mapping = aes(x = NumLeaf, y=ratio)) +
  geom_boxplot() + 
  geom_point()+
  ylim(0,2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_grid(.~contrast)

######create a group with NumLeaf and contrast
df_final$group <- paste0(df_final$NumLeaf,'_', df_final$contrast)

######transform group variable as factor
df_final$group <- df_final$group %>% factor

#######fit the ANOVA by contrast
df_fsnb <- df_final %>% subset(contrast == 'No Bacteria vs Full SynCom')

#######Fit the ANOVA model (one-way model)
m1_fsnb <- aov(data = df_fsnb, 
          formula = ratio ~ group)
summary(m1_fsnb)

######Elaborate the letter data frame
msum_fsnb <- glht(m1_fsnb, linfct = mcp(group='Tukey')) %>% summary
df_fsnb <- emmeans(m1_fsnb,specs = "group")  %>% as.data.frame
cld_fsnb <- cld(msum_fsnb)

cld_fsnb <- cld_fsnb$mcletters$Letters %>% as.data.frame
cld_fsnb$group <- rownames(cld_fsnb)
rownames(cld_fsnb) <- NULL
colnames(cld_fsnb)[1] <- 'letters'

####merge the datasets
df_fsnb <- merge(df_fsnb, cld_fsnb, by='group')

######merge the letters
df_letters <- rbind(df_fsnb)

######create the numLeaf and the contrast column
df_letters$NumLeaf <- sapply(strsplit(as.character(df_letters$group),
                                  '_'),`[`, 1)
df_letters$contrast <- sapply(strsplit(as.character(df_letters$group),
                                      '_'),`[`, 2)

######paleta leaf
paleta_leaf <- c('#d01c8b','#f1b6da','#7b3294',
                 '#c2a5cf','#a6dba0','#008837')
names(paleta_leaf) <- c('L1', 'L2', 'L3', 'L4', 'L5', 'L6')

######order the contrast
df_final$contrast <- df_final$contrast %>%
  factor(levels = c("No Bacteria vs Full SynCom"))

df_letters$contrast <- df_letters$contrast %>%
  factor(levels = c("No Bacteria vs Full SynCom"))

######select only the No Bacteria vs Full SynCom
df_final_final <- df_final %>%
  subset(contrast == 'No Bacteria vs Full SynCom')

df_letters_nb <- df_letters %>%
  subset(contrast == 'No Bacteria vs Full SynCom')

######plot the result
p_leaf <- ggplot(data=df_final_final, aes(x=NumLeaf, y=ratio, 
                                 color=NumLeaf))+
  geom_point(shape=21, size=4, alpha=1.5)+
  facet_grid(.~contrast)+
  geom_pointrange(data=df_letters_nb, aes(y = emmean, ymin = lower.CL, 
                                       ymax = upper.CL, 
                               color=NumLeaf), size=1.5)+
  geom_text(data = df_letters_nb, aes(x = NumLeaf,y = 2,
                           label = letters),
            inherit.aes = F,size = 6,family ="Arial",color = "black") +
  scale_color_manual(values = paleta_leaf)+
  xlab(NULL)+ 
  ylab('ratio cell size')+
  ylim(0.5,2)+
  clean +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size=5, face = 'bold'))

p_leaf

######save the figures
oh.save.pdf(p = p_leaf,
            outname = "figS3n_ratio_per_leaf_dropout_nbvssyncom.pdf",
            outdir = "../figures/", width = 5,height = 12)

ggsave(plot = p_leaf,'../figure/figS3K_10_ratio_per_leaf_dropout_compared_NB.png',
       units='cm',
       width=40, height=20)

