######load packages
library(ohchibi)
library(dplyr)
library(egg)
library(ggrepel)
library(car) 
library(rcompanion)#cldList function
library(FSA)
library(boot)

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
df <-  read.csv("../cleandata/leaf_biomass.csv", header = T, sep = ",",
                quote = "",comment.char = "")

######transform the factor
df$genotype <- df$genotype %>% factor
df$soil <- df$soil %>% factor
df$treatment <- df$treatment %>% factor
df$num_leaf <- df$num_leaf %>% factor
df$tray <- df$tray %>% factor
df$experiment <- df$experiment %>% factor
df$region <- df$region %>% factor
df$country <- df$country %>% factor

###Paletta soil
paleta_soil<- c("#FF7F00","#A6CEE3", '#1F78B4')
names(paleta_soil) <- c( "Sao Domingos", "Tarrafal",
                         "Sutton Bonington")

######paleta leaf
paleta_leaf <- c('#d01c8b','#f1b6da','#7b3294',
                 '#c2a5cf','#a6dba0','#008837')
names(paleta_leaf) <- c('L1', 'L2', 'L3', 'L4', 'L5', 'L6')

######select the 50%FC
df_sub <- df %>%
  subset(treatment == '50%FC') %>% droplevels

df_sub <- df_sub %>%
  subset(num_leaf != 'Stem') %>% droplevels

df_sub <- df_sub %>%
  subset(dw_mg != 'NA') %>% droplevels

######change the na,e of sao domingos
df_sub$soil <- df_sub$soil %>%
  gsub(pattern = 'Sao Domingos', replacement = 'Sao Domingos')

######order the soil names
df_sub$soil <- df_sub$soil %>%
  factor(levels = c('Sao Domingos', 'Tarrafal','Sutton Bonington'))

######create a group column with num_pot and experiment
#df_sub$group <- paste(df_sub$num_pot, df_sub$experiment,sep='_')

######average all the leaves/pot
df_sum <- aggregate(dw_mg ~ soil+id,
                     df_sub, sum)

######statistical assumptions
#1.independence of the variables

#2.variance homogeneity
#perform Bartlett's test
bartlett.test(dw_mg ~ soil, data = df_sum)
#p<0.05. Do not reject the null hypothesis
#3.check the normality of the data
qqnorm(df_sum$dw_mg)
qqline(df_sum$dw_mg, col='red')

##Data exploration
ggplot(data = df_sum,mapping = aes(x = soil,y = dw_mg)) +
  geom_boxplot() + 
  geom_point()

######apply kruskal-wallis test
m1 <- kruskal.test(dw_mg ~ soil,
                   data = df_sum)

m1$p.value

DT = dunnTest(dw_mg ~ soil,data=df_sum,method="bh")

PT = DT$res

PT

df_letter <- cldList(P.adj ~ Comparison,
                     data = PT,
                     threshold = 0.05)

######calculate the mean and the confidence interval for each analysis
df_sb <- df_sum %>% subset(soil == 'Sutton Bonington') %>% droplevels
df_mean_sb <- mean(df_sb$dw_mg) %>% as.data.frame
df_mean_sb$soil <- 'SuttonBonington'
colnames(df_mean_sb)[1] <- 'mean'

df_sd <- df_sum %>% subset(soil == 'Sao Domingos') %>% droplevels
df_mean_sd <- mean(df_sd$dw_mg) %>% as.data.frame
df_mean_sd$soil <- 'SaoDomingos'
colnames(df_mean_sd)[1] <- 'mean'

df_t <- df_sum %>% subset(soil == 'Tarrafal') %>% droplevels
df_mean_t <- mean(df_t$dw_mg) %>% as.data.frame
df_mean_t$soil <- 'Tarrafal'
colnames(df_mean_t)[1] <- 'mean'

######combine all the datasets
Res_mean <- rbind(df_mean_sd, df_mean_t, df_mean_sb)

######calculate the ci
results_sb <- boot(data=df_sb$dw_mg, statistic=Bmean, R=1000)
results_boot_sb <- boot.ci(results_sb, type=c("norm", "basic", "perc", "bca"))
df_ci_sb <- results_boot_sb$bca %>% as.data.frame
df_ci_sb$soil <- 'SuttonBonington'

results_sd <- boot(data=df_sd$dw_mg, statistic=Bmean, R=1000)
results_boot_sd <- boot.ci(results_sd, type=c("norm", "basic", "perc", "bca"))
df_ci_sd <- results_boot_sd$bca %>% as.data.frame
df_ci_sd$soil <- 'SaoDomingos'

results_t <- boot(data=df_t$dw_mg, statistic=Bmean, R=1000)
results_boot_t <- boot.ci(results_t, type=c("norm", "basic", "perc", "bca"))
df_ci_t <- results_boot_t$bca %>% as.data.frame
df_ci_t$soil <- 'Tarrafal'

######combine ci
Res_ci <- rbind(df_ci_sd, df_ci_t, df_ci_sb)

###select column from 4to6
Res_ci <- Res_ci[,4:6]
###rename the first column
colnames(Res_ci)[1]<-'lower'
colnames(Res_ci)[2]<-'upper'

######merge the two datasets
res <- merge(Res_mean, Res_ci, by='soil')

######change the group title
colnames(df_letter)[1] <- 'soil'

######merge res and letter
df <- merge(res, df_letter, by='soil')

######correct the soil names
df$soil <- df$soil %>%
  gsub(pattern = 'SaoDomingos', replacement = 'São Domingos')

df$soil <- df$soil %>%
  gsub(pattern = 'SuttonBonington', replacement = 'Sutton Bonington')

######organize the df soil column
df$soil <- df$soil %>%
  factor(levels = c('São Domingoss', 'Tarrafal', 'Sutton Bonington'))
df_sum$soil <- df_sum$soil %>%
  factor(levels = c('São Domingos', 'Tarrafal', 'Sutton Bonington'))

######plot the result
p_leaf <- ggplot(data=df_sum, aes(x=soil, y=dw_mg, 
                                     shape=soil, color=soil))+
  geom_point(shape=21, size=4, alpha=0.3)+
  geom_pointrange(data=df, aes(y=mean, ymin=lower, ymax=upper, 
                               shape=soil), 
                  size=1.5)+
  geom_text(data = df, aes(x = soil,y = 615,
                           label = Letter),
            inherit.aes = F,size = 10,family ="Arial",color = "black") +
  scale_color_manual(values = paleta_soil)+
  #ggtitle('All leaves')+
  xlab(NULL)+ 
  ylab('Leaf weight (mg)')+
  ylim(0,620)+
  clean +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 30),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 30),
        legend.position = 'right',
        plot.title = element_text(hjust = 0.5, size=30, face = 'bold'))

p_leaf

######save figures
oh.save.pdf(p = p_leaf,outname = "fig1A.1_total_leaf_weight.pdf",
            outdir = "../figures/",width = 20,height = 15)
