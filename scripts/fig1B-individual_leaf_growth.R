######load library
library(ohchibi)
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
# function to obtain the mean
Bmean <- function(data, indices) {
  d <- data[indices] # allows boot to select sample 
  return(mean(d))
}

######upload data
df <-  read.csv("../cleandata/leaf_metadata_length.csv", header = T, sep = ",",
                quote = "",comment.char = "")

######correct the name in the first column
colnames(df)[1] <- 'UiD'

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
names(paleta_soil) <- c("Sao Domingos", "Tarrafal",
                        "Sutton Bonington")

######paleta leaf
paleta_leaf <- c('#d01c8b','#f1b6da','#7b3294',
                 '#c2a5cf','#a6dba0','#008837')
names(paleta_leaf) <- c('L1', 'L2', 'L3', 'L4', 'L5', 'L6')

######select the 50%FC
df_50 <- df %>% subset((treatment == '50%FC') &
                         (num_leaf != 'Stem')) %>% droplevels

df_50 <- df_50 %>% subset(leaf_lenght_cm != 'NA') %>% droplevels

#####order the soils levels 
df_50$soil <- df_50$soil %>%
  gsub(pattern = 'Sao Domingos', replacement = 'Sao Domingos')

df_50$soil <- df_50$soil %>%
  factor(levels = c('Sao Domingos', 'Tarrafal', 'Sutton Bonington'))

######create a group with NumLeaf and Soil
df_50$group <- paste(df_50$soil, df_50$num_leaf, sep='_')

######statistical assumptions
#1.independence of the variables

#2.variance homogeneity
#perform levene test
leveneTest(leaf_lenght_cm ~ group, data = df_50)
#p-value < 0.05. Do not reject the null hypoyhesis that the group have similar
#3.check the normality of the data
qqnorm(df_50$leaf_lenght_cm)
qqline(df_50$leaf_lenght_cm, col='red')

##Data exploration
ggplot(data = df_50,mapping = aes(x = soil,y =leaf_lenght_cm)) +
  geom_boxplot() + 
  geom_point()+
  facet_grid(.~num_leaf)

######define group as a factor
df_50$group <- df_50$group %>% factor

######perform kruskal walis for each leaves
l1 <- df_50 %>% subset(num_leaf == 'L1') %>% droplevels

######apply kruskal-wallis test
m1_l1 <- kruskal.test(leaf_lenght_cm ~ group,
                      data = l1)

DT_l1 = dunnTest(leaf_lenght_cm ~ group,data=l1,method="bh")

PT_l1 = DT_l1$res

PT_l1

df_letter_l1 <- cldList(P.adj ~ Comparison,
                        data = PT_l1,
                        threshold = 0.05)
###########subset leaf 2
l2 <- df_50 %>% subset(num_leaf == 'L2') %>% droplevels

######apply kruskal-wallis test
m1_l2 <- kruskal.test(leaf_lenght_cm ~ group,
                      data = l2)

DT_l2 = dunnTest(leaf_lenght_cm ~ group,data=l2,method="bh")

PT_l2 = DT_l2$res

PT_l2

df_letter_l2 <- cldList(P.adj ~ Comparison,
                        data = PT_l2,
                        threshold = 0.05)

###########subset leaf 3
l3 <- df_50 %>% subset(num_leaf == 'L3') %>% droplevels

######apply kruskal-wallis test
m1_l3 <- kruskal.test(leaf_lenght_cm ~ group,
                      data = l3)

DT_l3 = dunnTest(leaf_lenght_cm ~ group,data=l3,method="bh")

PT_l3 = DT_l3$res

PT_l3

df_letter_l3 <- cldList(P.adj ~ Comparison,
                        data = PT_l3,
                        threshold = 0.05)

###########subset leaf 4
l4 <- df_50 %>% subset(num_leaf == 'L4') %>% droplevels

######apply kruskal-wallis test
m1_l4 <- kruskal.test(leaf_lenght_cm ~ group,
                      data = l4)

DT_l4 = dunnTest(leaf_lenght_cm ~ group,data=l4,method="bh")

PT_l4 = DT_l4$res

PT_l4

df_letter_l4 <- cldList(P.adj ~ Comparison,
                        data = PT_l4,
                        threshold = 0.05)

###########subset leaf 5
l5 <- df_50 %>% subset(num_leaf == 'L5') %>% droplevels

######apply kruskal-wallis test
m1_l5 <- kruskal.test(leaf_lenght_cm ~ group,
                      data = l5)

DT_l5 = dunnTest(leaf_lenght_cm ~ group,data=l5,method="bh")

PT_l5 = DT_l5$res

PT_l5

df_letter_l5 <- cldList(P.adj ~ Comparison,
                        data = PT_l5,
                        threshold = 0.05)

###########subset leaf 6
l6 <- df_50 %>% subset(num_leaf == 'L6') %>% droplevels

######apply kruskal-wallis test
m1_l6 <- kruskal.test(leaf_lenght_cm ~ group,
                      data = l6)

DT_l6 = dunnTest(leaf_lenght_cm ~ group,data=l6,method="bh")

PT_l6 = DT_l6$res

PT_l6

df_letter_l6 <- cldList(P.adj ~ Comparison,
                        data = PT_l6,
                        threshold = 0.05)

######combine the data sets
df_letter <- rbind(df_letter_l1, df_letter_l2, df_letter_l3,
                   df_letter_l4, df_letter_l5, df_letter_l6)


######calculate the mean and the confidence interval for each analysis
###create the unique names for the for loop
mleaf <- df_50$group %>% unique
###Empty dataset
Res_mean <- NULL
Res_ci <- NULL

###perform the for loop
for (leaf in mleaf){
  df_leaf <- df_50 %>% subset(group == leaf) %>% droplevels
  #calculate means
  df_mean <- mean(df_leaf$leaf_lenght_cm) %>% as.data.frame
  ###calculate ci
  # bootstrapping with 1000 replications 
  results <- boot(data=df_leaf$leaf_lenght_cm, statistic=Bmean, R=1000)
  # get 95% confidence interval 
  results_boot <- boot.ci(results, type=c("norm", "basic", "perc", "bca"))
  ###create a data frame for ci
  df_ci <- results_boot$bca %>% as.data.frame
  #add the leaf name on data frame
  df_mean$group <- leaf
  df_ci$group <- leaf
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
colnames(df_letter)[1] <- 'group'
df_letter$group <- df_letter$group %>% gsub(pattern = ' ', replacement = '')
######prepare a new column with the location information
df_letter$soil <- sapply(strsplit(as.character(df_letter$group),
                                  '_'),`[`, 1)
df_letter$num_leaf <- sapply(strsplit(as.character(df_letter$group),
                                     '_'),`[`, 2)

######reorder the df_group
df_letter$soil <- df_letter$soil %>%
  factor(levels = c('SaoDomingos','Tarrafal','SuttonBonington'))

######merge res and letter
df <- merge(res, df_letter, by='group')

######correct the soil names
df$soil <- df$soil %>%
  gsub(pattern = 'SaoDomingos', replacement = 'Sao Domingos')

df$soil <- df$soil %>%
  gsub(pattern = 'SuttonBonington', replacement = 'Sutton Bonington')

######NumLeaf as factor
df$num_leaf <- df$num_leaf %>% factor

######plot the result
p_leaf <- ggplot(data=df_50, aes(x=soil, y=leaf_lenght_cm, shape=soil,
                                 color=num_leaf))+
  geom_point(shape=21, size=1.5, alpha=1)+
  facet_grid(.~num_leaf)+
  geom_pointrange(data=df, aes(y=mean, ymin=lower, ymax=upper, 
                               color=num_leaf), 
                  size=1)+
  geom_text(data = df, aes(x = soil,y = 80,
                           label = Letter),
            inherit.aes = F,size = 6,family ="Arial",color = "black") +
  scale_color_manual(values = paleta_leaf)+
  xlab(NULL)+ 
  scale_x_discrete(breaks=c("S?o Domingos","Tarrafal","Sutton Bonington"),
                   labels=c("S?o\nDomingos",
                            "Tarrafal", "Sutton\nBonington"))+
  ylab('Leaf length (cm)')+
  clean +
  theme(axis.text.x = element_text(size=12,vjust=1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.position = 'none')

p_leaf

######save figures
oh.save.pdf(p = p_leaf,outname = "fig1B_individual_leaf_growth.pdf",
            outdir = "../figures/",width = 20,height = 10)
