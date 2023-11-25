######load library
library(ohchibi)
library(tidyverse)
library(emmeans)
library(multcomp)
library(car)
library(FSA)
library(rcompanion)
library(boot)
library(rstatix)

######set seed
set.seed(12356)

######set work directory
setwd("C:/Users/Castrillo Lab/OneDrive - The University of Nottingham/Publication 2019/Paper preparation/manuscripts/scripts")

######helper functions
source('0-Clean_up_plots.R')
Bmean <- function(data, indices) {
  d <- data[indices] # allows boot to select sample 
  return(mean(d))
}

######open file
df_final <- read.csv(file = '../cleandata/df_removed_final_stem_repeat.csv')

######change the treatment to +SynCom
df_final$treatment[which(df_final$treatment == 'syncom')] <- '+SynCom'

#######remove outliers
df <- df_final %>%
  subset(pot != 53)
df <- df %>%
  subset(pot != 65)
df <- df %>%
  subset(pot != 110)
df <- df %>%
  subset(pot != 111)
df <- df %>%
  subset(pot != 120)
df <- df %>%
  subset(pot != 128)

######remove the line 18
df_1 <- df %>%
  subset(stock != 'UFMu-00214')

df <- df_1

######define treatment as factor
df$treatment <- df$treatment %>%
  factor(levels = c('NB', '+SynCom'))

df$gene <- df$gene %>%
  factor(levels = c('W22', 'zm00001eb022810',
                    'zm00001eb290350',
                    'zm00001eb149250_#1', 'zm00001eb149250_#2',
                    'zm00001eb149250_#3'))

######plot
ggplot(data = df,mapping = aes(x = gene,y =NormValue, 
                                     color=treatment)) +
  geom_boxplot() + 
  geom_point()+
  facet_grid(experiment~leaf) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

######remove the three alelles not selected
#df_final_1 <- df_final %>%
  #subset(gene != 'Zm00001eb149250_#1') %>% droplevels
#df_final_1 <- df_final_1 %>%
  #subset(gene != 'Zm00001eb149250_#2') %>% droplevels
#df_final_1 <- df_final_1 %>%
  #subset(gene != 'Zm00001eb149250_#3') %>% droplevels

######select L1
df_l1_syncom <- df %>%
  subset((leaf == 'L1') & (treatment == '+SynCom')) %>% droplevels

######apply Kruskal wallis
m1_l1_syncom <- kruskal.test(NormValue ~ gene,
                             data = df_l1_syncom)

m1_l1_syncom$p.value

DT_l1_syncom = dunnTest(NormValue ~ gene,
                        data = df_l1_syncom,method="bh")

PT_l1_syncom = DT_l1_syncom$res

# Use subset to select rows containing 'W22' in the Comparison column
selected_l1_syncom <- subset(PT_l1_syncom, grepl("W22", Comparison))

df_letter_l1_syncom <- cldList(P.adj ~ Comparison,
                               data = selected_l1_syncom,
                               remove.zero = F,
                               threshold = 0.1)

######add column 
df_letter_l1_syncom$comparison_syncom <- '+Syncom W22 vs +SynCom mutants'

######add a treatment column
df_letter_l1_syncom$treatment <- '+SynCom'

######add the leaf column
df_letter_l1_syncom$leaf <- 'L1'

######rename the first column
colnames(df_letter_l1_syncom)[1] <- 'gene'

######select L1
df_l1 <- df %>%
  subset(leaf == 'L1') %>% droplevels

######select the files
mtreat <- df_l1$treatment %>% as.character %>% unique
mgene <- df_l1$gene %>% as.character %>% unique

######create an empty 
Res_letter_l1 <- NULL

######perform the for loop for each leaf
for (i in mgene) {
  df_gene <- df_l1 %>% subset(gene == i) %>% droplevels
  m1 <- kruskal.test(NormValue ~ treatment, data = df_gene)
  DT = wilcox.test(NormValue ~ treatment, data = df_gene,method="bh")
  df_letter <- c(p.value = DT$p.value) %>% data.frame
  df_letter$leaf <- 'L1'
  df_letter$gene <- i
  Res_letter_l1 <- rbind(df_letter, Res_letter_l1)
}

######rename column 1
colnames(Res_letter_l1)[1] <- 'p.value'

######empty data frame
Res_mean <- NULL
Res_ci <- NULL

###perform the for loop
for (i in mgene) {
  df_gene <- df_l1 %>% subset(gene == i) %>% droplevels
  for (treat in mtreat) {
    df_treat <- df_gene %>% subset(treatment == treat) %>% droplevels
    #calculate means
    df_mean_l1 <- mean(df_treat$NormValue) %>% as.data.frame
    ###calculate ci
    # bootstrapping with 1000 replications 
    results_l1 <- boot(data=df_treat$NormValue, statistic=Bmean, R=1000)
    # get 95% confidence interval 
    results_boot_l1 <- boot.ci(results_l1, conf = 0.90, 
                               type=c("norm", "basic", "perc", "bca"))
    ###create a data frame for ci
    df_ci_l1 <- results_boot_l1$bca %>% as.data.frame
    df_mean_l1$gene <- i
    df_mean_l1$treatment <- treat
    df_ci_l1$gene <- i
    df_ci_l1$treatment <- treat
    #combine the results
    Res_mean <- rbind(Res_mean, df_mean_l1)
    Res_ci <- rbind(Res_ci, df_ci_l1)
  }
}

######combine the two datasets
colnames(Res_mean)[1]<-'mean'
###select column from 4to7
Res_ci <- Res_ci[,4:7]
###rename the first column
colnames(Res_ci)[1]<-'lower'
colnames(Res_ci)[2]<-'upper'

######combine the two datasets
Res_em <- merge(Res_mean, Res_ci, by=c('gene', 'treatment'))

#####create color for the dataset
paleta_syncom <- c('#80cdc1','#018571')
names(paleta_syncom) <- c('NB','+SynCom')

######order the treatment
df_l1$treatment <- df_l1$treatment %>%
  factor(levels = c('NB', '+SynCom'))

Res_em$treatment <- Res_em$treatment %>%
  factor(levels = c('NB', '+SynCom'))

######combine gene and treatment
df_l1$group1 <- paste(df_l1$gene, df_l1$treatment, sep='_')
Res_em$group1 <- paste(Res_em$gene, Res_em$treatment, sep='_')
df_letter_l1_syncom$group1 <- paste(df_letter_l1_syncom$gene,
                                    df_letter_l1_syncom$treatment, sep='_')

######order the factors
df_l1$group1 <- df_l1$group1 %>%
  factor(levels = c('W22_NB', 'W22_+SynCom',
                    'zm00001eb022810_NB', 'zm00001eb022810_+SynCom',
                    'zm00001eb290350_NB','zm00001eb290350_+SynCom',
                    'zm00001eb149250_#1_NB','zm00001eb149250_#1_+SynCom',
                    'zm00001eb149250_#2_NB','zm00001eb149250_#2_+SynCom',
                    'zm00001eb149250_#3_NB','zm00001eb149250_#3_+SynCom'))

Res_em$group1 <- Res_em$group1 %>%
  factor(levels = c('W22_NB', 'W22_+SynCom',
                    'zm00001eb022810_NB', 'zm00001eb022810_+SynCom',
                    'zm00001eb290350_NB','zm00001eb290350_+SynCom',
                    'zm00001eb149250_#1_NB','zm00001eb149250_#1_+SynCom',
                    'zm00001eb149250_#2_NB','zm00001eb149250_#2_+SynCom',
                    'zm00001eb149250_#3_NB','zm00001eb149250_#3_+SynCom'))

df_letter_l1_syncom$group1 <- df_letter_l1_syncom$group1 %>%
  factor(levels = c('W22_NB', 'W22_+SynCom',
                    'zm00001eb022810_NB', 'zm00001eb022810_+SynCom',
                    'zm00001eb290350_NB','zm00001eb290350_+SynCom',
                    'zm00001eb149250_#1_NB','zm00001eb149250_#1_+SynCom',
                    'zm00001eb149250_#2_NB','zm00001eb149250_#2_+SynCom',
                    'zm00001eb149250_#3_NB','zm00001eb149250_#3_+SynCom'))

######plot the results
p <- ggplot(data=df_l1, aes(x=group1, y=NormValue,
                            color=treatment))+
  geom_point(aes(group=treatment),shape=21, size=1.5, alpha=1)+
  geom_pointrange(data=Res_em, aes(y=mean,
                                   ymin=lower, ymax=upper, 
                                   color= treatment),
                  size=1)+
  geom_text(data = df_letter_l1_syncom, aes(x = group1,y = 70,
                                              label = Letter),
            inherit.aes = F,size = 8,family ="Arial",color = "black") +
  ylab('Leaf length (cm)')+
  ggtitle('L1')+
  ylim(0,70)+
  scale_color_manual(values = paleta_syncom)+
  clean +
  theme(axis.text.x = element_text(size=20,vjust=1, angle = 45, hjust=1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, face='bold'),
        axis.title.x = element_blank(),
        plot.title = element_text(size=20, hjust=0.5),
        legend.position = 'none')

p

######select L2
df_l2_syncom <- df %>%
  subset((leaf == 'L2') & (treatment == '+SynCom')) %>% droplevels

######apply Kruskal wallis
m1_l2_syncom <- kruskal.test(NormValue ~ gene,
                             data = df_l2_syncom)

m1_l2_syncom$p.value

DT_l2_syncom = dunnTest(NormValue ~ gene,
                        data = df_l2_syncom,method="bh")

PT_l2_syncom = DT_l2_syncom$res

# Use subset to select rows containing 'W22' in the Comparison column
selected_l2_syncom <- subset(PT_l2_syncom, grepl("W22", Comparison))

df_letter_l2_syncom <- cldList(P.adj ~ Comparison,
                               data = selected_l2_syncom,
                               remove.zero = F,
                               threshold = 0.1)

######add column 
df_letter_l2_syncom$comparison_syncom <- '+Syncom W22 vs +SynCom mutants'

######add a treatment column
df_letter_l2_syncom$treatment <- '+SynCom'

######add the leaf column
df_letter_l2_syncom$leaf <- 'L2'

######rename the first column
colnames(df_letter_l2_syncom)[1] <- 'gene'

######select l2
df_l2 <- df %>%
  subset(leaf == 'L2') %>% droplevels

######select the files
mtreat <- df_l2$treatment %>% as.character %>% unique
mgene <- df_l2$gene %>% as.character %>% unique

######create an empty 
Res_letter_l2 <- NULL

######perform the for loop for each leaf
for (i in mgene) {
  df_gene <- df_l2 %>% subset(gene == i) %>% droplevels
  m1 <- kruskal.test(NormValue ~ treatment, data = df_gene)
  DT = wilcox.test(NormValue ~ treatment, data = df_gene,method="bh")
  df_letter <- c(p.value = DT$p.value) %>% data.frame
  df_letter$leaf <- 'L2'
  df_letter$gene <- i
  Res_letter_l2 <- rbind(df_letter, Res_letter_l2)
}

######rename column 1
colnames(Res_letter_l2)[1] <- 'p.value'

######empty data frame
Res_mean <- NULL
Res_ci <- NULL

###perform the for loop
for (i in mgene) {
  df_gene <- df_l2 %>% subset(gene == i) %>% droplevels
  for (treat in mtreat) {
    df_treat <- df_gene %>% subset(treatment == treat) %>% droplevels
    #calculate means
    df_mean_l2 <- mean(df_treat$NormValue) %>% as.data.frame
    ###calculate ci
    # bootstrapping with 1000 replications 
    results_l2 <- boot(data=df_treat$NormValue, statistic=Bmean, R=1000)
    # get 95% confidence interval 
    results_boot_l2 <- boot.ci(results_l2, conf = 0.90, 
                               type=c("norm", "basic", "perc", "bca"))
    ###create a data frame for ci
    df_ci_l2 <- results_boot_l2$bca %>% as.data.frame
    df_mean_l2$gene <- i
    df_mean_l2$treatment <- treat
    df_ci_l2$gene <- i
    df_ci_l2$treatment <- treat
    #combine the results
    Res_mean <- rbind(Res_mean, df_mean_l2)
    Res_ci <- rbind(Res_ci, df_ci_l2)
  }
}

######combine the two datasets
colnames(Res_mean)[1]<-'mean'
###select column from 4to7
Res_ci <- Res_ci[,4:7]
###rename the first column
colnames(Res_ci)[1]<-'lower'
colnames(Res_ci)[2]<-'upper'

######combine the two datasets
Res_em <- merge(Res_mean, Res_ci, by=c('gene', 'treatment'))

#####create color for the dataset
paleta_syncom <- c('#80cdc1','#018571')
names(paleta_syncom) <- c('NB','+SynCom')

######order the treatment
df_l2$treatment <- df_l2$treatment %>%
  factor(levels = c('NB', '+SynCom'))

Res_em$treatment <- Res_em$treatment %>%
  factor(levels = c('NB', '+SynCom'))

######combine gene and treatment
df_l2$group1 <- paste(df_l2$gene, df_l2$treatment, sep='_')
Res_em$group1 <- paste(Res_em$gene, Res_em$treatment, sep='_')
df_letter_l2_syncom$group1 <- paste(df_letter_l2_syncom$gene,
                                    df_letter_l2_syncom$treatment, sep='_')

######order the factors
df_l2$group1 <- df_l2$group1 %>%
  factor(levels = c('W22_NB', 'W22_+SynCom',
                    'zm00001eb022810_NB', 'zm00001eb022810_+SynCom',
                    'zm00001eb290350_NB','zm00001eb290350_+SynCom',
                    'zm00001eb149250_#1_NB','zm00001eb149250_#1_+SynCom',
                    'zm00001eb149250_#2_NB','zm00001eb149250_#2_+SynCom',
                    'zm00001eb149250_#3_NB','zm00001eb149250_#3_+SynCom'))

Res_em$group1 <- Res_em$group1 %>%
  factor(levels = c('W22_NB', 'W22_+SynCom',
                    'zm00001eb022810_NB', 'zm00001eb022810_+SynCom',
                    'zm00001eb290350_NB','zm00001eb290350_+SynCom',
                    'zm00001eb149250_#1_NB','zm00001eb149250_#1_+SynCom',
                    'zm00001eb149250_#2_NB','zm00001eb149250_#2_+SynCom',
                    'zm00001eb149250_#3_NB','zm00001eb149250_#3_+SynCom'))

df_letter_l2_syncom$group1 <- df_letter_l2_syncom$group1 %>%
  factor(levels = c('W22_NB', 'W22_+SynCom',
                    'zm00001eb022810_NB', 'zm00001eb022810_+SynCom',
                    'zm00001eb290350_NB','zm00001eb290350_+SynCom',
                    'zm00001eb149250_#1_NB','zm00001eb149250_#1_+SynCom',
                    'zm00001eb149250_#2_NB','zm00001eb149250_#2_+SynCom',
                    'zm00001eb149250_#3_NB','zm00001eb149250_#3_+SynCom'))

######plot the results
p_l2 <- ggplot(data=df_l2, aes(x=group1, y=NormValue,
                               color=treatment))+
  geom_point(aes(group=treatment),shape=21, size=1.5, alpha=1)+
  geom_pointrange(data=Res_em, aes(y=mean,
                                   ymin=lower, ymax=upper, 
                                   color= treatment),
                  size=1)+
  geom_text(data = df_letter_l2_syncom, aes(x = group1,y = 70,
                                            label = Letter),
            inherit.aes = F,size = 8,family ="Arial",color = "black") +
  ylab('Leaf length (cm)')+
  ggtitle('L2')+
  ylim(0,70)+
  scale_color_manual(values = paleta_syncom)+
  clean +
  theme(axis.text.x = element_text(size=20,vjust=1, angle = 45, hjust=1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, face='bold'),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size=20, hjust=0.5),
        legend.position = 'none')

p_l2

######select l3
df_l3_syncom <- df %>%
  subset((leaf == 'L3') & (treatment == '+SynCom')) %>% droplevels

######apply Kruskal wallis
m1_l3_syncom <- kruskal.test(NormValue ~ gene,
                             data = df_l3_syncom)

m1_l3_syncom$p.value

DT_l3_syncom = dunnTest(NormValue ~ gene,
                        data = df_l3_syncom,method="bh")

PT_l3_syncom = DT_l3_syncom$res

# Use subset to select rows containing 'W22' in the Comparison column
selected_l3_syncom <- subset(PT_l3_syncom, grepl("W22", Comparison))

df_letter_l3_syncom <- cldList(P.adj ~ Comparison,
                               data = selected_l3_syncom,
                               remove.zero = F,
                               threshold = 0.1)

######add column 
df_letter_l3_syncom$comparison_syncom <- '+Syncom W22 vs +SynCom mutants'

######add a treatment column
df_letter_l3_syncom$treatment <- '+SynCom'

######add the leaf column
df_letter_l3_syncom$leaf <- 'L3'

######rename the first column
colnames(df_letter_l3_syncom)[1] <- 'gene'

######select l3
df_l3 <- df %>%
  subset(leaf == 'L3') %>% droplevels

######select the files
mtreat <- df_l3$treatment %>% as.character %>% unique
mgene <- df_l3$gene %>% as.character %>% unique

######create an empty 
Res_letter_l3 <- NULL

######perform the for loop for each leaf
for (i in mgene) {
  df_gene <- df_l3 %>% subset(gene == i) %>% droplevels
  m1 <- kruskal.test(NormValue ~ treatment, data = df_gene)
  DT = wilcox.test(NormValue ~ treatment, data = df_gene,method="bh")
  df_letter <- c(p.value = DT$p.value) %>% data.frame
  df_letter$leaf <- 'L3'
  df_letter$gene <- i
  Res_letter_l3 <- rbind(df_letter, Res_letter_l3)
}

######rename column 1
colnames(Res_letter_l3)[1] <- 'p.value'

######empty data frame
Res_mean <- NULL
Res_ci <- NULL

###perform the for loop
for (i in mgene) {
  df_gene <- df_l3 %>% subset(gene == i) %>% droplevels
  for (treat in mtreat) {
    df_treat <- df_gene %>% subset(treatment == treat) %>% droplevels
    #calculate means
    df_mean_l3 <- mean(df_treat$NormValue) %>% as.data.frame
    ###calculate ci
    # bootstrapping with 1000 replications 
    results_l3 <- boot(data=df_treat$NormValue, statistic=Bmean, R=1000)
    # get 95% confidence interval 
    results_boot_l3 <- boot.ci(results_l3, conf = 0.90, 
                               type=c("norm", "basic", "perc", "bca"))
    ###create a data frame for ci
    df_ci_l3 <- results_boot_l3$bca %>% as.data.frame
    df_mean_l3$gene <- i
    df_mean_l3$treatment <- treat
    df_ci_l3$gene <- i
    df_ci_l3$treatment <- treat
    #combine the results
    Res_mean <- rbind(Res_mean, df_mean_l3)
    Res_ci <- rbind(Res_ci, df_ci_l3)
  }
}

######combine the two datasets
colnames(Res_mean)[1]<-'mean'
###select column from 4to7
Res_ci <- Res_ci[,4:7]
###rename the first column
colnames(Res_ci)[1]<-'lower'
colnames(Res_ci)[2]<-'upper'

######combine the two datasets
Res_em <- merge(Res_mean, Res_ci, by=c('gene', 'treatment'))

#####create color for the dataset
paleta_syncom <- c('#80cdc1','#018571')
names(paleta_syncom) <- c('NB','+SynCom')

######order the treatment
df_l3$treatment <- df_l3$treatment %>%
  factor(levels = c('NB', '+SynCom'))

Res_em$treatment <- Res_em$treatment %>%
  factor(levels = c('NB', '+SynCom'))

######combine gene and treatment
df_l3$group1 <- paste(df_l3$gene, df_l3$treatment, sep='_')
Res_em$group1 <- paste(Res_em$gene, Res_em$treatment, sep='_')
df_letter_l3_syncom$group1 <- paste(df_letter_l3_syncom$gene,
                                    df_letter_l3_syncom$treatment, sep='_')

######order the factors
df_l3$group1 <- df_l3$group1 %>%
  factor(levels = c('W22_NB', 'W22_+SynCom',
                    'zm00001eb022810_NB', 'zm00001eb022810_+SynCom',
                    'zm00001eb290350_NB','zm00001eb290350_+SynCom',
                    'zm00001eb149250_#1_NB','zm00001eb149250_#1_+SynCom',
                    'zm00001eb149250_#2_NB','zm00001eb149250_#2_+SynCom',
                    'zm00001eb149250_#3_NB','zm00001eb149250_#3_+SynCom'))

Res_em$group1 <- Res_em$group1 %>%
  factor(levels = c('W22_NB', 'W22_+SynCom',
                    'zm00001eb022810_NB', 'zm00001eb022810_+SynCom',
                    'zm00001eb290350_NB','zm00001eb290350_+SynCom',
                    'zm00001eb149250_#1_NB','zm00001eb149250_#1_+SynCom',
                    'zm00001eb149250_#2_NB','zm00001eb149250_#2_+SynCom',
                    'zm00001eb149250_#3_NB','zm00001eb149250_#3_+SynCom'))

df_letter_l3_syncom$group1 <- df_letter_l3_syncom$group1 %>%
  factor(levels = c('W22_NB', 'W22_+SynCom',
                    'zm00001eb022810_NB', 'zm00001eb022810_+SynCom',
                    'zm00001eb290350_NB','zm00001eb290350_+SynCom',
                    'zm00001eb149250_#1_NB','zm00001eb149250_#1_+SynCom',
                    'zm00001eb149250_#2_NB','zm00001eb149250_#2_+SynCom',
                    'zm00001eb149250_#3_NB','zm00001eb149250_#3_+SynCom'))

######plot the results
p_l3 <- ggplot(data=df_l3, aes(x=group1, y=NormValue,
                               color=treatment))+
  geom_point(aes(group=treatment),shape=21, size=1.5, alpha=1)+
  geom_pointrange(data=Res_em, aes(y=mean,
                                   ymin=lower, ymax=upper, 
                                   color= treatment),
                  size=1)+
  geom_text(data = df_letter_l3_syncom, aes(x = group1,y = 70,
                                            label = Letter),
            inherit.aes = F,size = 8,family ="Arial",color = "black") +
  ylab('Leaf length (cm)')+
  ggtitle('L3')+
  ylim(0,70)+
  scale_color_manual(values = paleta_syncom)+
  clean +
  theme(axis.text.x = element_text(size=20,vjust=1, angle = 45, hjust=1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, face='bold'),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size=20, hjust=0.5),
        legend.position = 'none')

p_l3

######add significance
p_l3 <- p_l3 +
  geom_line(data=tibble(x=c(1,2), y=c(64,64)), aes(x=x, y=y),
            inherit.aes = FALSE) +
  geom_text(data=tibble(x=1.5, y=65), aes(x=x, y=y, label='*'),
            inherit.aes = FALSE) +
  geom_line(data=tibble(x=c(3,4), y=c(64,64)), aes(x=x, y=y),
            inherit.aes = FALSE) +
  geom_text(data=tibble(x=3.5, y=65), aes(x=x, y=y, label='*'),
            inherit.aes = FALSE) +
  geom_line(data=tibble(x=c(5,6), y=c(64,64)), aes(x=x, y=y),
            inherit.aes = FALSE) +
  geom_text(data=tibble(x=5.5, y=65), aes(x=x, y=y, label='*'),
            inherit.aes = FALSE)

p_l3

######select l4
df_l4_syncom <- df %>%
  subset((leaf == 'L4') & (treatment == '+SynCom')) %>% droplevels

######apply Kruskal wallis
m1_l4_syncom <- kruskal.test(NormValue ~ gene,
                             data = df_l4_syncom)

m1_l4_syncom$p.value

DT_l4_syncom = dunnTest(NormValue ~ gene,
                        data = df_l4_syncom,method="bh")

PT_l4_syncom = DT_l4_syncom$res

# Use subset to select rows containing 'W22' in the Comparison column
selected_l4_syncom <- subset(PT_l4_syncom, grepl("W22", Comparison))

df_letter_l4_syncom <- cldList(P.adj ~ Comparison,
                               data = selected_l4_syncom,
                               remove.zero = F,
                               threshold = 0.1)

######add column 
df_letter_l4_syncom$comparison_syncom <- '+Syncom W22 vs +SynCom mutants'

######add a treatment column
df_letter_l4_syncom$treatment <- '+SynCom'

######add the leaf column
df_letter_l4_syncom$leaf <- 'L4'

######rename the first column
colnames(df_letter_l4_syncom)[1] <- 'gene'

######select l4
df_l4 <- df %>%
  subset(leaf == 'L4') %>% droplevels

######select the files
mtreat <- df_l4$treatment %>% as.character %>% unique
mgene <- df_l4$gene %>% as.character %>% unique

######create an empty 
Res_letter_l4 <- NULL

######perform the for loop for each leaf
for (i in mgene) {
  df_gene <- df_l4 %>% subset(gene == i) %>% droplevels
  m1 <- kruskal.test(NormValue ~ treatment, data = df_gene)
  DT = wilcox.test(NormValue ~ treatment, data = df_gene,method="bh")
  df_letter <- c(p.value = DT$p.value) %>% data.frame
  df_letter$leaf <- 'L4'
  df_letter$gene <- i
  Res_letter_l4 <- rbind(df_letter, Res_letter_l4)
}

######rename column 1
colnames(Res_letter_l4)[1] <- 'p.value'

######empty data frame
Res_mean <- NULL
Res_ci <- NULL

###perform the for loop
for (i in mgene) {
  df_gene <- df_l4 %>% subset(gene == i) %>% droplevels
  for (treat in mtreat) {
    df_treat <- df_gene %>% subset(treatment == treat) %>% droplevels
    #calculate means
    df_mean_l4 <- mean(df_treat$NormValue) %>% as.data.frame
    ###calculate ci
    # bootstrapping with 1000 replications 
    results_l4 <- boot(data=df_treat$NormValue, statistic=Bmean, R=1000)
    # get 95% confidence interval 
    results_boot_l4 <- boot.ci(results_l4, conf = 0.90, 
                               type=c("norm", "basic", "perc", "bca"))
    ###create a data frame for ci
    df_ci_l4 <- results_boot_l4$bca %>% as.data.frame
    df_mean_l4$gene <- i
    df_mean_l4$treatment <- treat
    df_ci_l4$gene <- i
    df_ci_l4$treatment <- treat
    #combine the results
    Res_mean <- rbind(Res_mean, df_mean_l4)
    Res_ci <- rbind(Res_ci, df_ci_l4)
  }
}

######combine the two datasets
colnames(Res_mean)[1]<-'mean'
###select column from 4to7
Res_ci <- Res_ci[,4:7]
###rename the first column
colnames(Res_ci)[1]<-'lower'
colnames(Res_ci)[2]<-'upper'

######combine the two datasets
Res_em <- merge(Res_mean, Res_ci, by=c('gene', 'treatment'))

#####create color for the dataset
paleta_syncom <- c('#80cdc1','#018571')
names(paleta_syncom) <- c('NB','+SynCom')

######order the treatment
df_l4$treatment <- df_l4$treatment %>%
  factor(levels = c('NB', '+SynCom'))

Res_em$treatment <- Res_em$treatment %>%
  factor(levels = c('NB', '+SynCom'))

######combine gene and treatment
df_l4$group1 <- paste(df_l4$gene, df_l4$treatment, sep='_')
Res_em$group1 <- paste(Res_em$gene, Res_em$treatment, sep='_')
df_letter_l4_syncom$group1 <- paste(df_letter_l4_syncom$gene,
                                    df_letter_l4_syncom$treatment, sep='_')

######order the factors
df_l4$group1 <- df_l4$group1 %>%
  factor(levels = c('W22_NB', 'W22_+SynCom',
                    'zm00001eb022810_NB', 'zm00001eb022810_+SynCom',
                    'zm00001eb290350_NB','zm00001eb290350_+SynCom',
                    'zm00001eb149250_#1_NB','zm00001eb149250_#1_+SynCom',
                    'zm00001eb149250_#2_NB','zm00001eb149250_#2_+SynCom',
                    'zm00001eb149250_#3_NB','zm00001eb149250_#3_+SynCom'))

Res_em$group1 <- Res_em$group1 %>%
  factor(levels = c('W22_NB', 'W22_+SynCom',
                    'zm00001eb022810_NB', 'zm00001eb022810_+SynCom',
                    'zm00001eb290350_NB','zm00001eb290350_+SynCom',
                    'zm00001eb149250_#1_NB','zm00001eb149250_#1_+SynCom',
                    'zm00001eb149250_#2_NB','zm00001eb149250_#2_+SynCom',
                    'zm00001eb149250_#3_NB','zm00001eb149250_#3_+SynCom'))

df_letter_l4_syncom$group1 <- df_letter_l4_syncom$group1 %>%
  factor(levels = c('W22_NB', 'W22_+SynCom',
                    'zm00001eb022810_NB', 'zm00001eb022810_+SynCom',
                    'zm00001eb290350_NB','zm00001eb290350_+SynCom',
                    'zm00001eb149250_#1_NB','zm00001eb149250_#1_+SynCom',
                    'zm00001eb149250_#2_NB','zm00001eb149250_#2_+SynCom',
                    'zm00001eb149250_#3_NB','zm00001eb149250_#3_+SynCom'))

######plot the results
p_l4 <- ggplot(data=df_l4, aes(x=group1, y=NormValue,
                               color=treatment))+
  geom_point(aes(group=treatment),shape=21, size=1.5, alpha=1)+
  geom_pointrange(data=Res_em, aes(y=mean,
                                   ymin=lower, ymax=upper, 
                                   color= treatment),
                  size=1)+
  geom_text(data = df_letter_l4_syncom, aes(x = group1,y = 70,
                                            label = Letter),
            inherit.aes = F,size = 8,family ="Arial",color = "black") +
  ylab('Leaf length (cm)')+
  ggtitle('L4')+
  ylim(0,70)+
  scale_color_manual(values = paleta_syncom)+
  clean +
  theme(axis.text.x = element_text(size=20,vjust=1, angle = 45, hjust=1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, face='bold'),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size=20, hjust=0.5),
        legend.position = 'none')

p_l4

######add significance
p_l4 <- p_l4 +
  geom_line(data=tibble(x=c(1,2), y=c(64,64)), aes(x=x, y=y),
            inherit.aes = FALSE) +
  geom_text(data=tibble(x=1.5, y=65), aes(x=x, y=y, label='*'),
            inherit.aes = FALSE) +
  geom_line(data=tibble(x=c(3,4), y=c(64,64)), aes(x=x, y=y),
            inherit.aes = FALSE) +
  geom_text(data=tibble(x=3.5, y=65), aes(x=x, y=y, label='*'),
            inherit.aes = FALSE) +
  geom_line(data=tibble(x=c(11,12), y=c(64,64)), aes(x=x, y=y),
            inherit.aes = FALSE) +
  geom_text(data=tibble(x=11.5, y=65), aes(x=x, y=y, label='*'),
            inherit.aes = FALSE)

p_l4

######select l5
df_l5_syncom <- df %>%
  subset((leaf == 'L5') & (treatment == '+SynCom')) %>% droplevels

######apply Kruskal wallis
m1_l5_syncom <- kruskal.test(NormValue ~ gene,
                             data = df_l5_syncom)

m1_l5_syncom$p.value

DT_l5_syncom = dunnTest(NormValue ~ gene,
                        data = df_l5_syncom,method="bh")

PT_l5_syncom = DT_l5_syncom$res

# Use subset to select rows containing 'W22' in the Comparison column
selected_l5_syncom <- subset(PT_l5_syncom, grepl("W22", Comparison))

df_letter_l5_syncom <- cldList(P.adj ~ Comparison,
                               data = selected_l5_syncom,
                               remove.zero = F,
                               threshold = 0.1)

######add column 
df_letter_l5_syncom$comparison_syncom <- '+Syncom W22 vs +SynCom mutants'

######add a treatment column
df_letter_l5_syncom$treatment <- '+SynCom'

######add the leaf column
df_letter_l5_syncom$leaf <- 'L5'

######rename the first column
colnames(df_letter_l5_syncom)[1] <- 'gene'

######select l5
df_l5 <- df %>%
  subset(leaf == 'L5') %>% droplevels

######select the files
mtreat <- df_l5$treatment %>% as.character %>% unique
mgene <- df_l5$gene %>% as.character %>% unique

######create an empty 
Res_letter_l5 <- NULL

######perform the for loop for each leaf
for (i in mgene) {
  df_gene <- df_l5 %>% subset(gene == i) %>% droplevels
  m1 <- kruskal.test(NormValue ~ treatment, data = df_gene)
  DT = wilcox.test(NormValue ~ treatment, data = df_gene,method="bh")
  df_letter <- c(p.value = DT$p.value) %>% data.frame
  df_letter$leaf <- 'L5'
  df_letter$gene <- i
  Res_letter_l5 <- rbind(df_letter, Res_letter_l5)
}

######rename column 1
colnames(Res_letter_l5)[1] <- 'p.value'

######empty data frame
Res_mean <- NULL
Res_ci <- NULL

###perform the for loop
for (i in mgene) {
  df_gene <- df_l5 %>% subset(gene == i) %>% droplevels
  for (treat in mtreat) {
    df_treat <- df_gene %>% subset(treatment == treat) %>% droplevels
    #calculate means
    df_mean_l5 <- mean(df_treat$NormValue) %>% as.data.frame
    ###calculate ci
    # bootstrapping with 1000 replications 
    results_l5 <- boot(data=df_treat$NormValue, statistic=Bmean, R=1000)
    # get 95% confidence interval 
    results_boot_l5 <- boot.ci(results_l5, conf = 0.90, 
                               type=c("norm", "basic", "perc", "bca"))
    ###create a data frame for ci
    df_ci_l5 <- results_boot_l5$bca %>% as.data.frame
    df_mean_l5$gene <- i
    df_mean_l5$treatment <- treat
    df_ci_l5$gene <- i
    df_ci_l5$treatment <- treat
    #combine the results
    Res_mean <- rbind(Res_mean, df_mean_l5)
    Res_ci <- rbind(Res_ci, df_ci_l5)
  }
}

######combine the two datasets
colnames(Res_mean)[1]<-'mean'
###select column from 4to7
Res_ci <- Res_ci[,4:7]
###rename the first column
colnames(Res_ci)[1]<-'lower'
colnames(Res_ci)[2]<-'upper'

######combine the two datasets
Res_em <- merge(Res_mean, Res_ci, by=c('gene', 'treatment'))

#####create color for the dataset
paleta_syncom <- c('#80cdc1','#018571')
names(paleta_syncom) <- c('NB','+SynCom')

######order the treatment
df_l5$treatment <- df_l5$treatment %>%
  factor(levels = c('NB', '+SynCom'))

Res_em$treatment <- Res_em$treatment %>%
  factor(levels = c('NB', '+SynCom'))

######combine gene and treatment
df_l5$group1 <- paste(df_l5$gene, df_l5$treatment, sep='_')
Res_em$group1 <- paste(Res_em$gene, Res_em$treatment, sep='_')
df_letter_l5_syncom$group1 <- paste(df_letter_l5_syncom$gene,
                                    df_letter_l5_syncom$treatment, sep='_')

######order the factors
df_l5$group1 <- df_l5$group1 %>%
  factor(levels = c('W22_NB', 'W22_+SynCom',
                    'zm00001eb022810_NB', 'zm00001eb022810_+SynCom',
                    'zm00001eb290350_NB','zm00001eb290350_+SynCom',
                    'zm00001eb149250_#1_NB','zm00001eb149250_#1_+SynCom',
                    'zm00001eb149250_#2_NB','zm00001eb149250_#2_+SynCom',
                    'zm00001eb149250_#3_NB','zm00001eb149250_#3_+SynCom'))

Res_em$group1 <- Res_em$group1 %>%
  factor(levels = c('W22_NB', 'W22_+SynCom',
                    'zm00001eb022810_NB', 'zm00001eb022810_+SynCom',
                    'zm00001eb290350_NB','zm00001eb290350_+SynCom',
                    'zm00001eb149250_#1_NB','zm00001eb149250_#1_+SynCom',
                    'zm00001eb149250_#2_NB','zm00001eb149250_#2_+SynCom',
                    'zm00001eb149250_#3_NB','zm00001eb149250_#3_+SynCom'))

df_letter_l5_syncom$group1 <- df_letter_l5_syncom$group1 %>%
  factor(levels = c('W22_NB', 'W22_+SynCom',
                    'zm00001eb022810_NB', 'zm00001eb022810_+SynCom',
                    'zm00001eb290350_NB','zm00001eb290350_+SynCom',
                    'zm00001eb149250_#1_NB','zm00001eb149250_#1_+SynCom',
                    'zm00001eb149250_#2_NB','zm00001eb149250_#2_+SynCom',
                    'zm00001eb149250_#3_NB','zm00001eb149250_#3_+SynCom'))

######plot the results
p_l5 <- ggplot(data=df_l5, aes(x=group1, y=NormValue,
                               color=treatment))+
  geom_point(aes(group=treatment),shape=21, size=1.5, alpha=1)+
  geom_pointrange(data=Res_em, aes(y=mean,
                                   ymin=lower, ymax=upper, 
                                   color= treatment),
                  size=1)+
  geom_text(data = df_letter_l5_syncom, aes(x = group1,y = 70,
                                            label = Letter),
            inherit.aes = F,size = 8,family ="Arial",color = "black") +
  ylab('Leaf length (cm)')+
  ggtitle('L5')+
  ylim(0,70)+
  scale_color_manual(values = paleta_syncom, name='Treatment')+
  clean +
  theme(axis.text.x = element_text(size=20,vjust=1, angle = 45, hjust=1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, face='bold'),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size=20, hjust=0.5),
        legend.position = 'right')

p_l5

######add significance
p_l5 <- p_l5 +
  geom_line(data=tibble(x=c(1,2), y=c(64,64)), aes(x=x, y=y),
            inherit.aes = FALSE) +
  geom_text(data=tibble(x=1.5, y=65), aes(x=x, y=y, label='*'),
            inherit.aes = FALSE) +
  geom_line(data=tibble(x=c(3,4), y=c(64,64)), aes(x=x, y=y),
            inherit.aes = FALSE) +
  geom_text(data=tibble(x=3.5, y=65), aes(x=x, y=y, label='*'),
            inherit.aes = FALSE) +
  geom_line(data=tibble(x=c(7,8), y=c(64,64)), aes(x=x, y=y),
            inherit.aes = FALSE) +
  geom_text(data=tibble(x=7.5, y=65), aes(x=x, y=y, label='*'),
            inherit.aes = FALSE)

p_l5

######composition
composition <- egg::ggarrange(p, p_l2, p_l3, p_l4, p_l5,
                              nrow=1, ncol = 5)

######save the data
oh.save.pdf(p = composition,
            outname = "figS5K_individual_leaf_growth_maize_mutants_final_withoutgene18.pdf",
            outdir = "../figures/",width = 40,height = 20)
