######load packages
library(ohchibi)
library(multcomp)
library(emmeans)
library(dplyr)
library(egg)
library(ggrepel)

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

######Open the csv file
root <- read.csv(file = "../rawdata/Root_dry_wgt.csv")

######Transform some variables in factor
root$treatment <- root$treatment %>% factor
root$experiment <- root$experiment %>% factor
root$genotype <- root$genotype %>% 
  factor(levels =c("B73", "Matuba", "CVSt", "Zm523"))
######correct sao domingos name
root$soil <- root$soil %>%
  gsub(pattern = 'Sao Domingos', replacement = 'S?o Domingos')
root$soil <- root$soil %>% 
  factor(levels = c("S?o Domingos","Tarrafal", "Sutton Bonington"))
root$tray <- root$tray %>% factor
root$municipalities <- root$municipalities %>% factor
root$region <- root$region %>% factor

######select the 50%FC
root_50 <- root %>%
  subset(treatment == '50%FC') %>% droplevels

###Paletta soil
paleta_soil<- c("#FF7F00","#A6CEE3",'#1F78B4')
names(paleta_soil) <- c( "S?o Domingos", "Tarrafal",
                        "Sutton Bonington")

######paleta genotype
paleta_geno <- c('#e66101', '#fdb863', '#b2abd2', '#5e3c99')
names(paleta_geno) <- c('B73', 'Matuba', 'CVSt', 'Zm523')

######Explore the root dry weight
ggplot(data = root_50, mapping = aes(x = genotype, y = rdw_mg))+
  geom_boxplot() +
  geom_point()+
  facet_grid(.~soil)+
  scale_y_continuous(limits = c(0,800))

######statistical assumptions
#1.independence of the variables

#2.variance homogeneity
#perform levene test
leveneTest(rdw_mg ~ genotype, data = root_50)
#p<0.05. Do not reject null hypothesis.
#Soil with similar variance
#3.check the normality of the data
qqnorm(root_50$rdw_mg)
qqline(root_50$rdw_mg, col='red')

######Linear model
m1 <- aov(formula = rdw_mg ~ genotype + soil+tray+experiment,
          data = root_50)
summary(m1)

#mcp meand multiple comparisons
msum <- glht(m1, linfct=mcp(genotype="Tukey")) %>% summary 
df_geno <- emmeans(m1,specs = "genotype")  %>% as.data.frame
df_geno$p.value <-  c(msum$test$pvalues[1:4])
cld_irri <- cld(msum)

cld_loc <- cld_irri$mcletters$Letters %>% as.data.frame
cld_loc$genotype <- rownames(cld_loc)
rownames(cld_loc) <- NULL
colnames(cld_loc)[1] <- 'letters'

######combine df_geno and cld
df_geno <- merge(df_geno, cld_loc, by='genotype')

######plot the result
rdw <- ggplot(data = root_50,aes(x = genotype,y =rdw_mg, color=genotype)) +
  geom_point(shape = 21,size = 1)+
  geom_text(data = df_geno, aes(x = genotype,y = 900,
                                 label = letters),
            inherit.aes = F,size = 6,family ="Arial",color = "black") +
  geom_pointrange(data = df_geno,
                  aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
                  size = 0.5) + 
  ylab('Root dry weight (mg)')+
  clean +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        axis.title.x = element_blank(),
        legend.position = 'right')+
  scale_color_manual(values = paleta_geno)

rdw

######multiple comprarison for soil
msum <- glht(m1, linfct=mcp(soil="Tukey")) %>% summary 
df_soil <- emmeans(m1,specs = "soil")  %>% as.data.frame
df_soil$p.value <-  c(msum$test$pvalues[1:3])

######display letterrs
cld_soil <- cld(msum)

cld_soil <- cld_soil$mcletters$Letters %>% as.data.frame
cld_soil$soil <- rownames(cld_soil)
rownames(cld_soil) <- NULL
colnames(cld_soil)[1] <- 'letters'

######combine df_geno and cld
df_soil <- merge(df_soil, cld_soil, by='soil')

######plot the result
rdw_soil <- ggplot(data = root_50,aes(x = soil,y =rdw_mg, color=soil,
                                      shape=soil)) +
  geom_point(shape = 21,size = 4)+
  geom_text(data = df_soil, aes(x = soil,y = 900,
                                label = letters),
            inherit.aes = F,size = 6,family ="Arial",color = "black") +
  geom_pointrange(data = df_soil,
                  aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
                  size = 1.5) + 
  ylab('Root dry weight (mg)')+
  clean +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 30),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 30),
        plot.title = element_text(size = 30, hjust = 0.5, face = 'bold'))+
  scale_color_manual(values = paleta_soil)

rdw_soil

######save the plot as pdf
oh.save.pdf(p = rdw,outname = "figS1i-root_dry_wgt_genotype.pdf",
            outdir = "../figures/",width = 15,height = 10)

oh.save.pdf(p = rdw_soil,outname = "figS1i-root_dry_wgt_soil.pdf",
            outdir = "../figures/",width = 15,height = 10)
