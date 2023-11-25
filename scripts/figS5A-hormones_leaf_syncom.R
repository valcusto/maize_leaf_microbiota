######load library
library(ohchibi)
library(ggbreak)

######set work directory
setwd("C:/Users/Castrillo Lab/OneDrive - The University of Nottingham/Publication 2019/Paper preparation/manuscripts/scripts")

######helper functions
source('0-Clean_up_plots.R')
# function to obtain the mean
Bmean <- function(data, indices) {
  d <- data[indices] # allows boot to select sample 
  return(mean(d))
}

######create color for the dataset
paleta_marker <- c('#80cdc1','#018571')
names(paleta_marker) <- c('NB', '+SynCom')

######load data
df <- read.csv('../rawdata/maize_hormonomics.csv')

######correct the UiD
colnames(df)[1] <- 'UiD'

######select only the full syncom and bacteria
df_sub <- df %>%
  subset((Marker != 'L3') & (Marker != 'L4_up') & (Marker != 'L5')) %>% droplevels

###rename the treatment column
df_sub$Treatment[which(df_sub$Treatment == 'Full syncom')] <- '+SynCom'
df_sub$Treatment[which(df_sub$Treatment == 'No Bacteria')] <- 'NB'

df_sub$Marker[which(df_sub$Marker == 'Full syncom')] <- '+SynCom'
df_sub$Marker[which(df_sub$Marker == 'No Bacteria')] <- 'NB'

####organize the markers
df_sub$Marker <- df_sub$Marker %>%
  factor(levels=c('NB', '+SynCom'))

######select the ABA hormone
df_aba <- df_sub %>%
  subset(hormone_class == 'ABA') %>% droplevels

###### check the ABA
ggplot(data = df_aba,aes(x = Marker,y=value, color = Marker)) +
  geom_boxplot() + geom_point() +
  theme_ohchibi()

######find influential points on the dataset
model_aba <- lm(value ~ Marker, data = df_aba)

#find Cook's distance for each observation in the dataset
cooksD_aba <- cooks.distance(model_aba)

# Plot Cook's Distance with a horizontal line at 4/n to see which observations
#exceed this thresdhold
n_aba <- nrow(df_aba)
plot(cooksD_aba, main = "Cooks Distance for Influential Obs")
abline(h = 4/n_aba, lty = 2, col = "steelblue") # add cutoff line

#identify influential points
influential_aba <- as.numeric(names(cooksD_aba)[(cooksD_aba > (4/n_aba))])

#define new data frame with influential points removed
df_removed_aba <- df_aba[-influential_aba, ]

######no influential point
ggplot(data = df_removed_aba,aes(x = Marker,y=value, color = Marker)) +
  geom_boxplot() + geom_point() +
  #ylim(0,0.2)+
  theme_ohchibi()

###1.interdependence of the variables
###2.variance homogeneity
leveneTest(value ~ Marker, data = df_removed_aba)
#results: p-value>0.05. Do not reject the null hypothesis.
#All the groups have similar variance

###3.normality
qqnorm(df_removed_aba$value)
qqline(df_removed_aba$value, col='red')
hist(df_removed_aba$value)

#######Fit the ANOVA model (one-way model)
m1_aba <- aov(data = df_removed_aba, 
               formula = value ~ Marker)
summary(m1_aba)

######Elaborate the letter data frame
df_em_aba <- emmeans(m1_aba,specs = "Marker")  %>% as.data.frame

######add the letters 
df_em_aba$letters <- 'a'

######define the paleta for aba
paleta_aba <- c('#c51b7d', '#8e0152')
names(paleta_aba) <- c('NB', '+SynCom')

######plot the result
p_aba <- ggplot(data=df_removed_aba, aes(x=Marker, y=value, 
                                    color=Marker))+
  geom_point(shape=21, size=4, alpha=1.5)+
  geom_pointrange(data=df_em_aba, aes(y = emmean, ymin = lower.CL, 
                                       ymax = upper.CL, 
                                       color=Marker), size=1.5)+
  geom_text(data = df_em_aba, aes(x = Marker,y = 30,
                                   label = letters),
            inherit.aes = F,size = 6,family ="Arial",color = "black") +
  scale_color_manual(values = paleta_aba)+
  ggtitle('ABA') +
  xlab(NULL)+ 
  ylab(NULL)+
  ylim(0,30)+
  clean +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size=30, face = 'bold'))

p_aba

######select the iaa hormone
df_iaa <- df_sub %>%
  subset(hormone_class == 'IAA') %>% droplevels

###### check the iaa
ggplot(data = df_iaa,aes(x = Marker,y=value, color = Marker)) +
  geom_boxplot() + geom_point() +
  theme_ohchibi()

######find influential points on the dataset
model_iaa <- lm(value ~ Marker, data = df_iaa)

#find Cook's distance for each observation in the dataset
cooksD_iaa <- cooks.distance(model_iaa)

# Plot Cook's Distance with a horizontal line at 4/n to see which observations
#exceed this thresdhold
n_iaa <- nrow(df_iaa)
plot(cooksD_iaa, main = "Cooks Distance for Influential Obs")
abline(h = 4/n_iaa, lty = 2, col = "steelblue") # add cutoff line

#identify influential points
influential_iaa <- as.numeric(names(cooksD_iaa)[(cooksD_iaa > (4/n_iaa))])

###1.interdependence of the variables
###2.variance homogeneity
leveneTest(value ~ Marker, data = df_iaa)
#results: p-value>0.05. Do not reject the null hypothesis.
#All the groups have similar variance

###3.normality
qqnorm(df_iaa$value)
qqline(df_iaa$value, col='red')
hist(df_iaa$value)

#######Fit the ANOVA model (one-way model)
m1_iaa <- aov(data = df_iaa, 
              formula = value ~ Marker)
summary(m1_iaa)

######Elaborate the letter data frame
df_em_iaa <- emmeans(m1_iaa,specs = "Marker")  %>% as.data.frame

######add the letters
df_em_iaa$letters <- 'a'

######define the paleta for iaa
paleta_iaa <- c('#fee0b6', '#fdb863')
names(paleta_iaa) <- c('NB', '+SynCom')

######plot the result
p_iaa <- ggplot(data=df_iaa, aes(x=Marker, y=value, 
                                 color=Marker))+
  geom_point(shape=21, size=4, alpha=1.5)+
  geom_pointrange(data=df_em_iaa, aes(y = emmean, ymin = lower.CL, 
                                      ymax = upper.CL, 
                                      color=Marker), size=1.5)+
  geom_text(data = df_em_iaa, aes(x = Marker,y = 10,
                                  label = letters),
            inherit.aes = F,size = 6,family ="Arial",color = "black") +
  scale_color_manual(values = paleta_iaa)+
  ggtitle('IAA') +
  xlab(NULL)+ 
  ylab(NULL)+
  ylim(0,10)+
  clean +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size=30, face = 'bold'))

p_iaa

######select the ip hormone
df_ip <- df_sub %>%
  subset(hormones_name == 'iP') %>% droplevels

###### cheip the ip
ggplot(data = df_ip,aes(x = Marker,y=value, color = Marker)) +
  geom_boxplot() + geom_point() +
  theme_ohchibi()

######find influential points on the dataset
model_ip <- lm(value ~ Marker, data = df_ip)

#find Cook's distance for each observation in the dataset
cooksD_ip <- cooks.distance(model_ip)

# Plot Cook's Distance with a horizontal line at 4/n to see which observations
#exceed this thresdhold
n_ip <- nrow(df_ip)
plot(cooksD_ip, main = "Cooks Distance for Influential Obs")
abline(h = 4/n_ip, lty = 2, col = "steelblue") # add cutoff line

#identify influential points
influential_ip <- as.numeric(names(cooksD_ip)[(cooksD_ip > (4/n_ip))])

#define new data frame with influential points removed
df_removed_ip <- df_ip[-influential_ip, ]

######no influential point
ggplot(data = df_removed_ip,aes(x = Marker,y=value, color = Marker)) +
  geom_boxplot() + geom_point() +
  #ylim(0,0.2)+
  theme_ohchibi()

###1.interdependence of the variables
###2.variance homogeneity
leveneTest(value ~ Marker, data = df_removed_ip)
#results: p-value>0.05. Do not reject the null hypothesis.
#All the groups have similar variance

###3.normality
qqnorm(df_removed_ip$value)
qqline(df_removed_ip$value, col='red')
hist(df_removed_ip$value)

#######Fit the ANOVA model (one-way model)
m1_ip <- aov(data = df_removed_ip, 
             formula = value ~ Marker)
summary(m1_ip)

######Elaborate the letter data frame
df_em_ip <- emmeans(m1_ip,specs = "Marker")  %>% as.data.frame

#####add letters
df_em_ip$letters <- 'a'

######define the paleta for ip
paleta_ip <- c('#d8daeb', '#b2abd2')
names(paleta_ip) <- c('NB', '+SynCom')

######plot the result
p_ip <- ggplot(data=df_ip, aes(x=Marker, y=value, 
                               color=Marker))+
  geom_point(shape=21, size=4, alpha=1.5)+
  geom_pointrange(data=df_em_ip, aes(y = emmean, ymin = lower.CL, 
                                     ymax = upper.CL, 
                                     color=Marker), size=1.5)+
  geom_text(data = df_em_ip, aes(x = Marker,y = 0.01,
                                 label = letters),
            inherit.aes = F,size = 6,family ="Arial",color = "black") +
  scale_color_manual(values = paleta_ip)+
  ggtitle('iP') +
  xlab(NULL)+ 
  ylab('Hormone concentration (ngH/g FW)')+
  ylim(0,0.010)+
  clean +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size=30, face = 'bold'))

p_ip

######select the tz hormone
df_tz <- df_sub %>%
  subset(hormones_name == 'tZ') %>% droplevels

###### chetz the tz
ggplot(data = df_tz,aes(x = Marker,y=value, color = Marker)) +
  geom_boxplot() + geom_point() +
  theme_ohchibi()

######find influential points on the dataset
model_tz <- lm(value ~ Marker, data = df_tz)

#find Cook's distance for each observation in the dataset
cooksD_tz <- cooks.distance(model_tz)

# Plot Cook's Distance with a horizontal line at 4/n to see which observations
#exceed this thresdhold
n_tz <- nrow(df_tz)
plot(cooksD_tz, main = "Cooks Distance for Influential Obs")
abline(h = 4/n_tz, lty = 2, col = "steelblue") # add cutoff line

#identify influential points
influential_tz <- as.numeric(names(cooksD_tz)[(cooksD_tz > (4/n_tz))])

#define new data frame with influential points removed
df_removed_tz <- df_tz[-influential_tz, ]

######no influential point
ggplot(data = df_removed_tz,aes(x = Marker,y=value, color = Marker)) +
  geom_boxplot() + geom_point() +
  #ylim(0,0.2)+
  theme_ohchibi()

###1.interdependence of the variables
###2.variance homogeneity
leveneTest(value ~ Marker, data = df_removed_tz)
#results: p-value>0.05. Do not reject the null hypothesis.
#All the groups have similar variance

###3.normality
qqnorm(df_removed_tz$value)
qqline(df_removed_tz$value, col='red')
hist(df_removed_tz$value)

#######Fit the ANOVA model (one-way model)
m1_tz <- aov(data = df_removed_tz, 
             formula = value ~ Marker)
summary(m1_tz)

######Elaborate the letter data frame
df_em_tz <- emmeans(m1_tz,specs = "Marker")  %>% as.data.frame

######add the letters
df_em_tz$letters <- 'a'

######define the paleta for tz
paleta_tz <- c('#807dba', '#6a51a3')
names(paleta_tz) <- c('NB', '+SynCom')

######plot the result
p_tz <- ggplot(data=df_tz, aes(x=Marker, y=value, 
                               color=Marker))+
  geom_point(shape=21, size=4, alpha=1.5)+
  geom_pointrange(data=df_em_tz, aes(y = emmean, ymin = lower.CL, 
                                     ymax = upper.CL, 
                                     color=Marker), size=1.5)+
  geom_text(data = df_em_tz, aes(x = Marker,y = 0.35,
                                 label = letters),
            inherit.aes = F,size = 6,family ="Arial",color = "black") +
  scale_color_manual(values = paleta_tz)+
  ggtitle('tZ') +
  xlab(NULL)+ 
  ylab(NULL)+
  ylim(0,0.35)+
  clean +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size=30, face = 'bold'))

p_tz

######select the dhz hormone
df_dhz <- df_sub %>%
  subset(hormones_name == 'DHZ') %>% droplevels

###### chedhz the dhz
ggplot(data = df_dhz,aes(x = Marker,y=value, color = Marker)) +
  geom_boxplot() + geom_point() +
  theme_ohchibi()

######find influential points on the dataset
model_dhz <- lm(value ~ Marker, data = df_dhz)

#find Cook's distance for each observation in the dataset
cooksD_dhz <- cooks.distance(model_dhz)

# Plot Cook's Distance with a horizontal line at 4/n to see which observations
#exceed this thresdhold
n_dhz <- nrow(df_dhz)
plot(cooksD_dhz, main = "Cooks Distance for Influential Obs")
abline(h = 4/n_dhz, lty = 2, col = "steelblue") # add cutoff line

#identify influential points
influential_dhz <- as.numeric(names(cooksD_dhz)[(cooksD_dhz > (4/n_dhz))])

#define new data frame with influential points removed
df_removed_dhz <- df_dhz %>%
  subset(UiD != 28) %>% droplevels

######no influential point
ggplot(data = df_removed_dhz,aes(x = Marker,y=value, color = Marker)) +
  geom_boxplot() + geom_point() +
  #ylim(0,0.2)+
  theme_ohchibi()

###1.interdependence of the variables
###2.variance homogeneity
leveneTest(value ~ Marker, data = df_removed_dhz)
#results: p-value>0.05. Do not reject the null hypothesis.
#All the groups have similar variance

###3.normality
qqnorm(df_removed_dhz$value)
qqline(df_removed_dhz$value, col='red')
hist(df_removed_dhz$value)

#######Fit the ANOVA model (one-way model)
m1_dhz <- aov(data = df_removed_dhz, 
              formula = value ~ Marker)
summary(m1_dhz)

######Elaborate the letter data frame
df_em_dhz <- emmeans(m1_dhz,specs = "Marker")  %>% as.data.frame

######add the letters
df_em_dhz$letters <- 'a'

######define the paleta for dhz
paleta_dhz <- c('#54278f', '#3f007d')
names(paleta_dhz) <- c('NB', '+SynCom')

######plot the result
p_dhz <- ggplot(data=df_removed_dhz, aes(x=Marker, y=value, 
                                         color=Marker))+
  geom_point(shape=21, size=4, alpha=1.5)+
  geom_pointrange(data=df_em_dhz, aes(y = emmean, ymin = lower.CL, 
                                      ymax = upper.CL, 
                                      color=Marker), size=1.5)+
  geom_text(data = df_em_dhz, aes(x = Marker,y = 0.03,
                                  label = letters),
            inherit.aes = F,size = 6,family ="Arial",color = "black") +
  scale_color_manual(values = paleta_dhz)+
  ggtitle('DHZ') +
  xlab(NULL)+ 
  ylab(NULL)+
  ylim(0,0.03)+
  clean +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size=30, face = 'bold'))

p_dhz

######select the ga4 hormone
df_ga4 <- df_sub %>%
  subset(hormones_name == 'GA4') %>% droplevels

###### chega4 the ga4
ggplot(data = df_ga4,aes(x = Marker,y=value, color = Marker)) +
  geom_boxplot() + geom_point() +
  theme_ohchibi()

######find influential points on the dataset
model_ga4 <- lm(value ~ Marker, data = df_ga4)

#find Cook's distance for each observation in the dataset
cooksD_ga4 <- cooks.distance(model_ga4)

# Plot Cook's Distance with a horizontal line at 4/n to see which observations
#exceed this thresdhold
n_ga4 <- nrow(df_ga4)
plot(cooksD_ga4, main = "Cooks Distance for Influential Obs")
abline(h = 4/n_ga4, lty = 2, col = "steelblue") # add cutoff line

#identify influential points
influential_ga4 <- as.numeric(names(cooksD_ga4)[(cooksD_ga4 > (4/n_ga4))])

#define new data frame with influential points removed
df_removed_ga4 <- df_ga4 %>%
  subset(UiD != 178) %>% droplevels

######no influential point
ggplot(data = df_removed_ga4,aes(x = Marker,y=value, color = Marker)) +
  geom_boxplot() + geom_point() +
  #ylim(0,0.2)+
  theme_ohchibi()

###1.interdependence of the variables
###2.variance homogeneity
leveneTest(value ~ Marker, data = df_removed_ga4)
#results: p-value>0.05. Do not reject the null hypothesis.
#All the groups have similar variance

###3.normality
qqnorm(df_removed_ga4$value)
qqline(df_removed_ga4$value, col='red')
hist(df_removed_ga4$value)

#######Fit the ANOVA model (one-way model)
m1_ga4 <- aov(data = df_removed_ga4, 
              formula = value ~ Marker)
summary(m1_ga4)

######Elaborate the letter data frame
df_em_ga4 <- emmeans(m1_ga4,specs = "Marker")  %>% as.data.frame

######add letters
df_em_ga4$letters <- 'a'

######define the paleta for ga4
paleta_ga4 <- c('#1b7837', '#00441b')
names(paleta_ga4) <- c('NB', '+SynCom')

######plot the result
p_ga4 <- ggplot(data=df_removed_ga4, aes(x=Marker, y=value, 
                                         color=Marker))+
  geom_point(shape=21, size=4, alpha=1.5)+
  geom_pointrange(data=df_em_ga4, aes(y = emmean, ymin = lower.CL, 
                                      ymax = upper.CL, 
                                      color=Marker), size=1.5)+
  geom_text(data = df_em_ga4, aes(x = Marker,y = 0.035,
                                  label = letters),
            inherit.aes = F,size = 6,family ="Arial",color = "black") +
  scale_color_manual(values = paleta_ga4)+
  ggtitle('GA4') +
  xlab(NULL)+ 
  ylab(NULL)+
  ylim(0,0.035)+
  clean +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size=30, face = 'bold'))

p_ga4

######select the ga12 hormone
df_ga12 <- df_sub %>%
  subset(hormones_name == 'GA12') %>% droplevels

###### chega12 the ga12
ggplot(data = df_ga12,aes(x = Marker,y=value, color = Marker)) +
  geom_boxplot() + geom_point() +
  theme_ohchibi()

######find influential points on the dataset
model_ga12 <- lm(value ~ Marker, data = df_ga12)

#find Cook's distance for each observation in the dataset
cooksD_ga12 <- cooks.distance(model_ga12)

# Plot Cook's Distance with a horizontal line at 12/n to see which observations
#exceed this thresdhold
n_ga12 <- nrow(df_ga12)
plot(cooksD_ga12, main = "Cooks Distance for Influential Obs")
abline(h = 4/n_ga12, lty = 2, col = "steelblue") # add cutoff line

#identify influential points
influential_ga12 <- as.numeric(names(cooksD_ga12)[(cooksD_ga12 > (4/n_ga12))])

#define new data frame with influential points removed
df_removed_ga12 <- df_ga12 %>%
  subset(UiD != 88) %>% droplevels

######no influential point
ggplot(data = df_removed_ga12,aes(x = Marker,y=value, color = Marker)) +
  geom_boxplot() + geom_point() +
  #ylim(0,0.2)+
  theme_ohchibi()

###1.interdependence of the variables
###2.variance homogeneity
leveneTest(value ~ Marker, data = df_removed_ga12)
#results: p-value>0.05. Do not reject the null hypothesis.
#All the groups have similar variance

###3.normality
qqnorm(df_removed_ga12$value)
qqline(df_removed_ga12$value, col='red')
hist(df_removed_ga12$value)

#######Fit the ANOVA model (one-way model)
m1_ga12 <- aov(data = df_removed_ga12, 
               formula = value ~ Marker)
summary(m1_ga12)

######Elaborate the letter data frame
df_em_ga12 <- emmeans(m1_ga12,specs = "Marker")  %>% as.data.frame

######add the letters
df_em_ga12$letters <- 'a'

######define the paleta for ga12
paleta_ga12 <- c('#a6dba0', '#5aae61')
names(paleta_ga12) <- c('NB', '+SynCom')

######plot the result
p_ga12 <- ggplot(data=df_removed_ga12, aes(x=Marker, y=value, 
                                           color=Marker))+
  geom_point(shape=21, size=4, alpha=1.5)+
  geom_pointrange(data=df_em_ga12, aes(y = emmean, ymin = lower.CL, 
                                       ymax = upper.CL, 
                                       color=Marker), size=1.5)+
  geom_text(data = df_em_ga12, aes(x = Marker,y = 2.3,
                                   label = letters),
            inherit.aes = F,size = 6,family ="Arial",color = "black") +
  scale_color_manual(values = paleta_ga12)+
  ggtitle('GA12') +
  xlab(NULL)+ 
  ylab(NULL)+
  ylim(0,2.3)+
  clean +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size=30, face = 'bold'))

p_ga12

#Create blank plot
p_blank <- ggplot() +theme_void()

######composition
composition <- egg::ggarrange(p_aba, p_iaa, p_blank,
                              p_ip, p_tz, p_dhz,
                              p_ga4, p_ga12,p_blank,
                              nrow = 3, ncol = 3)

######save the figure
oh.save.pdf(p = composition,
            outname = "figS5A_hormonomics_nb_syncom.pdf",
            outdir = "../figures/", width = 15,height = 15)
