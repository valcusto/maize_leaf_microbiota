######Load packages
library(ohchibi)
library(dplyr)
library(egg)
library(ggrepel)

######Set work directory
setwd("C:/Users/Castrillo Lab/OneDrive - The University of Nottingham/Publication 2019/Paper preparation/manuscripts/scripts")

######helper functions
source('0_1-Permanova_plots.R')

######Open Dat object
Dat <- readRDS(file = "../cleandata/dat_leafionome_with_without_syncom.RDS")

#####Exploring the data
Tab <- Dat$Tab  %>% t

###Transform in z-score matrix 
#A z-score is a measure of how many standard deviations 
#(how far my observation is from the mean) below or above
#the population mean a raw score is.
Tab_z <- Tab %>% scale

###PCA
colMeans(Tab_z) %>% sort %>% plot

mpca_z <- prcomp(x = Tab_z,center = F,scale. = F,retx = T)
summary(mpca_z)

mpca_z$ rotation %>% head

mpca_z$x %>% head

#Transfrom into a dataframe
#This is just to create a data frame to plot using ggplot
df_x <- mpca_z$x %>% as.data.frame
df_x$UiD <- rownames(df_x)

summary(mpca_z)

#Merge df_x with the original metadata
merged <- merge(df_x,Dat$Map, by = "UiD")

#Clean graph
size=25
alpha=1
stroke=1.5
col_shape_background="white"
alpha_shape_background=0
lines_zero = TRUE
size_axis_line=2
type_axis_line = "longdash"
ratio_size_shape_background=1.3
y_vjust=0.5
x_hjust=0.5
size_axis_text=20
size_axis_title=30
size_legend_text=20
size_title_text = 30
legend_proportion_size=0.5
size_lines_panel = 0
size_panel_border = 2
font_family = "Helvetica"

clean <- theme(axis.line = element_blank(),
               panel.background = element_rect(fill = 'white'),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(fill=NA,color =  "black",size = size_panel_border),
               axis.ticks = element_line(colour = "black",size = 2.5),
               axis.text.x = element_text(family = font_family,size =size_axis_text,colour="black",hjust = x_hjust),
               axis.text.y = element_text(family = font_family,size=size_axis_text,colour="black",vjust = y_vjust),
               axis.title.x = element_text(family = font_family,size = size_axis_title,colour = "black"),
               axis.title.y = element_text(family = font_family,size=size_axis_title,colour="black"),
               legend.background = element_blank(),legend.key.size = unit(legend_proportion_size,"line"),
               legend.title=element_text(size=size_title_text,
                                         family = font_family,face = "bold",colour = "black"),            
               #legend.title = element_blank(),
               legend.key = element_blank(),
               legend.text = element_text(size=size_legend_text,
                                          family = font_family,face = "bold",colour = "black"),
               legend.position ="right")

# function that draws key without gap
draw_key_polygon2 <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)
  
  grid::rectGrob(
    width = grid::unit(1, "npc"),
    height = grid::unit(1, "npc"),
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}

# register new key drawing function, effect is global and persistent
# throughout R session!
GeomBar$draw_key = draw_key_polygon2

#Define color for number of leaf
paleta_leaf <- c('#d01c8b','#f1b6da','#7b3294',
                 '#c2a5cf','#a6dba0','#008837')
names(paleta_leaf) <- c('L1', 'L2', 'L3', 'L4', 'L5', 'L6')

paleta_nutrient <- c('#762a83','#1b7837')
names(paleta_nutrient) <- c('Full Hoagland','1/4 Hoagland')

######plot the results
pleaf <- ggplot(data = merged,mapping = aes(x = PC1, y = PC2), color='black') +
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_point(aes(color = NumLeaf,shape = Treatment),size = 10)+
  xlab(label = "PC1 47.21% Variance Explained") +
  ylab(label = "PC2 19.21% Variance Explained")+
  clean +
  scale_color_manual(values = paleta_leaf)+
  labs(color = "Leaf number",
       shape = 'Treatment')+
  theme(legend.position = 'none')

pleaf

pleaf_leg <- ggplot(data = merged,mapping = aes(x = PC1, y = PC2), color='black') +
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_point(aes(color = NumLeaf,shape = Treatment),size = 10)+
  xlab(label = "PC1 47.21% Variance Explained") +
  ylab(label = "PC2 19.21% Variance Explained")+
  clean +
  scale_color_manual(values = paleta_leaf)+
  labs(color = "Leaf number",
       shape = 'Treatment')+
  theme(legend.position = 'right')

pleaf_leg

pnutrient <- ggplot(data = merged,mapping = aes(x = PC1, y = PC2), color='black') +
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_point(aes(color = Nutrient,shape = Treatment),size = 10)+
  xlab(label = "PC1 47.21% Variance Explained") +
  ylab(label = "PC2 19.21% Variance Explained")+
  clean +
  scale_color_manual(values = paleta_nutrient)+
  labs(color = "Nutrient",
       shape = 'Treatment')+
  theme(legend.position = 'none')

pnutrient

pnutrient_leg <- ggplot(data = merged,mapping = aes(x = PC1, y = PC2), color='black') +
  geom_vline(xintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_hline(yintercept = 0,size = 2,color = "#D9D9D9",linetype = "longdash")+
  geom_point(aes(color = Nutrient,shape = Treatment),size = 10)+
  xlab(label = "PC1 47.21% Variance Explained") +
  ylab(label = "PC2 19.21% Variance Explained")+
  clean +
  scale_color_manual(values = paleta_nutrient)+
  labs(color = "Nutrient",
       shape = 'Treatment')+
  theme(legend.position = 'right')

pnutrient_leg

#Quantitative way of measuring importance of variables
#PERMANOVA 
Tab_z %>% head

Tab_dist <- dist(x = Tab_z,method = "euclidean")

mypermanova <- adonis2(formula = Tab_dist ~  Treatment+NumLeaf+Nutrient,
                      data = Dat$Map,permutations = 9999)
res <- kretxeu.permanova(mypermanova = mypermanova)

pperm <- res$p + 
  scale_fill_manual(labels = c("Residual", "Nutrient","Leaf number"),
                    values = c("#7f3b08","#fee0b6","white") %>% rev)+
  guides(fill = guide_legend(
   override.aes = list(color = "black", size = 0.5)))+
  theme(axis.title.x =element_blank(),
        legend.position = "none")

pperm_leg <- res$p + 
  scale_fill_manual(labels = c("Residual", "Nutrient","Leaf number"),
                    values = c("#7f3b08","#fee0b6","white") %>% rev)+
  guides(fill = guide_legend(
    override.aes = list(color = "black", size = 0.5)))+
  theme(axis.title.x =element_blank(),
        legend.position = "left")

######composition
composition_noleg <- egg::ggarrange(pperm,pleaf,pnutrient, 
                                  nrow = 1,widths = c(0.05,1,1))

composition_leg <- egg::ggarrange(pperm_leg,pleaf_leg,pnutrient_leg, 
                                  nrow = 1,widths = c(0.05,1,1))

######save figures
oh.save.pdf(p = composition_noleg,outname = "figS3j-leaf_ionome_with_without_bacteria_nutrient_noleg.pdf",
            outdir = "../figures/",width = 20,height = 15)

oh.save.pdf(p = composition_leg,outname = "figS3H-leaf_ionome_with_without_bacteria_nutrient_leg.pdf",
            outdir = "../figures/",width = 30,height = 20)
