######load library
library(ohchibi)

######set seed
set.seed(130816)

######set work directory
setwd("C:/Users/Castrillo Lab/OneDrive - The University of Nottingham/Publication 2019/Paper preparation/manuscripts/scripts")

######helper functions
source('0_4-subset_Dataset.R')
source('0_1-Permanova_plots.R')
source('0-Clean_up_plots.R')
distfun <- function(x,method) vegan::vegdist(x = x, method = "bray")
size_legend_text <- 45
size_axis_title <- 35
size_legend_title <- 55
size_title_text <- 30
legend_proportion_size <- 4

######paleta fractions
paleta_fra <- c('#008837', '#fec44f', '#d95f0e')
names(paleta_fra) <- c('Shoot', 'Root', 'Substrate')

######paleta leaf
paleta_leaf <- c('#d01c8b','#f1b6da','#7b3294',
                 '#c2a5cf','#a6dba0','#008837')
names(paleta_leaf) <- c('L1', 'L2', 'L3',
                        'L4', 'L5', 'L6')

#### Using the full dataset 
Dat <- readRDS(file = "../cleandata/dat_syncomleaves_complete.end.RDS")

######select the relative abundance dataset
Dat_rar <- Dat$RelativeAbundance

######select only bacteria
Dat_syncom <- Dat_rar %>%
  subset.Dataset(Bacteria != "NB",drop = T,clean = T)

#######check the all sampkes together
Tab_bray_syncom <- distfun(t(Dat_syncom$Tab))

######calculate the permanova
mypermanova_syncom <- adonis2(Tab_bray_syncom ~  Fraction +Rep + LDepth+Nutrient,
                           data = Dat_syncom$Map,
                           permutations = 9999)
mypermanova_syncom

######plot mypermanova
p_perm_syncom <- kretxeu.permanova(mypermanova = mypermanova_syncom) %$% p + 
  scale_fill_manual(values = c("#f16913","#f7f7f7") %>% rev) +
  xlab(NULL) +
  theme(legend.key = element_rect(color = "black"),
        legend.title=element_text(size=size_title_text),
        legend.position = "left",
        panel.border = element_rect(fill=NA,color =  "black",
                                    size = 1))

######cap analysis
mcap <- oh.cap(Tab = Dat_syncom$Tab,Map = Dat_syncom$Map,
               formula = "Fraction + Condition(Nutrient) +Condition(Rep) + Condition(LDepth)",
               distfun = distfun,perms = 9999,sqrt = T)
mcap$perm_anova_terms

######caps plot
chibi.cap(list_ohpco = mcap,col_val = "Fraction", shape_val = 'Nutrient',
          mypch = 21,size = 20,alpha=1,comp_a = "CAP1",comp_b = "CAP2",
          col_shape_background = "black",
          size_legend_text = size_legend_text,
          size_axis_title = size_axis_title,size_axis_line = 2,
          size_title_text = size_legend_title,
          font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  #theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

######text
mtext <- 'R\u00B2=0.24
p-value=0.0001'

#######plot with the confifence interval
p_all <- ggplot(data = mcap$Map_cap,mapping = aes(CAP1,CAP2, color = Fraction,
                                                  fill = Fraction)) +
  geom_vline(xintercept = 0,linetype = "dashed",size = 2 , color = "#D9D9D9") +
  geom_hline(yintercept = 0,linetype = "dashed",size = 2 , color = "#D9D9D9") +
  #confidence interval = 95%
  #stat_ellipse(aes(group = Fraction, fill = Fraction),
               #geom = "polygon",alpha = 0.2) +
  geom_point(aes(fill = Fraction,shape = Nutrient),size = 12) +
  xlim(-1.5, 1.5)+
  scale_shape_manual(values = 21:24) + scale_fill_manual(values = paleta_fra)+
  scale_color_manual(values = paleta_fra)+
  theme(panel.grid.major.y = element_blank(),
                           panel.grid.minor.y = element_blank()) +
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank()) +
  xlab(label = "CAP1(79.94%)") + ylab(label = "CAP2(16.97%)") + 
  annotate(geom = "text",x = -1.3,y = 2.3,label = mtext,size = 10)+
  clean

p_all

########save the plot
oh.save.pdf(p = p_all,outname = 
              "figS3i-CAPs_maize_syncom_fractions.pdf",
            outdir = "../figures/",width = 20,height = 20)