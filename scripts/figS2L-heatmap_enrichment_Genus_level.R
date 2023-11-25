######load packages
library(ohchibi)
library(ggtree)
library(paletteer)
library(scales)

######set seed
set.seed(12356)

######set work directory
setwd("C:/Users/Castrillo Lab/OneDrive - The University of Nottingham/Publication 2019/Paper preparation/manuscripts/scripts")

######helper functions
source("0-Clean_up_plots.R")
source('0_11-theme_ohchibi_2.R')

##### Open the R dataset and the deseq data frame
df <- read.csv(file = "../cleandata/Enrichment_genera_level_exp_new.csv", 
               header = T, sep = "\t")

######Define some column as factor
df$contrast <- df$contrast %>% as.factor
df$Genus <- df$Genus %>% as.factor

######Add the significance column
pthres <- 0.05
df$Significance <- rep("NS", nrow(df))
df$Significance[which(df$padj < pthres)] <- "q < 0.05"
df$Significance <- df$Significance %>% factor

######Subset ASVs that are significant
uids <- df %>% subset(Significance == "q < 0.05") %$% Genus %>% unique

######create a log2foldchange matrix with significant ASV
res <- which(df$Genus %in% uids) %>%
  df[.,] %>% droplevels

######Change Log2FoldChange range
res$log2FoldChange[which(res$log2FoldChange > 5)] <- 5
res$log2FoldChange[which(res$log2FoldChange < -5)] <- -5

#####combine the soi and contrast
res$contrast_over <- paste(res$contrast, res$Soil, sep='_')

#####drop NA
res_final <- na.omit(res)

#######remove the uncultured from the dataframe
res_final <- res_final %>%
  subset(Genus != 'Bacillales_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Bacillaceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'BIrii41_ge') %>% droplevels
res_final <- res_final %>%
  subset(Genus != '[Aquaspirillum]_arcticum_group') %>% droplevels
res_final <- res_final %>%
  subset(Genus != '[Polyangium]_brachysporum_group') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Acetobacteraceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Family != 'Actinobacteria_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Alphaproteobacteria_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Archangiaceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Azospirillaceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Burkholderiaceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Candidatus_Levybacteria_ge') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Caulobacteraceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Chitinophagaceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Chloroflexaceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Cytophagales_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Deltaproteobacteria_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Devosiaceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Enterobacteriaceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'env.OPS_17_ge') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Gaiellales_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Gemmatimonadaceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Geodermatophilaceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Intrasporangiaceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Micrococcaceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Micromonosporaceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Microscillaceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Myxococcales_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Nocardioidaceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Paenibacillaceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Pedosphaeraceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Pseudonocardiaceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Rhizobiaceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'S0134_terrestrial_group_ge') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Saccharimonadales_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Solibacteraceae_(Subgroup_3)_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Sphingomonadaceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Streptomycetaceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Streptosporangiaceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Thermomonosporaceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Xanthobacteraceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Xanthomonadaceae_unclassified') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'uncultured') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'uncultured_ge') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'WD2101_soil_group_ge') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'possible_genus_04') %>% droplevels
res_final <- res_final %>%
  subset(Genus != 'Acidimicrobiia_unclassified') %>% droplevels

######Average the Contrast
Tab <- acast(data = res_final, 
             formula = Genus ~ contrast_over, 
             value.var = "log2FoldChange",fun.aggregate = mean)

Tab_final <- na.omit(Tab)

#######Cluster and order the family patterns
mclust <- hclust(d = (1-cor(Tab_final %>%t)) %>% as.dist,method = "ward.D2")
order <- mclust$order %>% mclust$labels[.]

mclust_samples <- hclust(d = (1-cor(Tab_final)) %>% as.dist,method = "ward.D2")
order_samples <- mclust_samples$order %>% mclust_samples$labels[.]

######Visualise the dendogram
mclust %>% plot
mclust_samples %>% plot

###### Melt the df_clust structures with the structure used to plot
melted <- Tab_final %>% melt
colnames(melted) <- c("Genus","contrast_over","log2FoldChange")

###Order genus according to the cluster dendogram
melted$Genus <- melted$Genus %>% factor(levels = order)

#####Divide the dendogram in cluster
df_samples <- mclust %>% dendextend::cutree(tree=., k=5) %>%
  data.frame(Genus=names(.), Clustsamples=paste0("Cl", .),
             row.names = NULL) 

#remove first column
df_samples <- df_samples[,-1]

######reorder the df_samples ASV column
df_samples$Clustsamples[which(df_samples$Clustsamples == 'Cl4')] <- 'Cl1_1'
df_samples$Clustsamples[which(df_samples$Clustsamples == 'Cl1')] <- 'Cl2_1'
df_samples$Clustsamples[which(df_samples$Clustsamples == 'Cl5')] <- 'Cl3_1'
df_samples$Clustsamples[which(df_samples$Clustsamples == 'Cl2')] <- 'Cl4_1'
df_samples$Clustsamples[which(df_samples$Clustsamples == 'Cl3')] <- 'Cl5_1'

df_samples$Clustsamples[which(df_samples$Clustsamples == 'Cl5_1')] <- 'Cl5'
df_samples$Clustsamples[which(df_samples$Clustsamples == 'Cl4_1')] <- 'Cl4'
df_samples$Clustsamples[which(df_samples$Clustsamples == 'Cl3_1')] <- 'Cl3'
df_samples$Clustsamples[which(df_samples$Clustsamples == 'Cl2_1')] <- 'Cl2'
df_samples$Clustsamples[which(df_samples$Clustsamples == 'Cl1_1')] <- 'Cl1'

######merge the melted and df_samples
melted_sum <- merge(melted, df_samples, by='Genus')

######check the change in the log2foldchange
melted_sum$log2FoldChange %>% sort %>% plot

#####Merge melted object with the metadata
melted_1 <- merge(melted_sum, df, by = c("Genus", "log2FoldChange"))

######correct the soil names
melted_1$Soil <- melted_1$Soil %>%
  gsub(pattern = 'SaoDomingos', replacement = 'São\nDomingos')

melted_1$Soil <- melted_1$Soil %>%
  gsub(pattern = 'SuttonBonington', replacement = 'Sutton\nBonington')

#Order the soils
melted_1$Soil <- melted_1$Soil %>%
  factor(levels = c("São\nDomingos", "Tarrafal", "Sutton\nBonington"))

######reorder the contrast
melted_1$contrast <- melted_1$contrast %>%
  factor(levels = c('Leaf_vs_Soil', 'Leaf_vs_Root',
                    'Root_vs_Soil'))

######correct the contrast names
melted_1$contrast <- melted_1$contrast %>%
  gsub(pattern = 'Leaf_vs_Soil', 
       replacement = 'Leaf vs Soil')

melted_1$contrast <- melted_1$contrast %>%
  gsub(pattern = 'Leaf_vs_Root', 
       replacement = 'Leaf vs Root')

melted_1$contrast <- melted_1$contrast %>%
  gsub(pattern = 'Root_vs_Soil', 
       replacement = 'Root vs Soil')

######reorder the contrast
melted_1$contrast <- melted_1$contrast %>%
  factor(levels = c('Leaf vs Soil',
                    'Leaf vs Root',
                    'Root vs Soil'))

######create the heatmap
heatmap_soil <- ggplot(data = melted_1, mapping = aes(x=Soil,y=Genus))+
  geom_raster(aes(fill = log2FoldChange), stat = "identity")+
  geom_tile(aes(color = Significance),fill = '#00000000', 
            size = 1,width = 0.85,height = 0.85) +
  facet_grid(Clustsamples~contrast, scales = "free",space ="free")+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-5,5),oob = squish) +
  scale_color_manual(values = c("#00000000","#414141"))+
  xlab(NULL)+
  ylab(NULL)+
  clean +
  theme_ohchibi( size_axis_title.x = 22,
                 size_axis_title.y = 22,
                 legend_proportion_size = 2,
                 size_title_text = 30,
                 size_legend_text = 20,
                 size_panel_border = 1.5,
                 size_lines_panel = 0,
                 size_axis_text.y = 22,
                 size_axis_text.x = 22) +
  theme(panel.spacing = unit(0.1, 'lines'),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 30),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

heatmap_soil

######save the heatmap
oh.save.pdf(p = heatmap_soil,
            outname = "test_figS2E_model_enrichment_genus_level_heatmap.pdf",
            outdir = "../figures/",width = 20,height = 15)

######Define the color for phyla
names_phyla <-c("Acidobacteria",
                "Actinobacteria",
                "Bacteroidetes",
                "Chloroflexi",
                "Cyanobacteria",
                "Firmicutes",
                "Gemmatimonadetes",
                "Patescibacteria",
                "Proteobacteria",
                "Verrucomicrobia",
                "Other")
paleta <- c(RColorBrewer::brewer.pal(n = 10,name = "Paired"),"grey")
names(paleta) <- names_phyla
mphyla <- names(paleta)

#######transform phylum in factor
melted_1$Phylum <- melted_1$Phylum %>%
  factor()

######change the phyla name to others 
melted_1$Phylum <- melted_1$Phylum %>% as.character
melted_1$Phylum[which(melted_1$Phylum == 'Fibrobacteres')] <- 'Other'
melted_1$Phylum[which(melted_1$Phylum == 'Nitrospirae')] <- 'Other'
melted_1$Phylum[which(melted_1$Phylum == 'Planctomycetes')] <- 'Other'
melted_1$Phylum[which(melted_1$Phylum == 'Rokubacteria')] <- 'Other'

######order phylum in melted_1
melted_1$Phylum <- melted_1$Phylum %>%
  factor(levels = c("Acidobacteria",
                    "Actinobacteria",
                    "Bacteroidetes",
                    "Chloroflexi",
                    "Cyanobacteria",
                    "Firmicutes",
                    "Gemmatimonadetes",
                    "Patescibacteria",
                    "Proteobacteria",
                    "Verrucomicrobia",
                    "Other"))

######arrange the plot
tree <- mclust %>% as.phylo

p_tiplab <- ggtree(tree,ladderize = FALSE,size = 1.2) +
  geom_tiplab(size=2)+
  #The expand parth is fundamental to make the tree fit the composition with the other plots
  scale_y_continuous(expand = c(0.002,0))

p_tiplab

p <- ggtree(tree,ladderize = FALSE,size = 1.2) +
  #The expand parth is fundamental to make the tree fit the composition with the other plots
  scale_y_continuous(expand = c(0.0045,0.0045))
p_tree_legend <- p %<+% melted_1 + 
  ######colored the tips
  geom_tippoint(aes(fill=Phylum, color=Phylum),
                size=4, shape=21)+
  scale_fill_manual(values = paleta)+
  scale_color_manual(values = paleta)+
  theme(legend.position = 'left',
        legend.text = element_text(size = 20))

p_tree <- p %<+% melted_1 + 
  ######colored the tips
  geom_tippoint(aes(fill=Phylum, color=Phylum),
                size=4, shape=21)+
  scale_fill_manual(values = paleta)+
  scale_color_manual(values = paleta)+
  theme(legend.position = 'none')

######tree composition
composition <- egg::ggarrange(p_tree_legend,heatmap_soil,
                              nrow=1, ncol=2, widths = c(0.5,1))

composition_1 <- egg::ggarrange(p_tree,heatmap_soil,
                                nrow=1, ncol=2, widths = c( 0.09,0.5))

######save the plot as pdf
oh.save.pdf(p = composition,
            outname = "test_figS2E_model_enrichment_genera_level.pdf",
            outdir = "../figures/",width = 40,height = 35)

oh.save.pdf(p = composition_1,outname = "figS2E_model_enrichment_genera_level_nolegend.pdf",
            outdir = "../figures/",width = 20,height = 35)
