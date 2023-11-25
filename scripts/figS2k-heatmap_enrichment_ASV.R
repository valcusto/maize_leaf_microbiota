######load packages
library(ohchibi)
library(ggtree)

######set seed
set.seed(12356)

######set work directory
setwd("C:/Users/Castrillo Lab/OneDrive - The University of Nottingham/Publication 2019/Paper preparation/manuscripts/scripts")

##### Open the R dataset and the deseq data frame
df <- read.csv(file = "../cleandata/Enrichment_ASV_level_exp.csv", 
               header = T, sep = "\t")

######Define some column as factor
df$contrast <- df$contrast %>% as.factor
df$Overall <- df$Overall %>% as.factor
df$Root <- df$Root %>% as.factor
df$Kingdom <- df$Kingdom %>% as.factor
df$Phylum <- df$Phylum %>% as.factor
df$Class <- df$Class %>% as.factor
df$Order <- df$Order %>% as.factor
df$Family <- df$Family %>% as.factor
df$Genus <- df$Genus %>% as.factor

######Add the significance column
pthres <- 0.05
df$Significance <- rep("NS", nrow(df))
df$Significance[which(df$padj < pthres)] <- "q < 0.05"
df$Significance <- df$Significance %>% factor

######Change Log2FoldChange range
df$log2FoldChange[which(df$log2FoldChange > 5)] <- 5
df$log2FoldChange[which(df$log2FoldChange < -5)] <- -5

######Subset ASVs that are significant
uids <- df %>% subset(Significance == "q < 0.05") %$% ASV %>% unique

######create a log2foldchange matrix with significant ASV
res <- which(df$ASV %in% uids) %>%
  df[.,] %>% droplevels

#####combine the soi and contrast
res$contrast_over <- paste(res$contrast, res$Soil, sep='_')

#####drop NA
res_final <- na.omit(res)

######Average the Contrast
Tab <- acast(data = res_final, 
             formula = ASV ~ contrast_over, 
             value.var = "log2FoldChange")

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
colnames(melted) <- c("ASV","contrast_over","log2FoldChange")

###Order family according to the cluster dendogram
melted$ASV <- melted$ASV %>% factor(levels = order)

#####Divide the dendogram in cluster
df_samples <- mclust %>% dendextend::cutree(tree=., k=5) %>%
  data.frame(ASV=names(.), Clustsamples=paste0("Cl", .),
             row.names = NULL) 

#remove first column
df_samples <- df_samples[,-1]

#####order the df_samples
#df_samples <- df_samples %>%
  #slice(match(order, ASV) %>% rev)

######reorder the df_samples ASV column
df_samples$Clustsamples[which(df_samples$Clustsamples == 'Cl2')] <- 'Cl1_1'
df_samples$Clustsamples[which(df_samples$Clustsamples == 'Cl3')] <- 'Cl2_1'
df_samples$Clustsamples[which(df_samples$Clustsamples == 'Cl1')] <- 'Cl3_1'
df_samples$Clustsamples[which(df_samples$Clustsamples == 'Cl5')] <- 'Cl4_1'
df_samples$Clustsamples[which(df_samples$Clustsamples == 'Cl4')] <- 'Cl5_1'

df_samples$Clustsamples[which(df_samples$Clustsamples == 'Cl1_1')] <- 'Cl1'
df_samples$Clustsamples[which(df_samples$Clustsamples == 'Cl2_1')] <- 'Cl2'
df_samples$Clustsamples[which(df_samples$Clustsamples == 'Cl3_1')] <- 'Cl3'
df_samples$Clustsamples[which(df_samples$Clustsamples == 'Cl4_1')] <- 'Cl4'
df_samples$Clustsamples[which(df_samples$Clustsamples == 'Cl5_1')] <- 'Cl5'

######merge the melted and df_samples
melted_sum <- merge(melted, df_samples, by='ASV')

######check the change in the log2foldchange
melted_sum$log2FoldChange %>% sort %>% plot

#####Merge melted object with the metadata
melted_1 <- merge(melted_sum, df, by = c("ASV", "log2FoldChange"))

######correct the soil names
melted_1$Soil <- melted_1$Soil %>%
  gsub(pattern = 'SaoDomingos', replacement = 'S?o\nDomingos')

melted_1$Soil <- melted_1$Soil %>%
  gsub(pattern = 'SuttonBonington', replacement = 'Sutton\nBonington')

#Order the soils
melted_1$Soil <- melted_1$Soil %>%
  factor(levels = c("S?o\nDomingos", "Tarrafal", "Sutton\nBonington"))

######correct the contrast names
melted_1$contrast <- melted_1$contrast %>%
  gsub(pattern = 'Leaf_vs_Root', replacement = 'Leaf vs\nRoot')
melted_1$contrast <- melted_1$contrast %>%
  gsub(pattern = 'Leaf_vs_Soil', replacement = 'Leaf vs\nSoil')
melted_1$contrast <- melted_1$contrast %>%
  gsub(pattern = 'Root_vs_Soil', replacement = 'Root vs\nSoil')

######create the heatmap
heatmap_soil <- ggplot(data = melted_1, mapping = aes(x=Soil,y=ASV))+
  geom_raster(aes(fill = log2FoldChange), stat = "identity")+
  geom_tile(aes(color = Significance),fill = '#00000000', size = 1,width = 0.85,height = 0.85) +
  facet_grid(Clustsamples~contrast, scales = "free", 
             space ="free")+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",limits = c(-5,5),oob = squish) +
  scale_color_manual(values = c("#00000000","#414141"))+
  xlab(NULL)+
  ylab(NULL)+
  theme_ohchibi(
    size_axis_title.x = 22,
    size_axis_title.y = 22,
    legend_proportion_size = 2,
    size_title_text = 30,
    size_legend_text = 20,
    size_panel_border = 1.5,
    size_lines_panel = 0) +
  theme(panel.spacing = unit(0.1, 'lines'),
        strip.background = element_blank())

heatmap_soil

######Define the color for phyla
mphyla <- palettesPM::pm.names.phyla()
paleta <- palettesPM::pm.colors.phyla()

######arrange the plot
tree <- mclust %>% as.phylo
p <- ggtree(tree,ladderize = FALSE,size = 1.2) +
  #The expand parth is fundamental to make the tree fit the composition with the other plots
  scale_y_continuous(expand = c(0.009,0.009))
p_tree_legend <- p %<+% melted_1 + 
  ######colored the tips
  geom_tippoint(aes(fill=Phylum, color=Phylum),
                size=4, shape=21)+
  scale_fill_manual(values = paleta)+
  scale_color_manual(values = paleta)+
  theme(legend.position = 'right',
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
                              nrow=1, ncol=2, widths = c( 0.09,0.5))

composition_1 <- egg::ggarrange(p_tree,heatmap_soil,
                              nrow=1, ncol=2, widths = c( 0.09,0.5))

######save the plot as pdf
oh.save.pdf(p = composition,outname = "figureS2E_model_enrichment_asv_level.pdf",
            outdir = "../figures/",width = 20,height = 35)

oh.save.pdf(p = composition_1,outname = "figureS2E_model_enrichment_asv_level_nolegend.pdf",
            outdir = "../figures/",width = 20,height = 35)