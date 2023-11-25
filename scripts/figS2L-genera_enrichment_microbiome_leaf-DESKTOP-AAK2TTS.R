######load packages
library(ohchibi)
library(DESeq2)

######set seed
set.seed(12356)

######set work directory
setwd("C:/Users/Castrillo Lab/OneDrive - The University of Nottingham/Publication 2019/Paper preparation/manuscripts/scripts")

######helper functions
source('0_4-subset_Dataset.R')
source('0_12-measurable_taxa.R')
source('0_8-collapse_by_taxonomy.R')

######helper functions
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

######upload dataset
dat <- readRDS('../cleandata/dat_asvs_maize_leaves_rep2.RDS')

######select rarefied dataset
dat_rar <- dat$Rarefied

######Select the data for analysis
dat_sub <- dat_rar %>%
  subset.Dataset(Treatment == '50FC', drop = T, clean = T)

######create the fractions column
dat_sub$Map$fractions <- 'Leaf'
dat_sub$Map$fractions[which(dat_sub$Map$NumLeaf == 'Soil')] <- 'Soil'
dat_sub$Map$fractions[which(dat_sub$Map$NumLeaf == 'Root')] <- 'Root'
dat_sub$Map$fractions <- dat_sub$Map$fractions %>%
  factor

#Make soil the reference level so all comparisons are made against it
dat_sub$Map$fractions <- dat_sub$Map$fractions %>% relevel(ref = "Soil")

######Select ASVs that appears at least in 3 samples and with 10 reads min
total <- dat_sub$Tab %>% sum
dat_filter <- measurable_taxa.Dataset(Dat = dat_sub, min_samples_otu = 3,
                                      min_reads_otu = 10)

#####how many ASVs we retained after the filtering
#calculate the percentage prevalence
(dat_filter$Tab %>% sum)/total
###we retained 57% of reads (prevalence)

######rename the dat file
dat_model <- dat_sub

#####Collapse the ASV data at the family levels
dat_fam <- dat_model %>% collapse_by_taxonomy.Dataset(Dat = ., level = 7)

######subset the soils
dat_sd <- dat_fam %>%
  subset.Dataset(Soil == 'SaoDomingos', drop = T, clean = T)

######Model
#DESeq
dds <- DESeqDataSetFromMatrix(countData = dat_sd$Tab, colData = dat_sd$Map,
                              design=~fractions)
geoMeans <- apply(counts(dds), 1, gm_mean)
dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
dds <- DESeq(dds, fitType="local")

######check the pairwise in the analysis
resultsNames(dds)

######recover the results
res <- results(dds) %>% data.frame

######create the datasets for each comparison
res_a <- results(dds, contrast=c('fractions','Leaf', 'Soil')) %>%
  data.frame
res_a$contrast <- rep("Leaf_vs_Soil", nrow(res_a))
res_a$Overall <- rep("Leaf vs Soil")
res_a$Soil <- 'SaoDomingos'
res_a$taxonomy <- rownames(res_a)
res_leaf_soil <- res_a

res_b <- results(dds, contrast=c('fractions','Leaf', 'Root')) %>%
  data.frame
res_b$contrast <- rep("Leaf_vs_Root", nrow(res_b))
res_b$Overall <- rep("Leaf vs Root")
res_b$Soil <- 'SaoDomingos'
res_b$taxonomy <- rownames(res_b)
res_leaf_root <- res_b

res_c <- results(dds, contrast=c('fractions','Root','Soil')) %>%
  data.frame
res_c$contrast <- rep("Root_vs_Soil", nrow(res_c))
res_c$Overall <- rep("Root vs Soil")
res_c$Soil <- 'SaoDomingos'
res_c$taxonomy <- rownames(res_c)
res_root_soil <- res_c

######combine the three dataset
Res_SD <- rbind(res_leaf_root, res_leaf_soil, res_root_soil)

######Prepare the tax file
df_tax_sd <- Res_SD$taxonomy %>% as.character %>%
  strsplit(split  = "\\;") %>%
  unlist %>%  gsub(pattern = "[a-z]__",replacement = "",) %>%
  gsub(pattern = " ",replacement = "") %>%
  matrix(data = .,ncol = 7,byrow = T) %>% as.data.frame
colnames(df_tax_sd) <-c("Root","Kingdom","Phylum","Class",
                        "Order", "Family", "Genus")
df_tax_sd <- cbind(Res_SD$taxonomy,df_tax_sd)
colnames(df_tax_sd)[1] <- "taxonomy"
df_tax_sd$Root <- df_tax_sd$Root %>% factor
df_tax_sd$Kingdom <- df_tax_sd$Kingdom %>% factor
df_tax_sd$Phylum <- df_tax_sd$Phylum %>% factor
df_tax_sd$Class <- df_tax_sd$Class %>% factor
df_tax_sd$Order <- df_tax_sd$Order %>% factor
df_tax_sd$Family <- df_tax_sd$Family %>% factor
df_tax_sd$Genus <- df_tax_sd$Genus %>% factor

#######combine the datasets
Res_SD <- cbind(Res_SD, df_tax_sd)

######subset the soil Tarrafal
dat_t <- dat_fam %>%
  subset.Dataset(Soil == 'Tarrafal', drop = T, clean = T)

######Model
#DESeq
dds <- DESeqDataSetFromMatrix(countData = dat_t$Tab, colData = dat_t$Map,
                              design=~fractions)
geoMeans <- apply(counts(dds), 1, gm_mean)
dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
dds <- DESeq(dds, fitType="local")

######check the pairwise in the analysis
resultsNames(dds)

######recover the results
res <- results(dds) %>% data.frame

######create the datasets for each comparison
res_a <- results(dds, contrast=c('fractions','Leaf', 'Soil')) %>%
  data.frame
res_a$contrast <- rep("Leaf_vs_Soil", nrow(res_a))
res_a$Overall <- rep("Leaf vs Soil")
res_a$Soil <- 'Tarrafal'
res_a$taxonomy <- rownames(res_a)
res_leaf_soil <- res_a

res_b <- results(dds, contrast=c('fractions','Leaf', 'Root')) %>%
  data.frame
res_b$contrast <- rep("Leaf_vs_Root", nrow(res_b))
res_b$Overall <- rep("Leaf vs Root")
res_b$Soil <- 'Tarrafal'
res_b$taxonomy <- rownames(res_b)
res_leaf_root <- res_b

res_c <- results(dds, contrast=c('fractions','Root','Soil')) %>%
  data.frame
res_c$contrast <- rep("Root_vs_Soil", nrow(res_c))
res_c$Overall <- rep("Root vs Soil")
res_c$Soil <- 'Tarrafal'
res_c$taxonomy <- rownames(res_c)
res_root_soil <- res_c

######combine the three dataset
Res_T <- rbind(res_leaf_root, res_leaf_soil, res_root_soil)

######Prepare the tax file
df_tax_t <- Res_T$taxonomy %>% as.character %>%
  strsplit(split  = "\\;") %>%
  unlist %>%  gsub(pattern = "[a-z]__",replacement = "",) %>%
  gsub(pattern = " ",replacement = "") %>%
  matrix(data = .,ncol = 7,byrow = T) %>% as.data.frame
colnames(df_tax_t) <-c("Root","Kingdom","Phylum","Class",
                       "Order", "Family", "Genus")
df_tax_t <- cbind(Res_T$taxonomy,df_tax_t)
colnames(df_tax_t)[1] <- "taxonomy"
df_tax_t$Root <- df_tax_t$Root %>% factor
df_tax_t$Kingdom <- df_tax_t$Kingdom %>% factor
df_tax_t$Phylum <- df_tax_t$Phylum %>% factor
df_tax_t$Class <- df_tax_t$Class %>% factor
df_tax_t$Order <- df_tax_t$Order %>% factor
df_tax_t$Family <- df_tax_t$Family %>% factor
df_tax_t$Genus <- df_tax_t$Genus %>% factor

#######combine the datasets
Res_T <- cbind(Res_T, df_tax_t)

######subset the soil SuttonBonington
dat_sb <- dat_fam %>%
  subset.Dataset(Soil == 'SuttonBonington', drop = T, clean = T)

######Model
#DESeq
dds <- DESeqDataSetFromMatrix(countData = dat_sb$Tab, colData = dat_sb$Map,
                              design=~fractions)
geoMeans <- apply(counts(dds), 1, gm_mean)
dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
dds <- DESeq(dds, fitType="local")

######check the pairwise in the analysis
resultsNames(dds)

######recover the results
res <- results(dds) %>% data.frame

######create the datasets for each comparison
res_a <- results(dds, contrast=c('fractions','Leaf', 'Soil')) %>%
  data.frame
res_a$contrast <- rep("Leaf_vs_Soil", nrow(res_a))
res_a$Overall <- rep("Leaf vs Soil")
res_a$Soil <- 'SuttonBonington'
res_a$taxonomy <- rownames(res_a)
res_leaf_soil <- res_a

res_b <- results(dds, contrast=c('fractions','Leaf', 'Root')) %>%
  data.frame
res_b$contrast <- rep("Leaf_vs_Root", nrow(res_b))
res_b$Overall <- rep("Leaf vs Root")
res_b$Soil <- 'SuttonBonington'
res_b$taxonomy <- rownames(res_b)
res_leaf_root <- res_b

res_c <- results(dds, contrast=c('fractions','Root','Soil')) %>%
  data.frame
res_c$contrast <- rep("Root_vs_Soil", nrow(res_c))
res_c$Overall <- rep("Root vs Soil")
res_c$Soil <- 'SuttonBonington'
res_c$taxonomy <- rownames(res_c)
res_root_soil <- res_c

######combine the three dataset
Res_sb <- rbind(res_leaf_root, res_leaf_soil, res_root_soil)

######Prepare the tax file
df_tax_sb <- Res_sb$taxonomy %>% as.character %>%
  strsplit(split  = "\\;") %>%
  unlist %>%  gsub(pattern = "[a-z]__",replacement = "",) %>%
  gsub(pattern = " ",replacement = "") %>%
  matrix(data = .,ncol = 7,byrow = T) %>% as.data.frame
colnames(df_tax_sb) <-c("Root","Kingdom","Phylum","Class",
                        "Order", "Family", "Genus")
df_tax_sb <- cbind(Res_sb$taxonomy,df_tax_sb)
colnames(df_tax_sb)[1] <- "taxonomy"
df_tax_sb$Root <- df_tax_sb$Root %>% factor
df_tax_sb$Kingdom <- df_tax_sb$Kingdom %>% factor
df_tax_sb$Phylum <- df_tax_sb$Phylum %>% factor
df_tax_sb$Class <- df_tax_sb$Class %>% factor
df_tax_sb$Order <- df_tax_sb$Order %>% factor
df_tax_sb$Family <- df_tax_sb$Family %>% factor
df_tax_sb$Genus <- df_tax_sb$Genus %>% factor

#######combine the datasets
Res_sb <- cbind(Res_sb, df_tax_sb)

####Create the data frame with all data
Res <- rbind(Res_SD, Res_T, Res_sb)
rownames(Res) <- NULL

######Save the object as csv file
write.table(x=Res, file="../cleandata/Enrichment_genera_level_exp_new.csv",
            append = F, quote = F, sep="\t", row.names = F, col.names = T)

######count the number of genera per fractions and soil
melted_sd %>%
  group_by(soil) %>% summarise(no_rows=length(soil))

melted_sd %>%
  group_by(genotype, num_leaf, Ion) %>% 
  summarise(no_rows=length(genotype))

melted_sd %>%
  group_by(num_leaf, Ion) %>% 
  summarise(no_rows=length(genotype))