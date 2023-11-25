######load packages
library(ohchibi)
library(DESeq2)

######set seed
set.seed(12356)

######set work directory
setwd("C:/Users/Castrillo Lab/OneDrive - The University of Nottingham/Publication 2019/Paper preparation/manuscripts/scripts")

######helper functions
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

######upload dataset
######upload data
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
(dat_filter$Tab %>% sum)/total
###we retained 57% of reads

######rename the dat file
dat_model <- dat_sub

######subset the soils
dat_sd <- dat_model %>%
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
res_a$ASV <- rownames(res_a)
res_leaf_soil <- res_a

res_b <- results(dds, contrast=c('fractions','Leaf', 'Root')) %>%
  data.frame
res_b$contrast <- rep("Leaf_vs_Root", nrow(res_b))
res_b$Overall <- rep("Leaf vs Root")
res_b$Soil <- 'SaoDomingos'
res_b$ASV <- rownames(res_b)
res_leaf_root <- res_b

res_c <- results(dds, contrast=c('fractions','Root','Soil')) %>%
  data.frame
res_c$contrast <- rep("Root_vs_Soil", nrow(res_c))
res_c$Overall <- rep("Root vs Soil")
res_c$Soil <- 'SaoDomingos'
res_c$ASV <- rownames(res_c)
res_root_soil <- res_c

######combine the three dataset
Res_SD <- rbind(res_leaf_root, res_leaf_soil, res_root_soil)

######Prepare the tax file
mdf_sd <- dat_sd$Tax %>% data.frame
colnames(mdf_sd)[1] <- "ASV"
Res_SD <- merge(Res_SD,mdf_sd, by = "ASV")

######subset the soil Tarrafal
dat_t <- dat_model %>%
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
res_a$ASV <- rownames(res_a)
res_leaf_soil <- res_a

res_b <- results(dds, contrast=c('fractions','Leaf', 'Root')) %>%
  data.frame
res_b$contrast <- rep("Leaf_vs_Root", nrow(res_b))
res_b$Overall <- rep("Leaf vs Root")
res_b$Soil <- 'Tarrafal'
res_b$ASV <- rownames(res_b)
res_leaf_root <- res_b

res_c <- results(dds, contrast=c('fractions','Root','Soil')) %>%
  data.frame
res_c$contrast <- rep("Root_vs_Soil", nrow(res_c))
res_c$Overall <- rep("Root vs Soil")
res_c$Soil <- 'Tarrafal'
res_c$ASV <- rownames(res_c)
res_root_soil <- res_c

######combine the three dataset
Res_T <- rbind(res_leaf_root, res_leaf_soil, res_root_soil)

######Prepare the tax file
mdf_t <- dat_t$Tax %>% data.frame
colnames(mdf_t)[1] <- "ASV"
Res_T <- merge(Res_T,mdf_t, by = "ASV")

######subset the soil SuttonBonington
dat_sb <- dat_model %>%
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
res_a$ASV <- rownames(res_a)
res_leaf_soil <- res_a

res_b <- results(dds, contrast=c('fractions','Leaf', 'Root')) %>%
  data.frame
res_b$contrast <- rep("Leaf_vs_Root", nrow(res_b))
res_b$Overall <- rep("Leaf vs Root")
res_b$Soil <- 'SuttonBonington'
res_b$ASV <- rownames(res_b)
res_leaf_root <- res_b

res_c <- results(dds, contrast=c('fractions','Root','Soil')) %>%
  data.frame
res_c$contrast <- rep("Root_vs_Soil", nrow(res_c))
res_c$Overall <- rep("Root vs Soil")
res_c$Soil <- 'SuttonBonington'
res_c$ASV <- rownames(res_c)
res_root_soil <- res_c

######combine the three dataset
Res_sb <- rbind(res_leaf_root, res_leaf_soil, res_root_soil)

######Prepare the tax file
mdf <- dat_sb$Tax %>% data.frame
colnames(mdf)[1] <- "ASV"
Res_sb <- merge(Res_sb,mdf, by = "ASV")

####Create the data frame with all data
Res <- rbind(Res_SD, Res_T, Res_sb)
rownames(Res) <- NULL

######taxonomy name
df_tax <- Res$Taxonomy %>% as.character %>%
  strsplit(split  = "\\;") %>%
  unlist %>%  gsub(pattern = "[a-z]__",replacement = "",) %>%
  gsub(pattern = " ",replacement = "") %>%
  matrix(data = .,ncol = 7,byrow = T) %>% as.data.frame
colnames(df_tax) <-c("Root","Kingdom","Phylum",
                     "Class","Order","Family","Genus")
df_tax <- cbind(Res$ASV,df_tax)
colnames(df_tax)[1] <- "ASV"
df_tax$ASV <- df_tax$ASV %>% factor
df_tax$Root <- df_tax$Root %>% factor
df_tax$Kingdom <- df_tax$Kingdom %>% factor
df_tax$Phylum <- df_tax$Phylum %>% factor
df_tax$Class <- df_tax$Class %>% factor
df_tax$Order <- df_tax$Order %>% factor
df_tax$Family <- df_tax$Family %>% factor
df_tax$Genus <- df_tax$Genus %>% factor
#rownames(df_tax) <- df_tax$ASV
#Dat_raw$Tax <- df_tax

#######combine the datasets
Res_final <- cbind(Res, df_tax)

######Save the object as csv file
write.table(x=Res_final, file="../cleandata/Enrichment_ASV_level_exp.csv",
            append = F, quote = F, sep="\t", row.names = F, col.names = T)
