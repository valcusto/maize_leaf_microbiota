######load packages
library(ohchibi)
library(tidyverse)

######set seed
set.seed(130816)

######set work directory
setwd("C:/Users/Castrillo Lab/OneDrive - The University of Nottingham/Publication 2019/Paper preparation/manuscripts/scripts")

######helper functions
source('0-Clean_up_plots.R')
source('0_4-subset_Dataset.R')

######upload data
dat <- readRDS('../cleandata/dat_asvs_maize_leaves_rep2.RDS')


#################################STEP01: DETERMINE THE PREVALENCE OF EACH ASV
#Rarefied version
Dat_raw <- dat$RawCounts

######################select only 50FC and SB##############################
Dat_raw_sb <- Dat_raw %>% subset.Dataset((Treatment == '50FC') &
                                           (NumLeaf != 'Soil') &
                                           (NumLeaf != 'Root') &
                                           (NumLeaf != 'L6') &
                                           (Soil == 'SuttonBonington'), 
                                         drop = T, clean = T)

####Define the RA
tab=Dat_raw_sb$Tab
#######Calculate the relative abundance
###understand sweep function
#sweep(x, MARGIN, STATS, FUN="-", check.margin=T, ...)
# x is the data
#Whether you operate by row (2) or column(1) is defined by MARGIN
# STATS refers to the summary statistics which you wish to 'sweep out'
# FUN is the function used to carry out the sweep, "-" is the default
tab_ra=sweep(tab,2,colSums(tab),"/")
##sweep is diving each row value by the column sums.In another word,
#for each row evaluated STATS. This value represents the 
#relative abundance
##Relative abundance is the percent composition of an ASV in 
#the samples relative to the total number of ASV in the samples
tab_c=tab_ra>0.001

######Melt the RA table
melted_prev <- tab_c %>% melt
colnames(melted_prev) <- c("ASV", "DADA2_Header","Present")

######remove all the present = False
melted_prev_sb <- melted_prev %>%
  subset(Present == TRUE) %>% droplevels

melted_prev <- merge(melted_prev_sb,Dat_raw_sb$Map, by = "DADA2_Header")

####Splits the data into subsets and computes statistics for each, and returns
##Count the number of Present ASVs per soil and per leaves
df_sum <- aggregate(Present ~ ASV+Soil+NumLeaf,
                    melted_prev,sum)

Mapa <- Dat_raw_sb$Map
Mapa$Count <- 1

###Calculate the number of samples/compartment
df_tot <- aggregate(Count~Soil+NumLeaf,Mapa,sum)

df_sum <- merge(df_sum,df_tot, by = c("Soil","NumLeaf"))

df_sum$Prevalence <- df_sum$Present/df_sum$Count

######create a column group
df_sum$group <- paste(df_sum$Soil, df_sum$NumLeaf, df_sum$ASV, sep='_')

######select all the ASVs with a prevalence > 30%
df_sum_03 <- df_sum %>%
  subset(Prevalence > 0.3) %>% droplevels

################STEP2:SELECT SAMPLES WITH READS > 1,000 AND PRESENT ASV PREVALENCE>30%
##########open the Rarefied version for Sutton Bonington
Dat_rar <- dat$Rarefied

######select only 50FC
Dat_rar_sb <- Dat_rar %>% subset.Dataset((Treatment == '50FC') &
                                           (NumLeaf != 'Soil') &
                                           (NumLeaf != 'Root') &
                                           (NumLeaf != 'L6') &
                                           (Soil == 'SuttonBonington'), 
                                         drop = T, clean = T)

#####Data exploration
Tab_sb <- Dat_rar_sb$Tab %>% t

######select ASVs with prevalence > 30%
usable_samples_sb <- df_sum_03$ASV %>% unique %>% as.character #unique ASV 133
Tab_sb <- Tab_sb[,colnames(Tab_sb) %in% usable_samples_sb]

##################################STEP03: RESCALE THE EACH VALUE IN THE DATASET
#####Transform the data in z_score
#A z-score is a measure of how many standard deviations 
#(how far my observation is from the mean) below or above
#the populatio mean a raw score is.
Tab_z_sb <- Tab_sb %>% scale

#####Starts to build the PCA
colMeans(Tab_z_sb) %>% sort %>% plot

#####create the Map datasets
Map_sb <- Dat_rar_sb$Map

######create a group column with DADA2_Header and Leaves
Map_sb$group <- paste(Map_sb$Soil, Map_sb$NumLeaf, sep='_')

######select the column to add the group column with the Tab
group_sb <- Map_sb[,c(7,10)]

######merge Tab with group
Tab1 <- data.frame(Tab_z_sb)

#######bind the Tab1 with group
Tab_z1 <- cbind(group_sb, Tab1)

######change the rowname for Tab_z1
rownames(Tab_z1) <- NULL

######remove the first 2 columns from Tab_z1
Tab_z1 <- Tab_z1[,-c(1)]

######melt the dataset. The RA is a zscore
melted <- Tab_z1 %>% melt
colnames(melted) <- c('group','ASV','RA')

######average the RA of ASVs
Tab_sum <- acast(data = melted, formula = ASV~group,
                 fun.aggregate = mean, value.var = 'RA')

######cluster ASVs
mclust_asv <- hclust(d=as.dist(1-cor(Tab_sum %>% t)), method = 'ward.D2')
mclust_samples <- hclust(d=as.dist(1-cor(Tab_sum)), method = 'ward.D2')

######plot cluster
plot(mclust_asv)
plot(mclust_samples)

#####Divide the dendogram in cluster
df <- mclust_asv %>% dendextend::cutree(tree=., k=6) %>%
  data.frame(ASV=names(.), ClustASV=paste0("Cl", .),
             row.names = NULL)

#remove first column
df <- df[,-1]

######order the ASVs
order <- mclust_asv$order %>% mclust_asv$labels[.]

######melt Tab_sum object
melted2 <- Tab_sum %>% melt
colnames(melted2) <- c('ASV', 'group', 'RA_zscore')

#####create a column with soil and NumLeaf 
melted2$group <- melted2$group %>% as.character
melted2$Soil <- sapply(strsplit(melted2$group, '_'), `[`, 1)
melted2$NumLeaf <- sapply(strsplit(melted2$group, '_'), `[`, 2)

#####change the zscore values to <0 when the value is negative
melted2$RA_zscore[which(melted2$RA_zscore < 0)] <- 0

######create a new dataset
merged3_sb <- merge(melted2, df, by='ASV')

#######################STEP04: SELECT ASV WITH ZSCORE>0.5
######select the ASV with zscore higher than 0.3
merged4_sb <- merged3_sb %>%
  subset(RA_zscore > 0.5) %>% droplevels

#################################STEP05: ADD THE TAXONOMY TO THE DATA
######merge merged4_sb with Dat_sb$Tax
df_tax_sb <- Dat_rar_sb$Tax

######rename the first column as ASV
colnames(df_tax_sb)[1] <- 'ASV'

######merge the two objects
merged4_sb <- merge(merged4_sb, df_tax_sb, by='ASV')

######taxonomy name
df_tax_genus_sb <- merged4_sb$Taxonomy %>% as.character %>%
  strsplit(split  = "\\;") %>%
  unlist %>%  gsub(pattern = "[a-z]__",replacement = "",) %>%
  gsub(pattern = " ",replacement = "") %>%
  matrix(data = .,ncol = 7,byrow = T) %>% as.data.frame
colnames(df_tax_genus_sb) <-c("Root","Kingdom","Phylum","Class",
                     "Order", "Family", "Genus")
df_tax_genus_sb <- cbind(merged4_sb$Taxonomy,df_tax_genus_sb)
colnames(df_tax_genus_sb)[1] <- "Taxonomy"
df_tax_genus_sb$Root <- df_tax_genus_sb$Root %>% factor
df_tax_genus_sb$Kingdom <- df_tax_genus_sb$Kingdom %>% factor
df_tax_genus_sb$Phylum <- df_tax_genus_sb$Phylum %>% factor
df_tax_genus_sb$Class <- df_tax_genus_sb$Class %>% factor
df_tax_genus_sb$Order <- df_tax_genus_sb$Order %>% factor
df_tax_genus_sb$Family <- df_tax_genus_sb$Family %>% factor
df_tax_genus_sb$Genus <- df_tax_genus_sb$Genus %>% factor

#######combine the datasets
Res_final_sb <- cbind(merged4_sb, df_tax_genus_sb)

######################select only 50FC and SD##############################
Dat_raw_sd <- Dat_raw %>% subset.Dataset((Treatment == '50FC') &
                                           (NumLeaf != 'Soil') &
                                           (NumLeaf != 'Root') &
                                           (NumLeaf != 'L6') &
                                           (Soil == 'SaoDomingos'), 
                                         drop = T, clean = T)

####Define the RA
tab=Dat_raw_sd$Tab
#######Calculate the relative abundance
###understand sweep function
#sweep(x, MARGIN, STATS, FUN="-", check.margin=T, ...)
# x is the data
#Whether you operate by row (2) or column(1) is defined by MARGIN
# STATS refers to the summary statistics which you wish to 'sweep out'
# FUN is the function used to carry out the sweep, "-" is the default
tab_ra=sweep(tab,2,colSums(tab),"/")
##sweep is diving each row value by the column sums.In another word,
#for each row evaluated STATS. This value represents the 
#relative abundance
##Relative abundance is the percent composition of an ASV in 
#the samples relative to the total number of ASV in the samples
tab_c=tab_ra>0.001

######Melt the RA table
melted_prev <- tab_c %>% melt
colnames(melted_prev) <- c("ASV", "DADA2_Header","Present")

######remove all the present = False
melted_prev_sd <- melted_prev %>%
  subset(Present == TRUE) %>% droplevels

melted_prev_sd <- merge(melted_prev_sd,Dat_raw_sd$Map, by = "DADA2_Header")

####Splits the data into subsets and computes statistics for each, and returns
##Count the number of Present ASVs per soil and per leaves
df_sum_sd <- aggregate(Present ~ ASV+Soil+NumLeaf,
                    melted_prev_sd,sum)

Mapa_sd <- Dat_raw_sd$Map
Mapa_sd$Count <- 1

###Calculate the number of samples/compartment
df_tot_sd <- aggregate(Count~Soil+NumLeaf,Mapa_sd,sum)

df_sum_sd <- merge(df_sum_sd,df_tot_sd, by = c("Soil","NumLeaf"))

df_sum_sd$Prevalence <- df_sum_sd$Present/df_sum_sd$Count

######create a column group
df_sum_sd$group <- paste(df_sum_sd$Soil, df_sum_sd$NumLeaf, df_sum_sd$ASV, sep='_')

######select all the ASVs with a prevalence > 0.3
df_sum_03_sd <- df_sum_sd %>%
  subset(Prevalence > 0.3) %>% droplevels

#Rarefied version
Dat_rar <- dat$Rarefied

######select only 50FC
Dat_rar_sd <- Dat_rar %>% subset.Dataset((Treatment == '50FC') &
                                           (NumLeaf != 'Soil') &
                                           (NumLeaf != 'Root') &
                                           (NumLeaf != 'L6') &
                                           (Soil == 'SaoDomingos'), 
                                         drop = T, clean = T)

#####Data exploration
Tab_sd <- Dat_rar_sd$Tab %>% t

######select the prevalence ASV
usable_samples_sd <- df_sum_03_sd$ASV %>% unique %>% as.character #unique ASV 118
Tab_sd <- Tab_sd[,colnames(Tab_sd) %in% usable_samples_sd]

#####Transform the data in z_score
#A z-score is a measure of how many standard deviations 
#(how far my observation is from the mean) below or above
#the populatio mean a raw score is.
Tab_z_sd <- Tab_sd %>% scale

#####Starts to build the PCA
colMeans(Tab_z_sd) %>% sort %>% plot

#####create the Map datasets
Map_sd <- Dat_rar_sd$Map

######create a group column with DADA2_Header and Leaves
Map_sd$group <- paste(Map_sd$Soil, Map_sd$NumLeaf, sep='_')

######select the column to add the group column with the Tab
group_sd <- Map_sd[,c(7,10)]

######merge Tab with group
Tab1_sd <- data.frame(Tab_z_sd)

#######bind the Tab1 with group
Tab_z1_sd <- cbind(group_sd, Tab1_sd)

######change the rowname for Tab_z1
rownames(Tab_z1_sd) <- NULL

######remove the first 2 columns from Tab_z1
Tab_z1_sd <- Tab_z1_sd[,-c(1)]

######melt the dataset. The RA is a zscore
melted_sd <- Tab_z1_sd %>% melt
colnames(melted_sd) <- c('group','ASV','RA')

######average the ASVs
Tab_sum_sd <- acast(data = melted_sd, formula = ASV~group,
                 fun.aggregate = mean, value.var = 'RA')

######cluster ASVs
mclust_asv_sd <- hclust(d=as.dist(1-cor(Tab_sum_sd %>% t)), method = 'ward.D2')
mclust_samples_sd <- hclust(d=as.dist(1-cor(Tab_sum_sd)), method = 'ward.D2')

######plot cluster
plot(mclust_asv_sd)
plot(mclust_samples_sd)

#####Divide the dendogram in cluster
df_sd <- mclust_asv_sd %>% dendextend::cutree(tree=., k=6) %>%
  data.frame(ASV=names(.), ClustASV=paste0("Cl", .),
             row.names = NULL)

#remove first column
df_sd <- df_sd[,-1]

######order the ASVs
order_sd <- mclust_asv_sd$order %>% mclust_asv_sd$labels[.]

######melt Tab_sum object
melted2_sd <- Tab_sum_sd %>% melt
colnames(melted2_sd) <- c('ASV', 'group', 'RA_zscore')

#####create a column with soil and NumLeaf 
melted2_sd$group <- melted2_sd$group %>% as.character
melted2_sd$Soil <- sapply(strsplit(melted2_sd$group, '_'), `[`, 1)
melted2_sd$NumLeaf <- sapply(strsplit(melted2_sd$group, '_'), `[`, 2)

#####change the zscore values to <0 when the value is negative
melted2_sd$RA_zscore[which(melted2_sd$RA_zscore < 0)] <- 0

######create a new dataset
merged3_sd <- merge(melted2_sd, df_sd, by='ASV')

######select the ASV with zscore higher than 0.5
merged4_sd <- merged3_sd %>%
  subset(RA_zscore > 0.5) %>% droplevels

######merge merged4_sd with Dat_sd$Tax
df_tax_sd <- Dat_rar_sd$Tax

######rename the first column as ASV
colnames(df_tax_sd)[1] <- 'ASV'

######merge the two objects
merged4_sd <- merge(merged4_sd, df_tax_sd, by='ASV')

######taxonomy name
df_tax_genus_sd <- merged4_sd$Taxonomy %>% as.character %>%
  strsplit(split  = "\\;") %>%
  unlist %>%  gsub(pattern = "[a-z]__",replacement = "",) %>%
  gsub(pattern = " ",replacement = "") %>%
  matrix(data = .,ncol = 7,byrow = T) %>% as.data.frame
colnames(df_tax_genus_sd) <-c("Root","Kingdom","Phylum","Class",
                              "Order", "Family", "Genus")
df_tax_genus_sd <- cbind(merged4_sd$Taxonomy,df_tax_genus_sd)
colnames(df_tax_genus_sd)[1] <- "Taxonomy"
df_tax_genus_sd$Root <- df_tax_genus_sd$Root %>% factor
df_tax_genus_sd$Kingdom <- df_tax_genus_sd$Kingdom %>% factor
df_tax_genus_sd$Phylum <- df_tax_genus_sd$Phylum %>% factor
df_tax_genus_sd$Class <- df_tax_genus_sd$Class %>% factor
df_tax_genus_sd$Order <- df_tax_genus_sd$Order %>% factor
df_tax_genus_sd$Family <- df_tax_genus_sd$Family %>% factor
df_tax_genus_sd$Genus <- df_tax_genus_sd$Genus %>% factor

#######combine the datasets
Res_final_sd <- cbind(merged4_sd, df_tax_genus_sd)

######################select only 50FC and t##############################
Dat_raw_t <- Dat_raw %>% subset.Dataset((Treatment == '50FC') &
                                          (NumLeaf != 'Soil') &
                                          (NumLeaf != 'Root') &
                                          (NumLeaf != 'L6') &
                                          (Soil == 'Tarrafal'), 
                                        drop = T, clean = T)

####Define the RA
tab=Dat_raw_t$Tab
#######Calculate the relative abundance
###understand sweep function
#sweep(x, MARGIN, STATS, FUN="-", check.margin=T, ...)
# x is the data
#Whether you operate by row (2) or column(1) is defined by MARGIN
# STATS refers to the summary statistics which you wish to 'sweep out'
# FUN is the function used to carry out the sweep, "-" is the default
tab_ra=sweep(tab,2,colSums(tab),"/")
##sweep is diving each row value by the column sums.In another word,
#for each row evaluated STATS. This value represents the 
#relative abundance
##Relative abundance is the percent composition of an ASV in 
#the samples relative to the total number of ASV in the samples
tab_c=tab_ra>0.001

######Melt the RA table
melted_prev <- tab_c %>% melt
colnames(melted_prev) <- c("ASV", "DADA2_Header","Present")

######remove all the present = False
melted_prev_t <- melted_prev %>%
  subset(Present == TRUE) %>% droplevels

melted_prev_t <- merge(melted_prev_t,Dat_raw_t$Map, by = "DADA2_Header")

####Splits the data into subsets and computes statistics for each, and returns
##Count the number of Present ASVs per soil and per leaves
df_sum_t <- aggregate(Present ~ ASV+Soil+NumLeaf,
                    melted_prev_t,sum)

Mapa_t <- Dat_raw_t$Map
Mapa_t$Count <- 1

###Calculate the number of samples/compartment
df_tot_t <- aggregate(Count~Soil+NumLeaf,Mapa_t,sum)

df_sum_t <- merge(df_sum_t,df_tot_t, by = c("Soil","NumLeaf"))

df_sum_t$Prevalence <- df_sum_t$Present/df_sum_t$Count

######create a column group
df_sum_t$group <- paste(df_sum_t$Soil, df_sum_t$NumLeaf, df_sum_t$ASV, sep='_')

######select all the ASVs with a prevalence > 0.3
df_sum_03_t <- df_sum_t %>%
  subset(Prevalence > 0.3) %>% droplevels

#Rarefied version
Dat_rar <- dat$Rarefied

######select only 50FC
Dat_rar_t <- Dat_rar %>% subset.Dataset((Treatment == '50FC') &
                                          (NumLeaf != 'Soil') &
                                          (NumLeaf != 'Root') &
                                          (NumLeaf != 'L6') &
                                          (Soil == 'Tarrafal'), 
                                        drop = T, clean = T)

#####Data exploration
Tab_t <- Dat_rar_t$Tab %>% t

######select the prevalence ASV
usable_samples_t <- df_sum_03_t$ASV %>% unique %>% as.character #unique ASV 106
Tab_t <- Tab_t[,colnames(Tab_t) %in% usable_samples_t]

#####Transform the data in z_score
#A z-score is a measure of how many standard deviations 
#(how far my observation is from the mean) below or above
#the populatio mean a raw score is.
Tab_z_t <- Tab_t %>% scale

#####Starts to build the PCA
colMeans(Tab_z_t) %>% sort %>% plot

#####create the Map datasets
Map_t <- Dat_rar_t$Map

######create a group column with DADA2_Header and Leaves
Map_t$group <- paste(Map_t$Soil, Map_t$NumLeaf, sep='_')

######select the column to add the group column with the Tab
group_t <- Map_t[,c(7,10)]

######merge Tab with group
Tab1_t <- data.frame(Tab_z_t)

#######bind the Tab1 with group
Tab_z1_t <- cbind(group_t, Tab1_t)

######change the rowname for Tab_z1
rownames(Tab_z1_t) <- NULL

######remove the first 2 columns from Tab_z1
Tab_z1_t <- Tab_z1_t[,-c(1)]

######melt the dataset. The RA is a zscore
melted_t <- Tab_z1_t %>% melt
colnames(melted_t) <- c('group','ASV','RA')

######average the ASVs
Tab_sum_t <- acast(data = melted_t, formula = ASV~group,
                 fun.aggregate = mean, value.var = 'RA')

######cluster ASVs
mclust_asv_t <- hclust(d=as.dist(1-cor(Tab_sum_t %>% t)), method = 'ward.D2')
mclust_samples_t <- hclust(d=as.dist(1-cor(Tab_sum_t)), method = 'ward.D2')

######plot cluster
plot(mclust_asv_t)
plot(mclust_samples_t)

#####Divide the dendogram in cluster
df_t <- mclust_asv_t %>% dendextend::cutree(tree=., k=6) %>%
  data.frame(ASV=names(.), ClustASV=paste0("Cl", .),
             row.names = NULL)

#remove first column
df_t <- df_t[,-1]

######order the ASVs
order_t <- mclust_asv_t$order %>% mclust_asv_t$labels[.]

######melt Tab_sum object
melted2_t <- Tab_sum_t %>% melt
colnames(melted2_t) <- c('ASV', 'group', 'RA_zscore')

#####create a column with soil and NumLeaf 
melted2_t$group <- melted2_t$group %>% as.character
melted2_t$Soil <- sapply(strsplit(melted2_t$group, '_'), `[`, 1)
melted2_t$NumLeaf <- sapply(strsplit(melted2_t$group, '_'), `[`, 2)

#####change the zscore values to <0 when the value is negative
melted2_t$RA_zscore[which(melted2_t$RA_zscore < 0)] <- 0

######create a new dataset
merged3_t <- merge(melted2_t, df_t, by='ASV')

######select the ASV with zscore higher than 0.5
merged4_t <- merged3_t %>%
  subset(RA_zscore > 0.5) %>% droplevels

######merge merged4_t with Dat_t$Tax
df_tax_t <- Dat_rar_t$Tax

######rename the first column as ASV
colnames(df_tax_t)[1] <- 'ASV'

######merge the two objects
merged4_t <- merge(merged4_t, df_tax_t, by='ASV')

######taxonomy name
df_tax_genus_t <- merged4_t$Taxonomy %>% as.character %>%
  strsplit(split  = "\\;") %>%
  unlist %>%  gsub(pattern = "[a-z]__",replacement = "",) %>%
  gsub(pattern = " ",replacement = "") %>%
  matrix(data = .,ncol = 7,byrow = T) %>% as.data.frame
colnames(df_tax_genus_t) <-c("Root","Kingdom","Phylum","Class",
                             "Order", "Family", "Genus")
df_tax_genus_t <- cbind(merged4_t$Taxonomy,df_tax_genus_t)
colnames(df_tax_genus_t)[1] <- "Taxonomy"
df_tax_genus_t$Root <- df_tax_genus_t$Root %>% factor
df_tax_genus_t$Kingdom <- df_tax_genus_t$Kingdom %>% factor
df_tax_genus_t$Phylum <- df_tax_genus_t$Phylum %>% factor
df_tax_genus_t$Class <- df_tax_genus_t$Class %>% factor
df_tax_genus_t$Order <- df_tax_genus_t$Order %>% factor
df_tax_genus_t$Family <- df_tax_genus_t$Family %>% factor
df_tax_genus_t$Genus <- df_tax_genus_t$Genus %>% factor

#######combine the datasets
Res_final_t <- cbind(merged4_t, df_tax_genus_t)

######drop the zeros
merged4_sb <- Res_final_sb %>%
  subset(RA_zscore >0) %>% droplevels

merged4_sd <- Res_final_sd %>%
  subset(RA_zscore >0) %>% droplevels

merged4_t <- Res_final_t %>%
  subset(RA_zscore >0) %>% droplevels

######select each leaf
df_sb_l1 <- merged4_sb %>%
  subset(NumLeaf == 'L1') %>% droplevels
df_sb_l2 <- merged4_sb %>%
  subset(NumLeaf == 'L2') %>% droplevels
df_sb_l3 <- merged4_sb %>%
  subset(NumLeaf == 'L3') %>% droplevels
df_sb_l4 <- merged4_sb %>%
  subset(NumLeaf == 'L4') %>% droplevels
df_sb_l5 <- merged4_sb %>%
  subset(NumLeaf == 'L5') %>% droplevels

df_sd_l1 <- merged4_sd %>%
  subset(NumLeaf == 'L1') %>% droplevels
df_sd_l2 <- merged4_sd %>%
  subset(NumLeaf == 'L2') %>% droplevels
df_sd_l3 <- merged4_sd %>%
  subset(NumLeaf == 'L3') %>% droplevels
df_sd_l4 <- merged4_sd %>%
  subset(NumLeaf == 'L4') %>% droplevels
df_sd_l5 <- merged4_sd %>%
  subset(NumLeaf == 'L5') %>% droplevels

df_t_l1 <- merged4_t %>%
  subset(NumLeaf == 'L1') %>% droplevels
df_t_l2 <- merged4_t %>%
  subset(NumLeaf == 'L2') %>% droplevels
df_t_l3 <- merged4_t %>%
  subset(NumLeaf == 'L3') %>% droplevels
df_t_l4 <- merged4_t %>%
  subset(NumLeaf == 'L4') %>% droplevels
df_t_l5 <- merged4_t %>%
  subset(NumLeaf == 'L5') %>% droplevels

######select the id per soil
ids_sb_l1_05 <- df_sb_l1$Genus %>% unique
ids_sd_l1_05 <- df_sd_l1$Genus %>% unique
ids_t_l1_05 <- df_t_l1$Genus %>% unique

######intersect the ids objects
intersect_sb_sd_05 <- intersect(ids_sb_l1_05, ids_sd_l1_05) 
intersect_sb_t_05 <- intersect(ids_sb_l1_05, ids_t_l1_05)
intersect_sd_t_05 <- intersect(ids_sd_l1_05, ids_t_l1_05)

######calculate proportion
prop_sb_sd <- ((length(intersect_sb_sd_05))/(length(ids_sb_l1_05)+length(ids_sd_l1_05)))*100
prop_sb_t <- ((length(intersect_sb_t_05))/(length(ids_sb_l1_05)+length(ids_t_l1_05)))*100
prop_sd_t <- ((length(intersect_sd_t_05))/(length(ids_sd_l1_05)+length(ids_t_l1_05)))*100

######create a dataframe for the result
overlap_l1_05 <- data.frame(numGenus = c(length(intersect_sb_sd_05), 
                                         length(intersect_sb_t_05),
                                         length(intersect_sd_t_05)),
                            prop_genus = c(prop_sb_sd,prop_sb_t,prop_sd_t),
                            contrast = c('Sutton Bonington vs São Domingos', 
                                         'Sutton Bonington vs Tarrafal', 
                                         'São Domingos vs Tarrafal'))

######add the NumLeaf column
overlap_l1_05$NumLeaf <- 'L1'

######the proportion is significant in L1
prop_l1 <- prop.test(x = prop_sb_sd, 
                     n = prop_sb_t,
                     alternative = "greater",
                     correct = TRUE)

######add the pvalue
prop_l1 <- prop_l1$p.value %>% data.frame
colnames(prop_l1)[1] <- 'pvalue'
prop_l1$contrast <- 'Sutton Bonington vs São Domingos'
prop_l1$NumLeaf <- 'L1'

######leaf 2
ids_sb_l2_05 <- df_sb_l2$Genus %>% unique
ids_sd_l2_05 <- df_sd_l2$Genus %>% unique
ids_t_l2_05 <- df_t_l2$Genus %>% unique

######intersect the ids objects
intersect_sb_sd_l2 <- intersect(ids_sb_l2_05, ids_sd_l2_05) 
intersect_sb_t_l2 <- intersect(ids_sb_l2_05, ids_t_l2_05)
intersect_sd_t_l2 <- intersect(ids_sd_l2_05, ids_t_l2_05)

######calculate proportion
prop_sb_sd_l2 <- ((length(intersect_sb_sd_l2))/(length(ids_sb_l2_05)+length(ids_sd_l2_05)))*100
prop_sb_t_l2 <- ((length(intersect_sb_t_l2))/(length(ids_sb_l2_05)+length(ids_t_l2_05)))*100
prop_sd_t_l2 <- ((length(intersect_sd_t_l2))/(length(ids_sd_l2_05)+length(ids_t_l2_05)))*100

######create a dataframe for the result
overlap_l2_05 <- data.frame(numGenus = c(length(intersect_sb_sd_l2), 
                                         length(intersect_sb_t_l2),
                                         length(intersect_sd_t_l2)),
                            prop_genus = c(prop_sb_sd_l2,prop_sb_t_l2,prop_sd_t_l2),
                            contrast = c('Sutton Bonington vs São Domingos', 
                                         'Sutton Bonington vs Tarrafal', 
                                         'São Domingos vs Tarrafal'))

######add the NumLeaf column
overlap_l2_05$NumLeaf <- 'L2'

######the proportion is significant in L2
prop_l2 <- prop.test(x = c(length(intersect_sb_sd_l2),
                           length(intersect_sb_t_l2)),
                     n = c(length(ids_sb_l2_05)+length(ids_sd_l2_05),
                           length(ids_sb_l2_05)+length(ids_t_l2_05)),
                           alternative = "greater", correct = FALSE)

######add the pvalue
prop_l2 <- prop_l2$p.value %>% data.frame
colnames(prop_l2)[1] <- 'pvalue'
prop_l2$contrast <- 'Sutton Bonington vs São Domingos'
prop_l2$NumLeaf <- 'L2'

######leaf 3
ids_sb_l3_05 <- df_sb_l3$Genus %>% unique
ids_sd_l3_05 <- df_sd_l3$Genus %>% unique
ids_t_l3_05 <- df_t_l3$Genus %>% unique

######intersect the ids objects
intersect_sb_sd_l3 <- intersect(ids_sb_l3_05, ids_sd_l3_05) 
intersect_sb_t_l3 <- intersect(ids_sb_l3_05, ids_t_l3_05)
intersect_sd_t_l3 <- intersect(ids_sd_l3_05, ids_t_l3_05)

######calculate proportion
prop_sb_sd_l3 <- ((length(intersect_sb_sd_l3))/(length(ids_sb_l3_05)+length(ids_sd_l3_05)))*100
prop_sb_t_l3 <- ((length(intersect_sb_t_l3))/(length(ids_sb_l3_05)+length(ids_t_l3_05)))*100
prop_sd_t_l3 <- ((length(intersect_sd_t_l3))/(length(ids_sd_l3_05)+length(ids_t_l3_05)))*100

######create a dataframe for the result
overlap_l3_05 <- data.frame(numGenus = c(length(intersect_sb_sd_l3), 
                                         length(intersect_sb_t_l3),
                                         length(intersect_sd_t_l3)),
                            prop_genus = c(prop_sb_sd_l3,prop_sb_t_l3,prop_sd_t_l3),
                            contrast = c('Sutton Bonington vs São Domingos', 
                                         'Sutton Bonington vs Tarrafal', 
                                         'São Domingos vs Tarrafal'))
library(stats)
fcr1 <- phyper(length(intersect_sb_sd_l3), 23,20,23,
               lower.tail = FALSE)

######the proportion is significant in L1
prop_l3 <- prop.test(x = c(length(intersect_sb_sd_l3),
                           length(intersect_sb_t_l3)),
                     n = c(length(ids_sb_l3_05)+length(ids_sd_l3_05),
                           length(ids_sb_l3_05)+length(ids_t_l3_05)),
                     alternative = "greater")

######add the pvalue
prop_l3 <- prop_l3$p.value %>% data.frame
colnames(prop_l3)[1] <- 'pvalue'
prop_l3$contrast <- 'Sutton Bonington vs São Domingos'
prop_l3$NumLeaf <- 'L3'

######add the NumLeaf column
overlap_l3_05$NumLeaf <- 'L3'

######leaf 4
ids_sb_l4_05 <- df_sb_l4$Genus %>% unique
ids_sd_l4_05 <- df_sd_l4$Genus %>% unique
ids_t_l4_05 <- df_t_l4$Genus %>% unique

######intersect the ids objects
intersect_sb_sd_l4 <- intersect(ids_sb_l4_05, ids_sd_l4_05) 
intersect_sb_t_l4 <- intersect(ids_sb_l4_05, ids_t_l4_05)
intersect_sd_t_l4 <- intersect(ids_sd_l4_05, ids_t_l4_05)

######calculate proportion
prop_sb_sd_l4 <- ((length(intersect_sb_sd_l4))/(length(ids_sb_l4_05)+length(ids_sd_l4_05)))*100
prop_sb_t_l4 <- ((length(intersect_sb_t_l4))/(length(ids_sb_l4_05)+length(ids_t_l4_05)))*100
prop_sd_t_l4 <- ((length(intersect_sd_t_l4))/(length(ids_sd_l4_05)+length(ids_t_l4_05)))*100

######create a dataframe for the result
overlap_l4_05 <- data.frame(numGenus = c(length(intersect_sb_sd_l4), 
                                         length(intersect_sb_t_l4),
                                         length(intersect_sd_t_l4)),
                            prop_genus = c(prop_sb_sd_l4,prop_sb_t_l4,prop_sd_t_l4),
                            contrast = c('Sutton Bonington vs São Domingos', 
                                         'Sutton Bonington vs Tarrafal', 
                                         'São Domingos vs Tarrafal'))

######the proportion is significant in L1
prop_l4 <- prop.test(x = c(length(intersect_sb_sd_l4),
                           length(intersect_sb_t_l4)),
                     n = c(length(ids_sb_l4_05)+length(ids_sd_l4_05),
                           length(ids_sb_l4_05)+length(ids_t_l4_05)),
                     alternative = "greater", conf.level = 0.90)

######add the pvalue
prop_l4 <- prop_l4$p.value %>% data.frame
colnames(prop_l4)[1] <- 'pvalue'
prop_l4$contrast <- 'Sutton Bonington vs São Domingos'
prop_l4$NumLeaf <- 'L4'

######add the NumLeaf column
overlap_l4_05$NumLeaf <- 'L4'

######leaf 5
ids_sb_l5_05 <- df_sb_l5$Genus %>% unique
ids_sd_l5_05 <- df_sd_l5$Genus %>% unique
ids_t_l5_05 <- df_t_l5$Genus %>% unique

######intersect the ids objects
intersect_sb_sd_l5 <- intersect(ids_sb_l5_05, ids_sd_l5_05) 
intersect_sb_t_l5 <- intersect(ids_sb_l5_05, ids_t_l5_05)
intersect_sd_t_l5 <- intersect(ids_sd_l5_05, ids_t_l5_05)

######calculate proportion
prop_sb_sd <- ((length(intersect_sb_sd_l5))/(length(ids_sb_l5_05)+length(ids_sd_l5_05)))*100
prop_sb_t <- ((length(intersect_sb_t_l5))/(length(ids_sb_l5_05)+length(ids_t_l5_05)))*100
prop_sd_t <- ((length(intersect_sd_t_l5))/(length(ids_sd_l5_05)+length(ids_t_l5_05)))*100

######create a dataframe for the result
overlap_l5_05 <- data.frame(numGenus = c('104', '106', '124'),
                            prop_genus = c('26.06','25.67','28.18'),
                            contrast = c('Sutton Bonington vs São Domingos', 
                                         'Sutton Bonington vs Tarrafal', 
                                         'São Domingos vs Tarrafal'))
######the proportion is significant in L1
prop_l5 <- prop.test(x = c(length(intersect_sb_sd_l5),
                           length(intersect_sb_t_l5)),
                     n = c(length(ids_sb_l5_05)+length(ids_sd_l5_05),
                           length(ids_sb_l5_05)+length(ids_t_l5_05)),
                     alternative = "greater")

######add the pvalue
prop_l5 <- prop_l5$p.value %>% data.frame
colnames(prop_l5)[1] <- 'pvalue'
prop_l5$contrast <- 'Sutton Bonington vs São Domingos'
prop_l5$NumLeaf <- 'L5'

######add the NumLeaf column
overlap_l5_05$NumLeaf <- 'L5'

######bind the columns
overlap_final_05 <- rbind(overlap_l1_05, overlap_l2_05, overlap_l3_05, 
                          overlap_l4_05, overlap_l5_05)

######bind the rows
prop <- rbind(prop_l1, prop_l2, prop_l3, prop_l4, prop_l5)

######as numeric
overlap_final_05$prop_genus <- overlap_final_05$prop_genus %>% as.numeric
overlap_final_05$numGenus <- overlap_final_05$numGenus %>% as.numeric

######add the significant column
pthres <- 0.05
overlap_final_05$significance <- rep("NoSignificant",nrow(overlap_final_05))
overlap_final_05$significance[which(overlap_final_05$NumLeaf == 'L4' &
                                      overlap_final_05$contrast == 'Sutton Bonington vs São Domingos')] <- "Significant"
overlap_final_05$significance <- overlap_final_05$significance %>% factor

######reorder the contrast 
overlap_final_05$contrast <- overlap_final_05$contrast %>%
  factor(levels = c('Sutton Bonington vs São Domingos',
                    'Sutton Bonington vs Tarrafal',
                    'São Domingos vs Tarrafal'))

######plot the results
p_overlap_05 <- ggplot(data = overlap_final_05, aes(x=NumLeaf, y=prop_genus,
                                                    fill = significance))+
  geom_bar(stat = 'identity')+
  facet_grid(.~contrast)+
  xlab(NULL)+
  ylab('Proportion of shared Genus (%)')+
  ylim(0,30)+
  ggtitle('RA_zscore > 0.5')+
  scale_fill_manual(labels = c('No significant','Significance in respect to Sutton Bonington vs Tarrafal proportion (p-value < 0.05)'),
                    values=c('#999999','#ef8a62'))+
  clean +
  theme(#legend.position = c(0, 73),
    legend.key.height= unit(0.3, 'cm'),
    legend.key.width = unit(0.6, 'cm'),
    legend.box.just = "right",
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    legend.position = 'bottom',
    strip.background = element_blank(),
    strip.text.x = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=20),
    plot.title = element_text(hjust = 0.5, size=20))

p_overlap_05

######save the plot as pdf
oh.save.pdf(p = p_overlap_05,outname = "figS2F-proportion_test_shared_genus_05.pdf",
            outdir = "../figures/",width = 18,height = 10)

