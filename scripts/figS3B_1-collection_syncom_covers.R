######Load packages
library(ohchibi)
library(seqinr)
library(dplyr)
library(tidyr)

######set seed
set.seed(12356)

######set work directory
setwd("C:/Users/Castrillo Lab/OneDrive - The University of Nottingham/Publication 2019/Paper preparation/manuscripts/scripts")

######helper functions
source('0-Clean_up_plots.R')
source('0_4-subset_Dataset.R')
source('0_8-collapse_by_taxonomy.R')

###Step1:Load the Dat file
Dat <- readRDS('../cleandata/dat_asvs_maize_leaves_rep2.RDS')

###Step2:Select the Relative abundance dataset
Dat_rab <- Dat$RelativeAbundance

######select only the 50% FC
Dat_rab_sub <- Dat_rab %>% subset.Dataset(Treatment == '50FC', 
                                         drop = T, clean = T)

######correct the name for the burkholderia and staphylococcus
mdf <- Dat_rab_sub$Tax
df_tax <- mdf$Taxonomy %>% as.character %>%
  strsplit(split  = "\\;") %>%
  unlist %>%  gsub(pattern = "[a-z]__",replacement = "",) %>%
  gsub(pattern = " ",replacement = "") %>%
  matrix(data = .,ncol = 7,byrow = T) %>% as.data.frame
colnames(df_tax) <-c("Root","Kingdom","Phylum",
                     "Class","Order","Family","Genus")

######add the mdf to the df_tax
df_tax <- cbind(mdf, df_tax)

######select the ASV and the rest
df_tax <- df_tax[,c(3:ncol(df_tax))]

######create column
df_tax$ASV <- row.names(df_tax)

######clean the row names
rownames(df_tax) <- NULL

######check burkholderia and staphylococcus on the dataset
burk <- df_tax %>% subset(Family == 'Burkholderiaceae') %>% droplevels
staphylo <- df_tax %>% subset(Genus == 'Staphylococcus') %>% droplevels

######correct the order for burkholderia and staphylococcus
df_tax$Order[which(df_tax$Family == 'Burkholderiaceae')] <- 'Burkholderiales'
df_tax$Order[which(df_tax$Genus == 'Staphylococcus')] <- 'Staphylococcales'

######merge the columns
df_tax <- df_tax %>%
  unite(tax, Kingdom, Phylum, Class, Order, Family, Genus, 
        sep = ";", remove = F)

######create the quime nomenclature
temp_new<-NULL
for(element in df_tax$tax){
  temp_vec<-unlist(strsplit(x = element,split = ";"))
  qiime_format<-paste("Root; k__",temp_vec[1],"; p__",temp_vec[2],"; c__",temp_vec[3],"; o__",temp_vec[4],"; f__",temp_vec[5],"; g__",temp_vec[6],sep="")
  temp_new<-c(temp_new,qiime_format)
}
df_tax$Taxonomy <- factor(temp_new)

######select the column 9 and 10
df_tax <- df_tax[,c(9,10)]

######select each part of the Dat_rab_sub dataset
Tab <- Dat_rab_sub$Tab %>% t
Map <- Dat_rab_sub$Map

######match the ids
df_tax <- df_tax[match(colnames(Tab),df_tax$ASV),]
row.names(df_tax) <- df_tax$ASV

######create the dataset
Dat <- create_dataset(Tab = Tab %>% t,Map = Map, Tax = df_tax)

######Collapse the ASVs at the order levels
Dat_order <- Dat %>% collapse_by_taxonomy.Dataset(Dat = ., 
                                                         level = 5)

######select the different fractions
Dat_soil <- Dat_order %>% subset.Dataset(NumLeaf == 'Soil', 
                                          drop = T, clean = T)

Dat_root <- Dat_order %>% subset.Dataset(NumLeaf == 'Root', 
                                           drop = T, clean = T)

Dat_leaf <- Dat_order %>% subset.Dataset((NumLeaf != 'Soil') & 
                                             (NumLeaf != 'Root'), 
                                           drop = T, clean = T)

######melt the datasets for each fractions
melted_soil <- Dat_soil$Tab %>% melt
melted_root <- Dat_root$Tab %>% melt
melted_leaf <- Dat_leaf$Tab %>% melt

######name the column in the datasets
colnames(melted_soil) <- c('order', 'DADA2_Header', 'RA')
colnames(melted_root) <- c('order', 'DADA2_Header', 'RA')
colnames(melted_leaf) <- c('order', 'DADA2_Header', 'RA')

######merge the Map
Map_soil <- Dat_soil$Map
Map_root <- Dat_root$Map
Map_leaf <- Dat_leaf$Map

#####merge the objects
merged_soil <- merge(melted_soil, Map_soil, by='DADA2_Header')%>%
  separate(order,sep=';',remove = FALSE, into = c("Root","Kingdom","Phylum",
                                                   "Class","Order")) %>%
  mutate(Kingdom=gsub(pattern = "[a-z]__| ",replacement = "",x=Kingdom)) %>%
  mutate(Phylum=gsub(pattern = "[a-z]__| ",replacement = "",x=Phylum)) %>%
  mutate(Class=gsub(pattern = "[a-z]__| ",replacement = "",x=Class)) %>%
  mutate(Order=gsub(pattern = "[a-z]__| ",replacement = "",x=Order))

merged_root <- merge(melted_root, Map_root, by='DADA2_Header')%>%
  separate(order,sep=';',remove = FALSE, into = c("Root","Kingdom","Phylum",
                                                  "Class","Order")) %>%
  mutate(Kingdom=gsub(pattern = "[a-z]__| ",replacement = "",x=Kingdom)) %>%
  mutate(Phylum=gsub(pattern = "[a-z]__| ",replacement = "",x=Phylum)) %>%
  mutate(Class=gsub(pattern = "[a-z]__| ",replacement = "",x=Class)) %>%
  mutate(Order=gsub(pattern = "[a-z]__| ",replacement = "",x=Order))

merged_leaf <- merge(melted_leaf, Map_leaf, by='DADA2_Header')%>%
  separate(order,sep=';',remove = FALSE, into = c("Root","Kingdom","Phylum",
                                                  "Class","Order")) %>%
  mutate(Kingdom=gsub(pattern = "[a-z]__| ",replacement = "",x=Kingdom)) %>%
  mutate(Phylum=gsub(pattern = "[a-z]__| ",replacement = "",x=Phylum)) %>%
  mutate(Class=gsub(pattern = "[a-z]__| ",replacement = "",x=Class)) %>%
  mutate(Order=gsub(pattern = "[a-z]__| ",replacement = "",x=Order))
  
######select for the SynCom
collection <- read.csv("../cleandata/TableS1.csv", 
                   sep = ",")

######select the order in the different fractions
soil_collection <- intersect(collection$Order, merged_soil$Order)
root_collection <- intersect(collection$Order, merged_root$Order)
leaf_collection <- intersect(collection$Order, merged_leaf$Order)

######select the order in each fractions
overlap_soil_collection <- which(merged_soil$Order %in% soil_collection) %>%
  merged_soil[.,]
overlap_root_collection <- which(merged_root$Order %in% root_collection) %>%
  merged_root[.,]
overlap_leaf_collection <- which(merged_leaf$Order %in% leaf_collection) %>%
  merged_leaf[.,]

######sum the relative abundance
df_soil_collection_sum <- aggregate(RA ~ Order, overlap_soil_collection, mean)
df_root_collection_sum <- aggregate(RA ~ Order, overlap_root_collection, mean)
df_leaf_collection_sum <- aggregate(RA ~ Order, overlap_leaf_collection, mean)

######add the fraction column
df_soil_collection_sum$Fraction <- 'Soil'
df_root_collection_sum$Fraction <- 'Root'
df_leaf_collection_sum$Fraction <- 'Leaf'

######combine all the fractions
df_fraction_collection_sum <- rbind(df_soil_collection_sum, 
                                df_root_collection_sum, df_leaf_collection_sum)

######sum the RA
df_collection_sum <- aggregate(RA~Fraction, df_fraction_collection_sum, sum)

######add the syncom column
df_collection_sum$Collection <- 'Collection'

######transform the fraction in factor
df_collection_sum$Fraction <- df_collection_sum$Fraction %>%
  factor(levels = c('Soil', 'Root', 'Leaf'))

######transform the RA in percentage
df_collection_sum$percentage <- (df_collection_sum$RA)*100

######combine the two datasets
#df <- rbind(df_collection, df_syncom)

######create the plot
p_cumulative_sum <- ggplot(data = df_collection_sum,aes(Fraction,percentage)) +
  geom_bar(stat = "identity") +
  ylim(0,100)+
  theme_ohchibi() +
  theme(
    axis.text.x = element_text(family = "Arial",size = 30),
    axis.text.y = element_text(family = "Arial",size = 30),
    axis.title.y = element_text(family = "Arial",size = 30)) +
  ylab(label = "Total cumulative relative abundance (%)") +
  xlab(label = element_blank())

p_cumulative_sum

######select for the SynCom
syncom <- read.csv("../cleandata/TableS2.csv", 
                   sep = ",")

######select the order in the different fractions
soil_syncom <- intersect(syncom$Order, merged_soil$Order)
root_syncom <- intersect(syncom$Order, merged_root$Order)
leaf_syncom <- intersect(syncom$Order, merged_leaf$Order)

######select the order in each fractions
overlap_soil_syncom <- which(merged_soil$Order %in% soil_syncom) %>%
  merged_soil[.,]
overlap_root_syncom <- which(merged_root$Order %in% root_syncom) %>%
  merged_root[.,]
overlap_leaf_syncom <- which(merged_leaf$Order %in% leaf_syncom) %>%
  merged_leaf[.,]

######average the relative abundance
df_soil_syncom <- aggregate(RA ~ Order, overlap_soil_syncom, mean)
df_root_syncom <- aggregate(RA ~ Order, overlap_root_syncom, mean)
df_leaf_syncom <- aggregate(RA ~ Order, overlap_leaf_syncom, mean)

######add the fraction column
df_soil_syncom$Fraction <- 'Soil'
df_root_syncom$Fraction <- 'Root'
df_leaf_syncom$Fraction <- 'Leaf'

######combine all the fractions
df_fraction_syncom <- rbind(df_soil_syncom, df_root_syncom, df_leaf_syncom)

######sum the RA
df_syncom <- aggregate(RA~Fraction, df_fraction_syncom, sum)

######add the syncom column
df_syncom$Collection <- 'SynCom'

######transform the fraction in factor
df_syncom$Fraction <- df_syncom$Fraction %>%
  factor(levels = c('Soil', 'Root', 'Leaf'))

######transform the RA in percentage
df_syncom$percentage <- (df_syncom$RA)*100

######combine the two datasets
#df <- rbind(df_collection, df_syncom)

######create the plot
p_cumulative <- ggplot(data = df_syncom,aes(Fraction,percentage)) +
  geom_bar(stat = "identity") +
  ylim(0,100)+
  theme_ohchibi() +
  theme(
    axis.text.x = element_text(family = "Arial",size = 30),
    axis.text.y = element_text(family = "Arial",size = 30),
    axis.title.y = element_text(family = "Arial",size = 30)) +
  ylab(label = "Total cumulative relative abundance (%)") +
  xlab(label = element_blank())

p_cumulative

oh.save.pdf(p = p_cumulative,
           outname = "figS3B-SynCom_cumulative_RA.pdf",
           outdir = "../figures/",width = 15,height = 25)