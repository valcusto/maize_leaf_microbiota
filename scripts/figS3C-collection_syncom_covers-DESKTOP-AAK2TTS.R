######Load packages
library(ohchibi)
library(seqinr)
library(dplyr)

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
Dat_rab_sb <- Dat_rab %>% subset.Dataset(Treatment == '50FC', 
                                         drop = T, clean = T)

######remove the soil samples
Dat_rab_sb <- Dat_rab_sb %>% subset.Dataset(NumLeaf != 'Soil', 
                                         drop = T, clean = T)

#Collapse the ASV data at order levels
Dat_order <- Dat_rab_sb %>% collapse_by_taxonomy.Dataset(Dat = ., 
                                                      level = 5)

###Step5:Average the relative abundance
df_freq <- Dat_order$Tab %>% rowMeans %>% data.frame(Id = names(.),RA = .,
                                                     row.names = NULL)

###Step6:Select the 100 most abundant ASV
#Strategy:
#1.Order the ASVs for the highest to the lowest
df_freq <- with(df_freq,order(-RA)) %>% df_freq[.,]
#2.Choose the most abundant ASVs
mchosen <- df_freq$Id[1:12] %>% as.character
Tab_chosen <- which(rownames(Dat_order$Tab) %in% mchosen) %>%
  Dat_order$Tab[.,] %>% as.matrix
#3.Define the low abundant ASV as others
Tab_other <- which(!(rownames(Dat_order$Tab) %in% mchosen)) %>%
  Dat_order$Tab[.,] %>% as.matrix
Other <- colSums(Tab_other)
#4.Merge the chosen and the others object
Tab_chosen <- rbind(Tab_chosen,Other)

melted <- Tab_chosen %>% melt
colnames(melted) <- c("Id","DADA2_Header","RA")

#####Prepare the data frame to merge
Dat_order$Map$DADA2_Header <- rownames(Dat_order$Map)

melted <- merge(melted,Dat_order$Map, by = "DADA2_Header")
melted$Approach <- "Natural\ncommunities"

##remove the Root and Bacteria character on the Id
melted$Id <- melted$Id %>% gsub(pattern = "Root; k__Bacteria;",
                                replacement = "")
melted$Id <- melted$Id %>% gsub(pattern = " p__Actinobacteria; c__Actinobacteria; o__",
                                replacement = "")
melted$Id <- melted$Id %>% gsub(pattern = " p__Bacteroidetes; c__Bacteroidia; o__",
                                replacement = "")
melted$Id <- melted$Id %>% gsub(pattern = " p__Firmicutes; c__Bacilli; o__",
                                replacement = "")
melted$Id <- melted$Id %>% gsub(pattern = " p__Proteobacteria; c__Alphaproteobacteria; o__",
                                replacement = "")
melted$Id <- melted$Id %>% gsub(pattern = " p__Proteobacteria; c__Gammaproteobacteria; o__",
                                replacement = "")

##Levels the ID
melted$Id <- melted$Id %>% factor(levels = c("Bacillales",
                                                     "Betaproteobacteriales",
                                                     "Chitinophagales", 
                                                     "Enterobacteriales",
                                                     "Micrococcales",
                                                     "Propionibacteriales",
                                                     "Pseudomonadales",
                                                     "Pseudonocardiales",
                                                     "Rhizobiales",
                                                     "Sphingomonadales",
                                                     "Streptomycetales",
                                                     "Xanthomonadales",
                                                     "Other"))
####Define the colors for each phyla or class
mOrder <- c("#1F77B4","#AEC7E8","#FF7F0E","#FFBB78","#2CA02C","#98DF8A",
            "#D62728","#FF9D9A","#9467BD","#C5B0D5","#8C564B", "#C49C94",
            "#F7B6D2","#C7C7C7", "#7F7F7F", "#E377C2")
names(mOrder) <- c("Bacillales","Betaproteobacteriales","Burkholderiales",
                   "Chitinophagales","Enterobacterales","Micrococcales",
                   "Propionibacteriales","Pseudomonadales","Pseudonocardiales",
                   "Rhizobiales","Sphingobacteriales","Sphingomonadales",
                   "Streptomycetales","Xanthomonadales", "Other", "Staphylococcales")

###Step7:Plot for the culture independent
p_natural <- ggplot(data = melted,aes(Approach,RA)) +
  geom_bar(stat = "identity",aes(fill = Id),position = "fill") +
  scale_fill_manual(values = mOrder,na.value = "#D9D9D9") +
  theme_ohchibi() +
  scale_y_continuous(breaks = seq(0,1,by = 0.1),expand = c(0,0), 
                     labels = scales::percent) +
  scale_x_discrete(expand = c(0,0)) +
  theme(legend.position = 'none',
        axis.text.x = element_text(family = "Arial",size = 20),
        axis.text.y = element_text(family = "Arial",size = 12),
        axis.title.y = element_text(family = "Arial",size = 30)) +
  xlab(label = element_blank())+
  ylab(label = "Relative abundance")

p_natural

###Step8:Count the numeber of ASVs and calculate the percentage
df_asv <- Dat_rab_sb$Tax
##Strategy:
#1.Total of ASV
Total <- nrow(df_asv) %>% as.numeric

######separate Taxonomy names in different column
mdf <- Dat_rab_sb$Tax
df_tax <- mdf$Taxonomy %>% as.character %>%
  strsplit(split  = "\\;") %>%
  unlist %>%  gsub(pattern = "[a-z]__",replacement = "",) %>%
  gsub(pattern = " ",replacement = "") %>%
  matrix(data = .,ncol = 7,byrow = T) %>% as.data.frame
colnames(df_tax) <-c("Root","Kingdom","Phylum",
                     "Class","Order","Family","Genus")
df_tax <- cbind(mdf,df_tax)
colnames(df_tax)[1] <- "ASV"
df_tax$ASV <- df_tax$ASV %>% factor
df_tax$Root <- df_tax$Root %>% factor
df_tax$Kingdom <- df_tax$Kingdom %>% factor
df_tax$Phylum <- df_tax$Phylum %>% factor
df_tax$Class <- df_tax$Class %>% factor
df_tax$Order <- df_tax$Order %>% factor
df_tax$Family <- df_tax$Family %>% factor
df_tax$Genus <- df_tax$Genus %>% factor

#4.Phyla_class as factor
df_tax$Order <- df_tax$Order %>% factor

#3.Create a new column with the new order
df_tax$new_order <- "Other"
df_tax$new_order[df_tax$Order == "Bacillales"] <- "Bacillales"
df_tax$new_order[df_tax$Order == "Betaproteobacteriales"] <- "Betaproteobacteriales"
df_tax$new_order[df_tax$Order== "Chitinophagales"] <- "Chitinophagales"
df_tax$new_order[df_tax$Order == "Enterobacteriales"] <- "Enterobacteriales"
df_tax$new_order[df_tax$Order == "Micrococcales"] <- "Micrococcales"
df_tax$new_order[df_tax$Order == "Propionibacteriales"] <- "Propionibacteriales"
df_tax$new_order[df_tax$Order == "Pseudomonadales"] <- "Pseudomonadales"
df_tax$new_order[df_tax$Order == "Pseudonocardiales"] <- "Pseudonocardiales"
df_tax$new_order[df_tax$Order == "Rhizobiales"] <- "Rhizobiales"
df_tax$new_order[df_tax$Order == "Sphingomonadales"] <-"Sphingomonadales"
df_tax$new_order[df_tax$Order == "Streptomycetales"] <-"Streptomycetales"
df_tax$new_order[df_tax$Order == "Xanthomonadales"] <-"Xanthomonadales"

#4.Phyla_class as factor
df_tax$new_order <- df_tax$new_order %>% factor

#5.Assign the number of ASVs in each class
df_bac <- df_tax %>%
  subset(new_order == "Bacillales") %>% droplevels
Total_bac <- nrow(df_bac) %>% as.numeric

df_beta <- df_tax %>%
  subset(new_order == "Betaproteobacteriales") %>% droplevels
Total_beta <- nrow(df_beta) %>% as.numeric

df_chiti <- df_tax %>%
  subset(new_order == "Chitinophagales") %>% droplevels
Total_chiti <- nrow(df_chiti) %>% as.numeric

df_ente <- df_tax %>%
  subset(new_order == "Enterobacteriales") %>% droplevels
Total_ente <- nrow(df_ente) %>% as.numeric

df_micro <- df_tax %>%
  subset(new_order == "Micrococcales") %>% droplevels
Total_micro <- nrow(df_micro) %>% as.numeric

df_other <- df_tax %>%
  subset(new_order == "Other") %>% droplevels
Total_other <- nrow(df_other) %>% as.numeric

df_propi <- df_tax %>%
  subset(new_order == "Propionibacteriales") %>% droplevels
Total_propi <- nrow(df_propi) %>% as.numeric

df_pseudo <- df_tax %>%
  subset(new_order == "Pseudomonadales") %>% droplevels
Total_pseudo <- nrow(df_pseudo) %>% as.numeric

df_pseudon <- df_tax %>%
  subset(new_order == "Pseudonocardiales") %>% droplevels
Total_pseudon <- nrow(df_pseudon) %>% as.numeric

df_rhizo <- df_tax %>%
  subset(new_order == "Rhizobiales") %>% droplevels
Total_rhizo <- nrow(df_rhizo) %>% as.numeric

df_sphingo <- df_tax %>%
  subset(new_order == "Sphingomonadales") %>% droplevels
Total_sphingo <- nrow(df_sphingo) %>% as.numeric

df_strep <- df_tax %>%
  subset(new_order == "Streptomycetales") %>% droplevels
Total_strep <- nrow(df_strep) %>% as.numeric

df_xantho <- df_tax %>%
  subset(new_order == "Xanthomonadales") %>% droplevels
Total_xantho <- nrow(df_xantho) %>% as.numeric

#6.Create the data frame for the counting
Order <- c("Bacillales","Betaproteobacteriales","Chitinophagales",
           "Enterobacterales",
           "Micrococcales","Propionibacteriales","Pseudomonadales",
           "Pseudonocardiales",
           "Rhizobiales","Sphingomonadales",
           "Streptomycetales","Xanthomonadales", "Other")
Percentage <- c(Total_bac/Total, Total_beta/Total,
                Total_chiti/Total, Total_ente/Total,
                Total_micro/Total, Total_propi/Total,
                Total_pseudo/Total, Total_pseudon/Total,
                Total_rhizo/Total, Total_sphingo/Total,
                Total_strep/Total, Total_xantho/Total, Total_other/Total)

df <- data.frame(Order, Percentage)

#7. Add the Approach in the data frame
df$Order <- df$Order %>% factor
df$Approach <- "Natural\ncommunities"

#########Prepare the data for culture dependent
collection <- read.csv("../cleandata/10.Zm_syncom.csv", 
                        sep = ",")

######factor the Genus
collection$Order <- collection$Order %>% factor

##Step2: Count the number of isolates
Total_collection <- nrow(collection) %>% as.numeric

#Count isolates
col_bac <- collection %>%
  subset(Order == "Bacillales") %>% droplevels
Tt_bac <- nrow(col_bac) %>% as.numeric

col_burk <- collection %>%
  subset(Order == "Burkholderiales") %>% droplevels
Tt_burk <- nrow(col_burk) %>% as.numeric

col_ente <- collection %>%
  subset(Order == "Enterobacterales") %>% droplevels
Tt_ente <- nrow(col_ente) %>% as.numeric

col_micro <- collection %>%
  subset(Order == "Micrococcales") %>% droplevels
Tt_micro <- nrow(col_micro) %>% as.numeric

col_propi <- collection %>%
  subset(Order == "Propionibacteriales") %>% droplevels
Tt_propi <- nrow(col_propi) %>% as.numeric

col_pseudo <- collection %>%
  subset(Order == "Pseudomonadales") %>% droplevels
Tt_pseudo <- nrow(col_pseudo) %>% as.numeric

col_rhizo <- collection %>%
  subset(Order == "Rhizobiales") %>% droplevels
Tt_rhizo <- nrow(col_rhizo) %>% as.numeric

col_sphingoba <- collection %>%
  subset(Order == "Sphingobacteriales") %>% droplevels
Tt_sphingoba <- nrow(col_sphingoba) %>% as.numeric

col_sphingo <- collection %>%
  subset(Order == "Sphingomonadales") %>% droplevels
Tt_sphingo <- nrow(col_sphingo) %>% as.numeric

col_staphy <- collection %>%
  subset(Order == "Staphylococcales") %>% droplevels
Tt_staphy <- nrow(col_staphy) %>% as.numeric

col_strep <- collection %>%
  subset(Order == "Streptomycetales") %>% droplevels
Tt_strep <- nrow(col_strep) %>% as.numeric

col_xantho <- collection %>%
  subset(Order == "Xanthomonadales") %>% droplevels
Tt_xantho <- nrow(col_xantho) %>% as.numeric

##Step3: Create data frame with all the results
Order <- c("Bacillales","Burkholderiales","Enterobacterales",
           "Micrococcales","Propionibacteriales","Pseudomonadales",
           "Rhizobiales","Sphingobacteriales","Sphingomonadales",
           "Staphylococcales","Streptomycetales","Xanthomonadales")
Percentage <- c(Tt_bac/Total_collection, Tt_burk/Total_collection,
                Tt_ente/Total_collection, Tt_micro/Total_collection,
                Tt_propi/Total_collection, Tt_pseudo/Total_collection,
                Tt_rhizo/Total_collection, Tt_sphingoba/Total_collection,
                Tt_sphingo/Total_collection, Tt_staphy/Total_collection,
                Tt_strep/Total_collection, Tt_xantho/Total_collection)

df_collection <- data.frame(Order, Percentage)

df_collection$Approach <- "Collection"

######Syncom analysis
syncom <- collection %>% subset(Syncom == 'Yes') %>% droplevels

######factor the Genus
syncom$Order <- syncom$Order %>% factor

##Step2: Count the number of isolates
Total_syncom <- nrow(syncom) %>% as.numeric

#Count isolates
syn_bac <- syncom %>%
  subset(Order == "Bacillales") %>% droplevels
T_bac <- nrow(syn_bac) %>% as.numeric

syn_burk <- syncom %>%
  subset(Order == "Burkholderiales") %>% droplevels
T_burk <- nrow(syn_burk) %>% as.numeric

syn_ente <- syncom %>%
  subset(Order == "Enterobacterales") %>% droplevels
T_ente <- nrow(syn_ente) %>% as.numeric

syn_micro <- syncom %>%
  subset(Order == "Micrococcales") %>% droplevels
T_micro <- nrow(syn_micro) %>% as.numeric

syn_propi <- syncom %>%
  subset(Order == "Propionibacteriales") %>% droplevels
T_propi <- nrow(syn_propi) %>% as.numeric

syn_pseudo <- syncom %>%
  subset(Order == "Pseudomonadales") %>% droplevels
T_pseudo <- nrow(syn_pseudo) %>% as.numeric

syn_rhizo <- syncom %>%
  subset(Order == "Rhizobiales") %>% droplevels
T_rhizo <- nrow(syn_rhizo) %>% as.numeric

syn_sphingoba <- syncom %>%
  subset(Order == "Sphingobacteriales") %>% droplevels
T_sphingoba <- nrow(syn_sphingoba) %>% as.numeric

syn_sphingo <- syncom %>%
  subset(Order == "Sphingomonadales") %>% droplevels
T_sphingo <- nrow(syn_sphingo) %>% as.numeric

syn_staphy <- syncom %>%
  subset(Order == "Staphylococcales") %>% droplevels
T_staphy <- nrow(syn_staphy) %>% as.numeric

syn_strep <- syncom %>%
  subset(Order == "Streptomycetales") %>% droplevels
T_strep <- nrow(syn_strep) %>% as.numeric

syn_xantho <- syncom %>%
  subset(Order == "Xanthomonadales") %>% droplevels
T_xantho <- nrow(syn_xantho) %>% as.numeric

##Step3: Create data frame with all the results
Order <- c("Bacillales","Burkholderiales","Enterobacterales",
           "Micrococcales","Propionibacteriales","Pseudomonadales",
           "Rhizobiales","Sphingobacteriales","Sphingomonadales",
           "Staphylococcales","Streptomycetales","Xanthomonadales")
Percentage <- c(T_bac/Total_syncom, T_burk/Total_syncom,
                T_ente/Total_syncom, T_micro/Total_syncom,
                T_propi/Total_syncom, T_pseudo/Total_syncom,
                T_rhizo/Total_syncom, T_sphingoba/Total_syncom,
                T_sphingo/Total_syncom,T_staphy/Total_syncom,
                T_strep/Total_syncom, T_xantho/Total_syncom)

df_syncom <- data.frame(Order, Percentage)

df_syncom$Approach <- "SynCom"

######Combine the two datasets
df_legend <- rbind(df,df_collection, df_syncom)
df <- rbind(df_collection, df_syncom)

#####Prepare the file for plot
#Orde the Approach
df$Approach <- df$Approach %>% 
  factor(levels = c("Natural\ncommunities", "Collection", "SynCom"))
df_legend$Approach <- df_legend$Approach %>% 
  factor(levels = c("Natural\ncommunities", "Collection", "SynCom"))

######organize the Order
df$Order <- df$Order %>% 
  factor(levels = c('Bacillales', 'Betaproteobacteriales', 'Burkholderiales',
                    'Chitinophagales', 'Enterobacterales', 'Micrococcales',
                    'Propionibacteriales','Pseudomonadales', 'Pseudonocardiales',
                    'Rhizobiales', 'Sphingobacteriales','Sphingomonadales',
                    'Staphylococcales','Streptomycetales','Xanthomonadales', 
                    'Other'))

df_legend$Order <- df_legend$Order %>% 
  factor(levels = c('Bacillales', 'Betaproteobacteriales', 'Burkholderiales',
                    'Chitinophagales', 'Enterobacterales', 'Micrococcales',
                    'Propionibacteriales','Pseudomonadales', 'Pseudonocardiales',
                    'Rhizobiales', 'Sphingobacteriales','Sphingomonadales',
                    'Staphylococcales','Streptomycetales','Xanthomonadales', 
                    'Other'))

######Compute the plot
p <- ggplot(data = df,aes(Approach,Percentage)) +
  geom_bar(stat = "identity",aes(fill = Order),position = "fill") +
  scale_fill_manual(values = mOrder, drop=FALSE) +
  theme_ohchibi() +
  scale_y_continuous(breaks = seq(0,1,by = 0.1),expand = c(0,0),
                     labels = scales::percent) +
  scale_x_discrete(expand = c(0,0)) +
  theme(
    axis.text.x = element_text(family = "Arial",size = 20),
    axis.text.y = element_text(family = "Arial",size = 12),
    axis.title.y = element_text(family = "Arial",size = 30)) +
  ylab(label = "% of Isolate") +
  xlab(label = element_blank())

p

#####composition
composition <- egg::ggarrange(p_natural, p, nrow = 1, widths = c(1,2))

######save the pdf
oh.save.pdf(p = composition,
            outname = "figS3B-collection_syncom_coverage.pdf",
            outdir = "../figures/",width = 15,height = 25)

# Set the plot margins
par(mar = c(10, 10, 8, 4) + 0.1)

######save as png
png(filename = "../figure/figS3B-collection_syncom_coverage.png", units='cm',
    width = 25, height = 30, res = 300)
composition <- egg::ggarrange(p_natural, p, nrow = 1, widths = c(1,2))
dev.off()
