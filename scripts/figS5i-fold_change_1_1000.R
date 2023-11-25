######Load packages
library(ohchibi)
library(multcomp)
library(emmeans)

######Set seed
set.seed(130816)

######Set work directory
setwd("C:/Users/Castrillo Lab/OneDrive - The University of Nottingham/Publication 2019/Microbiome/Mutants_screening/scripts")
source("0-Clean_up_plots.R")
######create color for the dataset
paleta_syncom <- c('#80cdc1','#018571')
names(paleta_syncom) <- c('NB','+SynCom')

######Load dataset
df_1000 <- read.csv('../rawdata/arabidopsis_mutants_1_1000_final_full_nb_fold_change_dunnet.csv',
                   sep=',')

########rename the tair
df_1000$tair[which((df_1000$sample_id == '14'))] <- 'AT5G65710_#1'
df_1000$tair[which((df_1000$sample_id == '35'))] <- 'AT5G65710_#2'

######select the tair 
tair_id <- df_1000$tair %>% unique %>% as.character

######create an empty dataset
Res_em <- NULL

######calculate the mean for each tair
for (i in tair_id) {
  df_col <- df_1000 %>% subset(tair == i) %>% droplevels
  m1_col <- lm(formula = Area_mm_2 ~ treatment,data = df_col)
  m1_em_col <- emmeans(m1_col,pairwise ~ treatment,ref = 1,adjust = "none")
  df_em_col <- m1_em_col$emmeans %>% as.data.frame
  df_em_col$tair  <- i
  Res_em <- rbind(Res_em, df_em_col)
}

######group tair and treatment
Res_em$group <- paste(Res_em$tair, Res_em$treatment, sep='_')

#######select each tair
df_col_nb <- Res_em %>%
  subset((tair == 'Col-0') & (treatment == 'NB')) %>% droplevels
df_col_sync <- Res_em %>%
  subset((tair == 'Col-0') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_col_nb)[2] <- 'emmeans_nb_col'
colnames(df_col_sync)[2] <- 'emmeans_syn_col'

######combine the two datasets by column
df_col <- cbind(df_col_nb, df_col_sync)

######calculate the foldchange
df_col$fold <- df_col$emmeans_syn_col/df_col$emmeans_nb_col

######select the columns
df_col_fold <- df_col[,c(7,17)]

###########################################select each tair
df_AT5G19200_nb <- Res_em %>%
  subset((tair == 'AT5G19200') & (treatment == 'NB')) %>% droplevels
df_AT5G19200_sync <- Res_em %>%
  subset((tair == 'AT5G19200') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT5G19200_nb)[2] <- 'emmeans_nb_AT5G19200'
colnames(df_AT5G19200_sync)[2] <- 'emmeans_syn_AT5G19200'

######combine the two datasets by column
df_AT5G19200 <- cbind(df_AT5G19200_nb, df_AT5G19200_sync)

######calculate the foldchange
df_AT5G19200$fold <- df_AT5G19200$emmeans_syn_AT5G19200/df_AT5G19200$emmeans_nb_AT5G19200

######select the columns
df_AT5G19200_fold <- df_AT5G19200[,c(7,17)]

###########################################select each tair
df_AT4G20110_nb <- Res_em %>%
  subset((tair == 'AT4G20110') & (treatment == 'NB')) %>% droplevels
df_AT4G20110_sync <- Res_em %>%
  subset((tair == 'AT4G20110') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT4G20110_nb)[2] <- 'emmeans_nb_AT4G20110'
colnames(df_AT4G20110_sync)[2] <- 'emmeans_syn_AT4G20110'

######combine the two datasets by column
df_AT4G20110 <- cbind(df_AT4G20110_nb, df_AT4G20110_sync)

######calculate the foldchange
df_AT4G20110$fold <- df_AT4G20110$emmeans_syn_AT4G20110/df_AT4G20110$emmeans_nb_AT4G20110

######select the columns
df_AT4G20110_fold <- df_AT4G20110[,c(7,17)]

##########################
df_AT4G27970_nb <- Res_em %>%
  subset((tair == 'AT4G27970') & (treatment == 'NB')) %>% droplevels
df_AT4G27970_sync <- Res_em %>%
  subset((tair == 'AT4G27970') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT4G27970_nb)[2] <- 'emmeans_nb_AT4G27970'
colnames(df_AT4G27970_sync)[2] <- 'emmeans_syn_AT4G27970'

######combine the two datasets by column
df_AT4G27970 <- cbind(df_AT4G27970_nb, df_AT4G27970_sync)

######calculate the foldchange
df_AT4G27970$fold <- df_AT4G27970$emmeans_syn_AT4G27970/df_AT4G27970$emmeans_nb_AT4G27970

######select the columns
df_AT4G27970_fold <- df_AT4G27970[,c(7,17)]

#############################AT4G18030
df_AT4G18030_nb <- Res_em %>%
  subset((tair == 'AT4G18030') & (treatment == 'NB')) %>% droplevels
df_AT4G18030_sync <- Res_em %>%
  subset((tair == 'AT4G18030') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT4G18030_nb)[2] <- 'emmeans_nb_AT4G18030'
colnames(df_AT4G18030_sync)[2] <- 'emmeans_syn_AT4G18030'

######combine the two datasets by column
df_AT4G18030 <- cbind(df_AT4G18030_nb, df_AT4G18030_sync)

######calculate the foldchange
df_AT4G18030$fold <- df_AT4G18030$emmeans_syn_AT4G18030/df_AT4G18030$emmeans_nb_AT4G18030

######select the columns
df_AT4G18030_fold <- df_AT4G18030[,c(7,17)]

#############################AT5G65710
df_AT5G65710_nb <- Res_em %>%
  subset((tair == 'AT5G65710') & (treatment == 'NB')) %>% droplevels
df_AT5G65710_sync <- Res_em %>%
  subset((tair == 'AT5G65710') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT5G65710_nb)[2] <- 'emmeans_nb_AT5G65710'
colnames(df_AT5G65710_sync)[2] <- 'emmeans_syn_AT5G65710'

######combine the two datasets by column
df_AT5G65710 <- cbind(df_AT5G65710_nb, df_AT5G65710_sync)

######calculate the foldchange
df_AT5G65710$fold <- df_AT5G65710$emmeans_syn_AT5G65710/df_AT5G65710$emmeans_nb_AT5G65710

######select the columns
df_AT5G65710_fold <- df_AT5G65710[,c(7,17)]

#############################AT3G53840
df_AT3G53840_nb <- Res_em %>%
  subset((tair == 'AT3G53840') & (treatment == 'NB')) %>% droplevels
df_AT3G53840_sync <- Res_em %>%
  subset((tair == 'AT3G53840') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT3G53840_nb)[2] <- 'emmeans_nb_AT3G53840'
colnames(df_AT3G53840_sync)[2] <- 'emmeans_syn_AT3G53840'

######combine the two datasets by column
df_AT3G53840 <- cbind(df_AT3G53840_nb, df_AT3G53840_sync)

######calculate the foldchange
df_AT3G53840$fold <- df_AT3G53840$emmeans_syn_AT3G53840/df_AT3G53840$emmeans_nb_AT3G53840

######select the columns
df_AT3G53840_fold <- df_AT3G53840[,c(7,17)]

#############################AT2G45830
df_AT2G45830_nb <- Res_em %>%
  subset((tair == 'AT2G45830') & (treatment == 'NB')) %>% droplevels
df_AT2G45830_sync <- Res_em %>%
  subset((tair == 'AT2G45830') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT2G45830_nb)[2] <- 'emmeans_nb_AT2G45830'
colnames(df_AT2G45830_sync)[2] <- 'emmeans_syn_AT2G45830'

######combine the two datasets by column
df_AT2G45830 <- cbind(df_AT2G45830_nb, df_AT2G45830_sync)

######calculate the foldchange
df_AT2G45830$fold <- df_AT2G45830$emmeans_syn_AT2G45830/df_AT2G45830$emmeans_nb_AT2G45830

######select the columns
df_AT2G45830_fold <- df_AT2G45830[,c(7,17)]

#############################AT3G09030
df_AT3G09030_nb <- Res_em %>%
  subset((tair == 'AT3G09030') & (treatment == 'NB')) %>% droplevels
df_AT3G09030_sync <- Res_em %>%
  subset((tair == 'AT3G09030') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT3G09030_nb)[2] <- 'emmeans_nb_AT3G09030'
colnames(df_AT3G09030_sync)[2] <- 'emmeans_syn_AT3G09030'

######combine the two datasets by column
df_AT3G09030 <- cbind(df_AT3G09030_nb, df_AT3G09030_sync)

######calculate the foldchange
df_AT3G09030$fold <- df_AT3G09030$emmeans_syn_AT3G09030/df_AT3G09030$emmeans_nb_AT3G09030

######select the columns
df_AT3G09030_fold <- df_AT3G09030[,c(7,17)]

#############################AT3G48980
df_AT3G48980_nb <- Res_em %>%
  subset((tair == 'AT3G48980') & (treatment == 'NB')) %>% droplevels
df_AT3G48980_sync <- Res_em %>%
  subset((tair == 'AT3G48980') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT3G48980_nb)[2] <- 'emmeans_nb_AT3G48980'
colnames(df_AT3G48980_sync)[2] <- 'emmeans_syn_AT3G48980'

######combine the two datasets by column
df_AT3G48980 <- cbind(df_AT3G48980_nb, df_AT3G48980_sync)

######calculate the foldchange
df_AT3G48980$fold <- df_AT3G48980$emmeans_syn_AT3G48980/df_AT3G48980$emmeans_nb_AT3G48980

######select the columns
df_AT3G48980_fold <- df_AT3G48980[,c(7,17)]

#############################AT2G45910
df_AT2G45910_nb <- Res_em %>%
  subset((tair == 'AT2G45910') & (treatment == 'NB')) %>% droplevels
df_AT2G45910_sync <- Res_em %>%
  subset((tair == 'AT2G45910') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT2G45910_nb)[2] <- 'emmeans_nb_AT2G45910'
colnames(df_AT2G45910_sync)[2] <- 'emmeans_syn_AT2G45910'

######combine the two datasets by column
df_AT2G45910 <- cbind(df_AT2G45910_nb, df_AT2G45910_sync)

######calculate the foldchange
df_AT2G45910$fold <- df_AT2G45910$emmeans_syn_AT2G45910/df_AT2G45910$emmeans_nb_AT2G45910

######select the columns
df_AT2G45910_fold <- df_AT2G45910[,c(7,17)]

#############################AT1G48880
df_AT1G48880_nb <- Res_em %>%
  subset((tair == 'AT1G48880') & (treatment == 'NB')) %>% droplevels
df_AT1G48880_sync <- Res_em %>%
  subset((tair == 'AT1G48880') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT1G48880_nb)[2] <- 'emmeans_nb_AT1G48880'
colnames(df_AT1G48880_sync)[2] <- 'emmeans_syn_AT1G48880'

######combine the two datasets by column
df_AT1G48880 <- cbind(df_AT1G48880_nb, df_AT1G48880_sync)

######calculate the foldchange
df_AT1G48880$fold <- df_AT1G48880$emmeans_syn_AT1G48880/df_AT1G48880$emmeans_nb_AT1G48880

######select the columns
df_AT1G48880_fold <- df_AT1G48880[,c(7,17)]

#############################AT3G61270
df_AT1G48880_nb <- Res_em %>%
  subset((tair == 'AT1G48880') & (treatment == 'NB')) %>% droplevels
df_AT1G48880_sync <- Res_em %>%
  subset((tair == 'AT1G48880') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT1G48880_nb)[2] <- 'emmeans_nb_AT1G48880'
colnames(df_AT1G48880_sync)[2] <- 'emmeans_syn_AT1G48880'

######combine the two datasets by column
df_AT1G48880 <- cbind(df_AT1G48880_nb, df_AT1G48880_sync)

######calculate the foldchange
df_AT1G48880$fold <- df_AT1G48880$emmeans_syn_AT1G48880/df_AT1G48880$emmeans_nb_AT1G48880

######select the columns
df_AT1G48880_fold <- df_AT1G48880[,c(7,17)]

#############################AT3G61270
df_AT3G61270_nb <- Res_em %>%
  subset((tair == 'AT3G61270') & (treatment == 'NB')) %>% droplevels
df_AT3G61270_sync <- Res_em %>%
  subset((tair == 'AT3G61270') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT3G61270_nb)[2] <- 'emmeans_nb_AT3G61270'
colnames(df_AT3G61270_sync)[2] <- 'emmeans_syn_AT3G61270'

######combine the two datasets by column
df_AT3G61270 <- cbind(df_AT3G61270_nb, df_AT3G61270_sync)

######calculate the foldchange
df_AT3G61270$fold <- df_AT3G61270$emmeans_syn_AT3G61270/df_AT3G61270$emmeans_nb_AT3G61270

######select the columns
df_AT3G61270_fold <- df_AT3G61270[,c(7,17)]

#############################AT4G19080
df_AT4G19080_nb <- Res_em %>%
  subset((tair == 'AT4G19080') & (treatment == 'NB')) %>% droplevels
df_AT4G19080_sync <- Res_em %>%
  subset((tair == 'AT4G19080') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT4G19080_nb)[2] <- 'emmeans_nb_AT4G19080'
colnames(df_AT4G19080_sync)[2] <- 'emmeans_syn_AT4G19080'

######combine the two datasets by column
df_AT4G19080 <- cbind(df_AT4G19080_nb, df_AT4G19080_sync)

######calculate the foldchange
df_AT4G19080$fold <- df_AT4G19080$emmeans_syn_AT4G19080/df_AT4G19080$emmeans_nb_AT4G19080

######select the columns
df_AT4G19080_fold <- df_AT4G19080[,c(7,17)]

#############################AT3G05050
df_AT3G05050_nb <- Res_em %>%
  subset((tair == 'AT3G05050') & (treatment == 'NB')) %>% droplevels
df_AT3G05050_sync <- Res_em %>%
  subset((tair == 'AT3G05050') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT3G05050_nb)[2] <- 'emmeans_nb_AT3G05050'
colnames(df_AT3G05050_sync)[2] <- 'emmeans_syn_AT3G05050'

######combine the two datasets by column
df_AT3G05050 <- cbind(df_AT3G05050_nb, df_AT3G05050_sync)

######calculate the foldchange
df_AT3G05050$fold <- df_AT3G05050$emmeans_syn_AT3G05050/df_AT3G05050$emmeans_nb_AT3G05050

######select the columns
df_AT3G05050_fold <- df_AT3G05050[,c(7,17)]

#############################AT3G02840
df_AT3G02840_nb <- Res_em %>%
  subset((tair == 'AT3G02840') & (treatment == 'NB')) %>% droplevels
df_AT3G02840_sync <- Res_em %>%
  subset((tair == 'AT3G02840') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT3G02840_nb)[2] <- 'emmeans_nb_AT3G02840'
colnames(df_AT3G02840_sync)[2] <- 'emmeans_syn_AT3G02840'

######combine the two datasets by column
df_AT3G02840 <- cbind(df_AT3G02840_nb, df_AT3G02840_sync)

######calculate the foldchange
df_AT3G02840$fold <- df_AT3G02840$emmeans_syn_AT3G02840/df_AT3G02840$emmeans_nb_AT3G02840

######select the columns
df_AT3G02840_fold <- df_AT3G02840[,c(7,17)]

#############################AT5G07100
df_AT5G07100_nb <- Res_em %>%
  subset((tair == 'AT5G07100') & (treatment == 'NB')) %>% droplevels
df_AT5G07100_sync <- Res_em %>%
  subset((tair == 'AT5G07100') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT5G07100_nb)[2] <- 'emmeans_nb_AT5G07100'
colnames(df_AT5G07100_sync)[2] <- 'emmeans_syn_AT5G07100'

######combine the two datasets by column
df_AT5G07100 <- cbind(df_AT5G07100_nb, df_AT5G07100_sync)

######calculate the foldchange
df_AT5G07100$fold <- df_AT5G07100$emmeans_syn_AT5G07100/df_AT5G07100$emmeans_nb_AT5G07100

######select the columns
df_AT5G07100_fold <- df_AT5G07100[,c(7,17)]

#############################AT5G10770
df_AT5G10770_nb <- Res_em %>%
  subset((tair == 'AT5G10770') & (treatment == 'NB')) %>% droplevels
df_AT5G10770_sync <- Res_em %>%
  subset((tair == 'AT5G10770') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT5G10770_nb)[2] <- 'emmeans_nb_AT5G10770'
colnames(df_AT5G10770_sync)[2] <- 'emmeans_syn_AT5G10770'

######combine the two datasets by column
df_AT5G10770 <- cbind(df_AT5G10770_nb, df_AT5G10770_sync)

######calculate the foldchange
df_AT5G10770$fold <- df_AT5G10770$emmeans_syn_AT5G10770/df_AT5G10770$emmeans_nb_AT5G10770

######select the columns
df_AT5G10770_fold <- df_AT5G10770[,c(7,17)]

#############################AT5G54430
df_AT5G54430_nb <- Res_em %>%
  subset((tair == 'AT5G54430') & (treatment == 'NB')) %>% droplevels
df_AT5G54430_sync <- Res_em %>%
  subset((tair == 'AT5G54430') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT5G54430_nb)[2] <- 'emmeans_nb_AT5G54430'
colnames(df_AT5G54430_sync)[2] <- 'emmeans_syn_AT5G54430'

######combine the two datasets by column
df_AT5G54430 <- cbind(df_AT5G54430_nb, df_AT5G54430_sync)

######calculate the foldchange
df_AT5G54430$fold <- df_AT5G54430$emmeans_syn_AT5G54430/df_AT5G54430$emmeans_nb_AT5G54430

######select the columns
df_AT5G54430_fold <- df_AT5G54430[,c(7,17)]

#############################AT2G02070
df_AT2G02070_nb <- Res_em %>%
  subset((tair == 'AT2G02070') & (treatment == 'NB')) %>% droplevels
df_AT2G02070_sync <- Res_em %>%
  subset((tair == 'AT2G02070') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT2G02070_nb)[2] <- 'emmeans_nb_AT2G02070'
colnames(df_AT2G02070_sync)[2] <- 'emmeans_syn_AT2G02070'

######combine the two datasets by column
df_AT2G02070 <- cbind(df_AT2G02070_nb, df_AT2G02070_sync)

######calculate the foldchange
df_AT2G02070$fold <- df_AT2G02070$emmeans_syn_AT2G02070/df_AT2G02070$emmeans_nb_AT2G02070

######select the columns
df_AT2G02070_fold <- df_AT2G02070[,c(7,17)]

#############################AT5G22380
df_AT5G22380_nb <- Res_em %>%
  subset((tair == 'AT5G22380') & (treatment == 'NB')) %>% droplevels
df_AT5G22380_sync <- Res_em %>%
  subset((tair == 'AT5G22380') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT5G22380_nb)[2] <- 'emmeans_nb_AT5G22380'
colnames(df_AT5G22380_sync)[2] <- 'emmeans_syn_AT5G22380'

######combine the two datasets by column
df_AT5G22380 <- cbind(df_AT5G22380_nb, df_AT5G22380_sync)

######calculate the foldchange
df_AT5G22380$fold <- df_AT5G22380$emmeans_syn_AT5G22380/df_AT5G22380$emmeans_nb_AT5G22380

######select the columns
df_AT5G22380_fold <- df_AT5G22380[,c(7,17)]

#############################AT5G25400
df_AT5G25400_nb <- Res_em %>%
  subset((tair == 'AT5G25400') & (treatment == 'NB')) %>% droplevels
df_AT5G25400_sync <- Res_em %>%
  subset((tair == 'AT5G25400') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT5G25400_nb)[2] <- 'emmeans_nb_AT5G25400'
colnames(df_AT5G25400_sync)[2] <- 'emmeans_syn_AT5G25400'

######combine the two datasets by column
df_AT5G25400 <- cbind(df_AT5G25400_nb, df_AT5G25400_sync)

######calculate the foldchange
df_AT5G25400$fold <- df_AT5G25400$emmeans_syn_AT5G25400/df_AT5G25400$emmeans_nb_AT5G25400

######select the columns
df_AT5G25400_fold <- df_AT5G25400[,c(7,17)]

#############################AT5G07120
df_AT5G07120_nb <- Res_em %>%
  subset((tair == 'AT5G07120') & (treatment == 'NB')) %>% droplevels
df_AT5G07120_sync <- Res_em %>%
  subset((tair == 'AT5G07120') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT5G07120_nb)[2] <- 'emmeans_nb_AT5G07120'
colnames(df_AT5G07120_sync)[2] <- 'emmeans_syn_AT5G07120'

######combine the two datasets by column
df_AT5G07120 <- cbind(df_AT5G07120_nb, df_AT5G07120_sync)

######calculate the foldchange
df_AT5G07120$fold <- df_AT5G07120$emmeans_syn_AT5G07120/df_AT5G07120$emmeans_nb_AT5G07120

######select the columns
df_AT5G07120_fold <- df_AT5G07120[,c(7,17)]

#############################AT5G13130
df_AT5G13130_nb <- Res_em %>%
  subset((tair == 'AT5G13130') & (treatment == 'NB')) %>% droplevels
df_AT5G13130_sync <- Res_em %>%
  subset((tair == 'AT5G13130') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT5G13130_nb)[2] <- 'emmeans_nb_AT5G13130'
colnames(df_AT5G13130_sync)[2] <- 'emmeans_syn_AT5G13130'

######combine the two datasets by column
df_AT5G13130 <- cbind(df_AT5G13130_nb, df_AT5G13130_sync)

######calculate the foldchange
df_AT5G13130$fold <- df_AT5G13130$emmeans_syn_AT5G13130/df_AT5G13130$emmeans_nb_AT5G13130

######select the columns
df_AT5G13130_fold <- df_AT5G13130[,c(7,17)]

#############################AT1G62020
df_AT1G62020_nb <- Res_em %>%
  subset((tair == 'AT1G62020') & (treatment == 'NB')) %>% droplevels
df_AT1G62020_sync <- Res_em %>%
  subset((tair == 'AT1G62020') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT1G62020_nb)[2] <- 'emmeans_nb_AT1G62020'
colnames(df_AT1G62020_sync)[2] <- 'emmeans_syn_AT1G62020'

######combine the two datasets by column
df_AT1G62020 <- cbind(df_AT1G62020_nb, df_AT1G62020_sync)

######calculate the foldchange
df_AT1G62020$fold <- df_AT1G62020$emmeans_syn_AT1G62020/df_AT1G62020$emmeans_nb_AT1G62020

######select the columns
df_AT1G62020_fold <- df_AT1G62020[,c(7,17)]

#############################AT4G24970
df_AT4G24970_nb <- Res_em %>%
  subset((tair == 'AT4G24970') & (treatment == 'NB')) %>% droplevels
df_AT4G24970_sync <- Res_em %>%
  subset((tair == 'AT4G24970') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT4G24970_nb)[2] <- 'emmeans_nb_AT4G24970'
colnames(df_AT4G24970_sync)[2] <- 'emmeans_syn_AT4G24970'

######combine the two datasets by column
df_AT4G24970 <- cbind(df_AT4G24970_nb, df_AT4G24970_sync)

######calculate the foldchange
df_AT4G24970$fold <- df_AT4G24970$emmeans_syn_AT4G24970/df_AT4G24970$emmeans_nb_AT4G24970

######select the columns
df_AT4G24970_fold <- df_AT4G24970[,c(7,17)]

#############################AT5G23850
df_AT5G23850_nb <- Res_em %>%
  subset((tair == 'AT5G23850') & (treatment == 'NB')) %>% droplevels
df_AT5G23850_sync <- Res_em %>%
  subset((tair == 'AT5G23850') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT5G23850_nb)[2] <- 'emmeans_nb_AT5G23850'
colnames(df_AT5G23850_sync)[2] <- 'emmeans_syn_AT5G23850'

######combine the two datasets by column
df_AT5G23850 <- cbind(df_AT5G23850_nb, df_AT5G23850_sync)

######calculate the foldchange
df_AT5G23850$fold <- df_AT5G23850$emmeans_syn_AT5G23850/df_AT5G23850$emmeans_nb_AT5G23850

######select the columns
df_AT5G23850_fold <- df_AT5G23850[,c(7,17)]

#############################AT1G55910
df_AT1G55910_nb <- Res_em %>%
  subset((tair == 'AT1G55910') & (treatment == 'NB')) %>% droplevels
df_AT1G55910_sync <- Res_em %>%
  subset((tair == 'AT1G55910') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT1G55910_nb)[2] <- 'emmeans_nb_AT1G55910'
colnames(df_AT1G55910_sync)[2] <- 'emmeans_syn_AT1G55910'

######combine the two datasets by column
df_AT1G55910 <- cbind(df_AT1G55910_nb, df_AT1G55910_sync)

######calculate the foldchange
df_AT1G55910$fold <- df_AT1G55910$emmeans_syn_AT1G55910/df_AT1G55910$emmeans_nb_AT1G55910

######select the columns
df_AT1G55910_fold <- df_AT1G55910[,c(7,17)]

#############################AT5G58440
df_AT5G58440_nb <- Res_em %>%
  subset((tair == 'AT5G58440') & (treatment == 'NB')) %>% droplevels
df_AT5G58440_sync <- Res_em %>%
  subset((tair == 'AT5G58440') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT5G58440_nb)[2] <- 'emmeans_nb_AT5G58440'
colnames(df_AT5G58440_sync)[2] <- 'emmeans_syn_AT5G58440'

######combine the two datasets by column
df_AT5G58440 <- cbind(df_AT5G58440_nb, df_AT5G58440_sync)

######calculate the foldchange
df_AT5G58440$fold <- df_AT5G58440$emmeans_syn_AT5G58440/df_AT5G58440$emmeans_nb_AT5G58440

######select the columns
df_AT5G58440_fold <- df_AT5G58440[,c(7,17)]

#############################AT2G34940
df_AT2G34940_nb <- Res_em %>%
  subset((tair == 'AT2G34940') & (treatment == 'NB')) %>% droplevels
df_AT2G34940_sync <- Res_em %>%
  subset((tair == 'AT2G34940') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT2G34940_nb)[2] <- 'emmeans_nb_AT2G34940'
colnames(df_AT2G34940_sync)[2] <- 'emmeans_syn_AT2G34940'

######combine the two datasets by column
df_AT2G34940 <- cbind(df_AT2G34940_nb, df_AT2G34940_sync)

######calculate the foldchange
df_AT2G34940$fold <- df_AT2G34940$emmeans_syn_AT2G34940/df_AT2G34940$emmeans_nb_AT2G34940

######select the columns
df_AT2G34940_fold <- df_AT2G34940[,c(7,17)]

#############################AT2G30250
df_AT2G30250_nb <- Res_em %>%
  subset((tair == 'AT2G30250') & (treatment == 'NB')) %>% droplevels
df_AT2G30250_sync <- Res_em %>%
  subset((tair == 'AT2G30250') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT2G30250_nb)[2] <- 'emmeans_nb_AT2G30250'
colnames(df_AT2G30250_sync)[2] <- 'emmeans_syn_AT2G30250'

######combine the two datasets by column
df_AT2G30250 <- cbind(df_AT2G30250_nb, df_AT2G30250_sync)

######calculate the foldchange
df_AT2G30250$fold <- df_AT2G30250$emmeans_syn_AT2G30250/df_AT2G30250$emmeans_nb_AT2G30250

######select the columns
df_AT2G30250_fold <- df_AT2G30250[,c(7,17)]

#############################AT1G30900
df_AT1G30900_nb <- Res_em %>%
  subset((tair == 'AT1G30900') & (treatment == 'NB')) %>% droplevels
df_AT1G30900_sync <- Res_em %>%
  subset((tair == 'AT1G30900') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT1G30900_nb)[2] <- 'emmeans_nb_AT1G30900'
colnames(df_AT1G30900_sync)[2] <- 'emmeans_syn_AT1G30900'

######combine the two datasets by column
df_AT1G30900 <- cbind(df_AT1G30900_nb, df_AT1G30900_sync)

######calculate the foldchange
df_AT1G30900$fold <- df_AT1G30900$emmeans_syn_AT1G30900/df_AT1G30900$emmeans_nb_AT1G30900

######select the columns
df_AT1G30900_fold <- df_AT1G30900[,c(7,17)]

#############################AT3G07480
df_AT3G07480_nb <- Res_em %>%
  subset((tair == 'AT3G07480') & (treatment == 'NB')) %>% droplevels
df_AT3G07480_sync <- Res_em %>%
  subset((tair == 'AT3G07480') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT3G07480_nb)[2] <- 'emmeans_nb_AT3G07480'
colnames(df_AT3G07480_sync)[2] <- 'emmeans_syn_AT3G07480'

######combine the two datasets by column
df_AT3G07480 <- cbind(df_AT3G07480_nb, df_AT3G07480_sync)

######calculate the foldchange
df_AT3G07480$fold <- df_AT3G07480$emmeans_syn_AT3G07480/df_AT3G07480$emmeans_nb_AT3G07480

######select the columns
df_AT3G07480_fold <- df_AT3G07480[,c(7,17)]

#############################AT5G38560
df_AT5G38560_nb <- Res_em %>%
  subset((tair == 'AT5G38560') & (treatment == 'NB')) %>% droplevels
df_AT5G38560_sync <- Res_em %>%
  subset((tair == 'AT5G38560') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT5G38560_nb)[2] <- 'emmeans_nb_AT5G38560'
colnames(df_AT5G38560_sync)[2] <- 'emmeans_syn_AT5G38560'

######combine the two datasets by column
df_AT5G38560 <- cbind(df_AT5G38560_nb, df_AT5G38560_sync)

######calculate the foldchange
df_AT5G38560$fold <- df_AT5G38560$emmeans_syn_AT5G38560/df_AT5G38560$emmeans_nb_AT5G38560

######select the columns
df_AT5G38560_fold <- df_AT5G38560[,c(7,17)]

#############################AT1G66160
df_AT1G66160_nb <- Res_em %>%
  subset((tair == 'AT1G66160') & (treatment == 'NB')) %>% droplevels
df_AT1G66160_sync <- Res_em %>%
  subset((tair == 'AT1G66160') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT1G66160_nb)[2] <- 'emmeans_nb_AT1G66160'
colnames(df_AT1G66160_sync)[2] <- 'emmeans_syn_AT1G66160'

######combine the two datasets by column
df_AT1G66160 <- cbind(df_AT1G66160_nb, df_AT1G66160_sync)

######calculate the foldchange
df_AT1G66160$fold <- df_AT1G66160$emmeans_syn_AT1G66160/df_AT1G66160$emmeans_nb_AT1G66160

######select the columns
df_AT1G66160_fold <- df_AT1G66160[,c(7,17)]

#############################AT5G24030
df_AT5G24030_nb <- Res_em %>%
  subset((tair == 'AT5G24030') & (treatment == 'NB')) %>% droplevels
df_AT5G24030_sync <- Res_em %>%
  subset((tair == 'AT5G24030') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT5G24030_nb)[2] <- 'emmeans_nb_AT5G24030'
colnames(df_AT5G24030_sync)[2] <- 'emmeans_syn_AT5G24030'

######combine the two datasets by column
df_AT5G24030 <- cbind(df_AT5G24030_nb, df_AT5G24030_sync)

######calculate the foldchange
df_AT5G24030$fold <- df_AT5G24030$emmeans_syn_AT5G24030/df_AT5G24030$emmeans_nb_AT5G24030

######select the columns
df_AT5G24030_fold <- df_AT5G24030[,c(7,17)]

#############################AT4G11650
df_AT4G11650_nb <- Res_em %>%
  subset((tair == 'AT4G11650') & (treatment == 'NB')) %>% droplevels
df_AT4G11650_sync <- Res_em %>%
  subset((tair == 'AT4G11650') & (treatment == '+SynCom')) %>% droplevels
######rename the column
colnames(df_AT4G11650_nb)[2] <- 'emmeans_nb_AT4G11650'
colnames(df_AT4G11650_sync)[2] <- 'emmeans_syn_AT4G11650'

######combine the two datasets by column
df_AT4G11650 <- cbind(df_AT4G11650_nb, df_AT4G11650_sync)

######calculate the foldchange
df_AT4G11650$fold <- df_AT4G11650$emmeans_syn_AT4G11650/df_AT4G11650$emmeans_nb_AT4G11650

######select the columns
df_AT4G11650_fold <- df_AT4G11650[,c(7,17)]

######combine the dataset
df_final <- rbind(df_col_fold, df_AT5G19200_fold,
                  df_AT4G20110_fold, df_AT4G27970_fold,
                  df_AT4G18030_fold, df_AT5G65710_fold,
                  df_AT3G53840_fold, df_AT2G45830_fold,
                  df_AT3G09030_fold, df_AT3G48980_fold,
                  df_AT2G45910_fold, df_AT1G48880_fold,
                  df_AT3G61270_fold, df_AT4G19080_fold,
                  df_AT3G05050_fold, df_AT3G02840_fold,
                  df_AT5G07100_fold, df_AT5G10770_fold,
                  df_AT5G54430_fold, df_AT2G02070_fold,
                  df_AT5G22380_fold, df_AT5G25400_fold,
                  df_AT5G07120_fold, df_AT5G13130_fold,
                  df_AT1G62020_fold, df_AT4G24970_fold,
                  df_AT5G23850_fold, df_AT1G55910_fold,
                  df_AT5G58440_fold, df_AT2G34940_fold,
                  df_AT2G30250_fold, df_AT1G30900_fold,
                  df_AT3G07480_fold, df_AT5G38560_fold,
                  df_AT1G66160_fold, df_AT5G24030_fold,
                  df_AT4G11650_fold)

######reorder the tair column
df_final$tair <- df_final$tair %>%
  factor(levels=c("Col-0",
                  "AT5G19200", "AT2G45910","AT5G10770","AT3G58490",
                  "AT2G30250","AT5G38560","AT1G66160","AT5G24030",
                  "AT4G11650","AT4G20110", "AT4G27970", "AT4G17950",
                  "AT4G18030", "AT5G65710_#1", "AT3G53840", "AT2G45830", 
                  "AT3G09030",
                  "AT3G48980","AT5G11230", "AT1G48880", "AT3G61270",
                  "AT4G19080", "AT3G05050", "AT2G17610", "AT3G02840",
                  "AT5G07100","AT5G54430","AT2G02070","AT5G22380",
                  "AT5G25400","AT5G65710_#2","AT5G07120","AT5G13130",   
                  "AT1G62020","AT4G24970","AT5G23850","AT1G55910",
                  "AT5G58440","AT2G34940","AT5G18520",
                  "AT1G30900","AT3G07480"))

######create the plot
ggplot(df_final, aes(x=tair, y=fold))+
  geom_col()+
  ylab('Fold change
      with respect to No bacteria')+
  ylim(0.0, 1.5)+
  clean+
  theme(axis.text.x = element_text(size=12,vjust=1, angle=45, hjust=1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.position = 'none')

########significant fold change
######Load dataset
df_fold <- read.csv('../rawdata/arabidopsis_mutants_1_1000_final_full_fold_change_final.csv',
                    sep=',')

######create the Col-0_+SynCom and Col-0_NB
df_fold$tair[which((df_fold$sample_id == '14'))] <- 'AT5G65710_#1'
df_fold$tair[which((df_fold$sample_id == '35'))] <- 'AT5G65710_#2'

######Reorder the tair
df_fold$tair <- df_fold$tair %>%
  factor(levels=c('Col-0',"AT5G19200", "AT2G45910","AT5G10770","AT3G58490",
                  "AT2G30250","AT5G38560","AT1G66160","AT5G24030",
                  "AT4G11650","AT4G20110", "AT4G27970",
                  "AT4G18030", "AT5G65710_#1", "AT3G53840", "AT2G45830", 
                  "AT3G09030",
                  "AT3G48980","AT5G11230", "AT1G48880", "AT3G61270",
                  "AT4G19080", "AT3G05050", "AT2G17610", "AT3G02840",
                  "AT5G07100","AT5G54430","AT2G02070","AT5G22380",
                  "AT5G25400","AT5G65710_#2","AT5G07120","AT5G13130",   
                  "AT1G62020","AT4G24970","AT5G23850","AT1G55910",
                  "AT5G58440","AT2G34940","AT5G18520",
                  "AT1G30900","AT3G07480"))
######change the NB and SynCom in numeric
df_fold$NB <- df_fold$NB %>% as.numeric
df_fold$SynCom <- df_fold$SynCom %>% as.numeric

######calculate the foldchange
df_fold$foldchange <- df_fold$SynCom/df_fold$NB

#######################significance analysis
df_fold_AT5G19200 <- df_fold %>%
  subset(tair == 'AT5G19200') %>% droplevels
t.test(df_fold_AT5G19200$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT2G45910 <- df_fold %>%
  subset(tair == 'AT2G45910') %>% droplevels
t.test(df_fold_AT2G45910$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT5G10770 <- df_fold %>%
  subset(tair == 'AT5G10770') %>% droplevels
t.test(df_fold_AT5G10770$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT2G30250 <- df_fold %>%
  subset(tair == 'AT2G30250') %>% droplevels
t.test(df_fold_AT2G30250$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT5G38560 <- df_fold %>%
  subset(tair == 'AT5G38560') %>% droplevels
t.test(df_fold_AT5G38560$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT1G66160 <- df_fold %>%
  subset(tair == 'AT1G66160') %>% droplevels
t.test(df_fold_AT1G66160$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT5G24030 <- df_fold %>%
  subset(tair == 'AT5G24030') %>% droplevels
t.test(df_fold_AT5G24030$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT4G11650 <- df_fold %>%
  subset(tair == 'AT4G11650') %>% droplevels
t.test(df_fold_AT4G11650$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT4G20110 <- df_fold %>%
  subset(tair == 'AT4G20110') %>% droplevels
t.test(df_fold_AT4G20110$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT4G27970 <- df_fold %>%
  subset(tair == 'AT4G27970') %>% droplevels
t.test(df_fold_AT4G27970$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT4G18030 <- df_fold %>%
  subset(tair == 'AT4G18030') %>% droplevels
t.test(df_fold_AT4G18030$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT5G65710_1 <- df_fold %>%
  subset(tair == 'AT5G65710_#1') %>% droplevels
t.test(df_fold_AT5G65710_1$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT3G53840 <- df_fold %>%
  subset(tair == 'AT3G53840') %>% droplevels
t.test(df_fold_AT3G53840$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT2G45830 <- df_fold %>%
  subset(tair == 'AT2G45830') %>% droplevels
t.test(df_fold_AT2G45830$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT3G09030 <- df_fold %>%
  subset(tair == 'AT3G09030') %>% droplevels
t.test(df_fold_AT3G09030$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT3G48980 <- df_fold %>%
  subset(tair == 'AT3G48980') %>% droplevels
t.test(df_fold_AT3G48980$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT1G48880 <- df_fold %>%
  subset(tair == 'AT1G48880') %>% droplevels
t.test(df_fold_AT1G48880$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT3G61270 <- df_fold %>%
  subset(tair == 'AT3G61270') %>% droplevels
t.test(df_fold_AT3G61270$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT4G19080 <- df_fold %>%
  subset(tair == 'AT4G19080') %>% droplevels
t.test(df_fold_AT4G19080$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT3G05050 <- df_fold %>%
  subset(tair == 'AT3G05050') %>% droplevels
t.test(df_fold_AT3G05050$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT3G02840 <- df_fold %>%
  subset(tair == 'AT3G02840') %>% droplevels
t.test(df_fold_AT3G02840$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT5G07100 <- df_fold %>%
  subset(tair == 'AT5G07100') %>% droplevels
t.test(df_fold_AT5G07100$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT5G54430 <- df_fold %>%
  subset(tair == 'AT5G54430') %>% droplevels
t.test(df_fold_AT5G54430$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT2G02070 <- df_fold %>%
  subset(tair == 'AT2G02070') %>% droplevels
t.test(df_fold_AT2G02070$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT5G22380 <- df_fold %>%
  subset(tair == 'AT5G22380') %>% droplevels
t.test(df_fold_AT5G22380$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT5G25400 <- df_fold %>%
  subset(tair == 'AT5G25400') %>% droplevels
t.test(df_fold_AT5G25400$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT5G65710_2 <- df_fold %>%
  subset(tair == 'AT5G65710_#2') %>% droplevels
t.test(df_fold_AT5G65710_2$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT5G07120 <- df_fold %>%
  subset(tair == 'AT5G07120') %>% droplevels
t.test(df_fold_AT5G07120$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT5G13130 <- df_fold %>%
  subset(tair == 'AT5G13130') %>% droplevels
t.test(df_fold_AT5G13130$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT1G62020 <- df_fold %>%
  subset(tair == 'AT1G62020') %>% droplevels
t.test(df_fold_AT1G62020$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT4G24970 <- df_fold %>%
  subset(tair == 'AT4G24970') %>% droplevels
t.test(df_fold_AT4G24970$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT5G23850 <- df_fold %>%
  subset(tair == 'AT5G23850') %>% droplevels
t.test(df_fold_AT5G23850$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT1G55910 <- df_fold %>%
  subset(tair == 'AT1G55910') %>% droplevels
t.test(df_fold_AT1G55910$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT5G58440 <- df_fold %>%
  subset(tair == 'AT5G58440') %>% droplevels
t.test(df_fold_AT5G58440$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT2G34940 <- df_fold %>%
  subset(tair == 'AT2G34940') %>% droplevels
t.test(df_fold_AT2G34940$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT1G30900 <- df_fold %>%
  subset(tair == 'AT1G30900') %>% droplevels
t.test(df_fold_AT1G30900$foldchange, mu = 1.008092, alternative = "greater")

df_fold_AT3G07480 <- df_fold %>%
  subset(tair == 'AT3G07480') %>% droplevels
t.test(df_fold_AT3G07480$foldchange, mu = 1.008092, alternative = "greater")

tair <- c("Col-0","AT5G19200", "AT2G45910","AT5G10770",
          "AT2G30250","AT5G38560","AT1G66160",
          "AT5G24030","AT4G11650","AT4G20110", 
          "AT4G27970","AT4G18030", "AT5G65710_#1",
          "AT3G53840", "AT2G45830", "AT3G09030",
          "AT3G48980","AT1G48880","AT3G61270",
          "AT4G19080", "AT3G05050", "AT3G02840",
          "AT5G07100","AT5G54430","AT2G02070",
          "AT5G22380","AT5G25400","AT5G65710_#2",
          "AT5G07120","AT5G13130","AT1G62020",
          "AT4G24970","AT5G23850","AT1G55910",
          "AT5G58440","AT2G34940","AT1G30900",
          "AT3G07480")
pvalue <- c(1,0.204, 2.431e-05, 0.1234, 
            0.00418, 0.8568, 0.1139,
            0.9904, 0.006039, 0.5671, 
            0.994, 0.9791, 0.8896,
            0.5722, 0.9408, 0.7585,
            0.9952,0.9851,0.2156,
            0.1185, 0.9011, 0.9997,
            0.6387,0.9757, 0.9872,
            1, 0.9544,0.2964,
            0.002202,0.1645,0.9602,
            1,0.3991,0.2056,
            0.4533,0.9077, 0.583,
            0.1102)
df_pvalue <- data.frame(tair, pvalue)

#Read significant results
pthres <- 0.05
df_pvalue$significance <- rep("NoSignificant",nrow(df_pvalue))
df_pvalue$significance[which(df_pvalue$pvalue < pthres)] <- "Significant"

######merge the two datasets
df <- merge(df_final, df_pvalue, by='tair')

#######create the final plot
p <- ggplot(df, aes(x=tair, y=fold, fill = significance))+
  geom_col()+
  ylab('Fold change with respect to No Bacteria (NB)')+
  ylim(0.0, 1.5)+
  scale_fill_manual(name = 'Significance',labels = c('No significant',
                               'Significant in respect to Col-0 (p < 0.05)'),
                    values=c('#999999', '#d8daeb'))+
  clean+
  theme(axis.text.x = element_text(size=12,vjust=1, angle=45, hjust=1,
                                   color = c('black', 'black', 
                                             'red', 'black',
                                             'red', 'black',
                                             'black', 'black',
                                             'black', 'black',
                                             'black', 'black',
                                             'black', 'black',
                                             'red', 'black',
                                             'black', 'black',
                                             'black', 'black',
                                             'black', 'black',
                                             'black', 'black',
                                             'black', 'black',
                                             'red', 'black',
                                             'black', 'black',
                                             'black', 'black',
                                             'black', 'black',
                                             'black', 'black')),
        legend.key.height= unit(0.3, 'cm'),
        legend.key.width = unit(0.6, 'cm'),
        axis.title.x = element_blank(),
        legend.position = 'right')

p

######save figures
oh.save.pdf(p = p,outname = "fig5_foldchange_arabidopsis_1_1000_final.pdf",
            outdir = "../figures/",width = 20,height = 10)
