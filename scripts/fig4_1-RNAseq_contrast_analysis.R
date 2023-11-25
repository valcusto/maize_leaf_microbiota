library(ohchibi)
library(dplyr)
library(DESeq2)
library(tibble)

set.seed(130816)


Dats_rnaseq <- readRDS(file = "../cleandata/dat_rnaseq_syncom_do.RDS")


dds <- Dats_rnaseq$dds_R1R2_bpFALSE

### Contrasts ####

list_comparisons <- list(
  
  c("group","L3_L3","L3_FullSynCom"),
  c("group","L3_NB","L3_FullSynCom"),
  c("group","L3_L4_up","L3_FullSynCom"),
  c("group","L3_L5","L3_FullSynCom"),
  
  c("group","L4_L3","L4_FullSynCom"),
  c("group","L4_NB","L4_FullSynCom"),
  c("group","L4_L4_up","L4_FullSynCom"),
  c("group","L4_L5","L4_FullSynCom"),
  
  c("group","L5_L3","L5_FullSynCom"),
  c("group","L5_NB","L5_FullSynCom"),
  c("group","L5_L4_up","L5_FullSynCom"),
  c("group","L5_L5","L5_FullSynCom")
  
)

Res_contrasts <- NULL

for(comparison in list_comparisons){
  
  cat("Working on ",comparison,"\n")
  
  Res_contrasts <- dds %>% results(contrast = comparison) %>% 
    as.data.frame %>% 
    dplyr::mutate(gene_id = rownames(.),
                  Contrast =comparison[2:3] %>% paste0(collapse = "_vs_")) %>%
    dplyr::relocate(.data = .,gene_id) %>%
    remove_rownames(.data = .) %>%
    rbind(Res_contrasts,.)
  
}

saveRDS(object = Res_contrasts,file = "../cleandata/res_rnaseq_contrasts_syncom_do_all.RDS")
rm(list=ls())
gc()
