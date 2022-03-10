# this script imports the raw data 

library(tidyverse)

# ADJUST THIS PATH IF NECESSARY
path_data <- "data_for_publication/raw_data/"

cohorts <- c("hugo" , "vanallen" ,"riaz_pre" ,"miao", "mcdermott_Atezo")
#mut_class <- c("snv", "indel", "fusion")

# SNVs
dat_transformed <- lapply(cohorts, function(x){
  f <- paste0(path_data,"neoantigen_candidates_snv_", x, ".tsv")
  read_delim(f, delim = "\t")
})

names(dat_transformed) <- cohorts

# INDELs
dat_indel_transformed<- lapply(cohorts, function(x){
  f <- paste0(path_data,"neoantigen_candidates_indel_", x, ".tsv")
  read_delim(f, delim = "\t")
})

names(dat_indel_transformed) <- cohorts

# INDELs
dat_fusions_transformed <- lapply(cohorts, function(x){
  f <- paste0(path_data,"neoantigen_candidates_fusion_", x, ".tsv")
  read_delim(f, delim = "\t")
})

names(dat_fusions_transformed) <- cohorts
