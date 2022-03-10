# CHANGES:
# - remove patients with technical difficulties when running iCaM from combined dataset --> only combined approach will be re-run
# Hugo Pt 8 
# McDermott Pt164370145747331

library(tidyverse)
setwd("/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/code/miles_publication/")
source( "import_data.R")
date_today <- gsub("-", "", Sys.Date())

path_snv_out <- "data_for_publication/MILES_input/SNV/"
path_indel_out <- "data_for_publication/MILES_input/INDEL/"
path_fusion_out <- "data_for_publication/MILES_input/fusion_gene//"
path_combined_out <- "data_for_publication/MILES_input/combined/"


# Neoantigen features -----------------------------------------------------


#neoantigen features
path.to.feature.file <- "data_for_publication/20211112_interest_feature_names.csv"

features <- read_csv2(path.to.feature.file)$col.names
features.nice <- read_csv2(path.to.feature.file)$names_nice
directions <- read_csv2(path.to.feature.file)$direction


features_full <- c("Best_rank_MHCI_score", "Best_rank_MHCII_score", "rnaExpression", "Pathogensimiliarity_MHCI_9mer", "Pathogensimiliarity_MHCII",
                     "Selfsimilarity_MHCI", "Selfsimilarity_MHCII", "rnaVariantAlleleFrequency", "Amplitude_MHCII_rank", "Amplitude_MHCI_affinity",
                     "DAI_MHCI_affinity", "Dissimilarity_MHCI", "Generator_rate_MHCI", "Generator_rate_MHCII", "IEDB_Immunogenicity_MHCI", 
                     "IEDB_Immunogenicity_MHCII","MixMHC2pred_best_rank", "MixMHCpred_best_rank", "PHBR_I", "PHBR_II", "vaxrank_binding_score",
                     "Hex_alignment_score_MHCI", "Hex_alignment_score_MHCII", "Neoag_immunogenicity" , "Priority_score", "Recognition_Potential_MHCI_9mer",
                   "Tcell_predictor_score", "vaxrank_total_score", "PRIME_best_rank")
directions_full <- directions[which(features %in% features_full)]



# filter cohorts 
# cohorts: hugo, vanallen, riaz_pre, McDermott
cohorts_of_interest <- c("hugo", "vanallen", "riaz_pre")
names(cohorts_of_interest) <- cohorts_of_interest

# fusion genes
dat_fusions_transformed <- lapply(cohorts_of_interest, function(x) dat_fusions_transformed[[x]])


# indel
dat_indel_transformed <- lapply(cohorts_of_interest, function(x) dat_indel_transformed[[x]])

# snv
dat_transformed <- lapply(cohorts_of_interest, function(x) dat_transformed[[x]])






# TRANSFORMATION FUNCTION -------------------------------------------------


# bring features into the same directions
transform_feature_direction <- function(x) {

  res <- x %>%
    mutate(Best_rank_MHCI_score = 100 - Best_rank_MHCI_score)%>%
    mutate(Best_rank_MHCII_score = 100 - Best_rank_MHCII_score) %>%
    mutate(Dissimilarity_MHCI = 1- Dissimilarity_MHCI)%>%
    mutate(MixMHC2pred_best_rank = 100 - MixMHC2pred_best_rank)%>%
    mutate(MixMHCpred_best_rank = 100 - MixMHCpred_best_rank)%>%
    mutate(Selfsimilarity_MHCI = 1 - Selfsimilarity_MHCI)%>%
    mutate(Selfsimilarity_MHCII = 1 - Selfsimilarity_MHCII)%>%
    mutate(PRIME_best_rank = 100 - PRIME_best_rank)%>%
    mutate(PHBR_I = 100 - PHBR_I)%>%
    mutate(PHBR_II = 100 - PHBR_II) %>%
    mutate(rnaVariantAlleleFrequency = ifelse(rnaVariantAlleleFrequency == -1, NA, rnaVariantAlleleFrequency))
  return(res)
}

# substitute NA values by "worst" value
substitute.NA.feat <- function(df, features, feature.directions, response.nams ){
  df1 <- mapply(function(f, d){
    if(d == "<"){
      max.val <- min(df[,f], na.rm = T)
    } else if (d == ">"){
      max.val <- max(df[,f], na.rm = T)
    }
    df[is.na(df[,f]),f] <- max.val
    df[,f]
  }, f = features, d = feature.directions, SIMPLIFY = F)
  df2 <- do.call(cbind.data.frame, df1)
  res <- cbind.data.frame(df2, df[,response.nams], df[,"id"])
  colnames(res) <- c(colnames(df2), response.nams, "id")
  return(res)
}  


transform_data <- function(list, features, directions){
  dat <- lapply(list, function(x){
    cols <- c("binary_response", features, "patientIdentifier")
    x1 <- x %>%
      dplyr::select(cols) #%>%
    return(x1)
    
  })
  
  # format data
  dat_features <- bind_rows(dat, .id = "dataset")
  dat_features <- dat_features %>%
    mutate(id = paste0(dataset, "_", patientIdentifier)) %>%
    select(-dataset, -patientIdentifier) %>%
    relocate(id, .after = binary_response) %>%
    mutate(binary_response = ifelse(binary_response == "yes", 1, 0))
  
  # substitute NA values 
  dat_features_subst <- substitute.NA.feat(dat_features, features, directions, "binary_response")
  
  as_tibble(dat_features_subst)
  
}

# generate a list of dataframes with a respective feature for examination of feature importance
randomize_feature <- function(feature, dat){
  
  col_ind <- which(colnames(dat) == feature)
  dat_rem <- dat
  
  set.seed(42)
  l_ids <- list()
  l_ids[[1]] <- dat_rem[sample(nrow(dat_rem), replace = F),col_ind]
  for(i in 2:100){
    set.seed(i)
    l_ids[[i]] <- sample(l_ids[[i - 1]], replace = F)
  }
  
  dat_rem[,col_ind] <- l_ids[[100]]
  return(as_tibble(dat_rem))
}



# SNV -------------------------------------------------------------------

# transform feature direction
dat_features <- lapply(dat_transformed, transform_feature_direction)

# transform for miles
dat_features_full <- transform_data(list = dat_features, features = features_full, directions = directions_full)    


# randomise features
dat_snv_random <- lapply(features_full, randomize_feature, dat= dat_features_full)
names(dat_snv_random) <- features_full
path_snv_imp <- paste0(path_snv_out, "/", date_today, "_feature_importance_mel")
dir.create(path_snv_imp)

dat_snv_random %>%
  map2(names(dat_snv_random), function(x,y){
    nam <- paste0(path_snv_imp,"/", y, ".csv")
    write_csv(x, nam)
  })

# INDEL -------------------------------------------------------------------

# transform feature direction
dat_features_indel <- lapply(dat_indel_transformed, transform_feature_direction)

# transform for miles
dat_features_indel_full <- transform_data(list = dat_features_indel, features = features_full, directions = directions_full)    
dat_features_indel_full <- dat_features_indel_full %>%
  select(-Neoag_immunogenicity, -Tcell_predictor_score)

# randomise features
dat_indel_random <- lapply(features_full, randomize_feature, dat= dat_features_indel_full)
names(dat_indel_random) <- features_full
path_indel_imp <- paste0(path_indel_out, "/", date_today, "_feature_importance_mel")
dir.create(path_indel_imp)

dat_indel_random %>%
  map2(names(dat_indel_random), function(x,y){
    nam <- paste0(path_indel_imp,"/", y, ".csv")
    write_csv(x, nam)
  })


# fusion gene -------------------------------------------------------------

# transform feature direction
dat_features_fusions <- lapply(dat_fusions_transformed, transform_feature_direction)

# transform for miles
dat_features_fusions_full <- transform_data(list = dat_features_fusions, features = features_full, directions = directions_full)    
dat_features_fusions_full <- dat_features_fusions_full %>%
  select(-Neoag_immunogenicity, -Tcell_predictor_score)


# randomise features
dat_fusion_random <- lapply(features_full, randomize_feature, dat= dat_features_fusions_full)
names(dat_fusion_random) <- features_full
path_fusion_imp <- paste0(path_fusion_out, "/", date_today, "_feature_importance_mel")
dir.create(path_fusion_imp)

dat_fusion_random %>%
  map2(names(dat_fusion_random), function(x,y){
    nam <- paste0(path_fusion_imp,"/", y, ".csv")
    write_csv(x, nam)
  })



# combination snv, fusion genes, indels -----------------------------------

patients_icam_technical_issues <- c("mcdermott_Atezo_Pt164370145747331", "hugo_Pt8")

combined_full <- bind_rows(dat_features_full, dat_features_indel_full, dat_features_fusions_full)
combined_full <- combined_full %>%
  select(-Neoag_immunogenicity, -Tcell_predictor_score)%>%
  filter(! id %in% patients_icam_technical_issues)


# randomise features
dat_combined_random <- lapply(features_full, randomize_feature, dat= combined_full)
names(dat_combined_random) <- features_full
path_combined_imp <- paste0(path_combined_out, "/", date_today, "_randomise_feature_importance_melanoma")
dir.create(path_combined_imp)

dat_combined_random %>%
  map2(names(dat_combined_random), function(x,y){
    nam <- paste0(path_combined_imp,"/", y, ".csv")
    write_csv(x, nam)
  })

