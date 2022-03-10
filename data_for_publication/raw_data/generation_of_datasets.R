
# generation of dataset files for publication on figshare that can be used by others to reproduce the results 

source("/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/code/helpers/config_cohort_files_v3.R")
source("/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/code/helpers/config_cohort_transformation.R")

out_path = "/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/data_for_publication/"

# SNVs
lapply(names(dat_transformed), function(x){
  f <- paste0(out_path, "/neoantigen_candidates_snv_", x, ".tsv")
  write_delim(dat_transformed[[x]] , f, delim = "\t" )
})

# INDELs
lapply(names(dat_indel_transformed), function(x){
  f <- paste0(out_path, "/neoantigen_candidates_indel_", x, ".tsv")
  write_delim(dat_indel_transformed[[x]] , f, delim = "\t" )
})

# Fusion genes 
lapply(names(dat_fusions_transformed), function(x){
  f <- paste0(out_path, "/neoantigen_candidates_fusion_", x, ".tsv")
  write_delim(dat_fusions_transformed[[x]] , f, delim = "\t" )
})
