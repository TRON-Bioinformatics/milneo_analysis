
# About

This repository provides code and data related to the manuscript "Multiple instance learning to predict immune checkpoint blockade efficacy using neoantigen candidates" (https://doi.org/10.1101/2022.05.06.490587)

R markdown scripts to generate the figures of the manuscript can be found [here](data_analysis/Generate_Manuscript_Figures).  

Furthermore, this repository contains raw, intermediate and downstream data that was generated in this work and related code:  

*  Detected neoantigen candidates from SNVs, INDELs and fusion genes for 5 ICB cohorts can be found [here](data_for_publication/raw_data)
*  All R-scripts related to data preparation for MILES can be found in the folder [data_preparation](data_preparation)
    * R-scripts for each tumor entity to generate datasets in the context of the mutation type to run MILES are [here](data_preparation). The suffix "wo_SD" indicates that patients with stable disease are not considered in this analysis.
    * R-scripts  for each tumor entity to generate data for feature importance analysis are [here](data_preparation/feature_importance)
    * R-scripts  for each tumor entity to generate randomized data are [here](data_preparation/randomized_candidates)
*  The MILES input data as it was used in this study can be found  [here](data_for_publication/MILES_input). There is a folder for each mutation type (here: combined = all mutation types). E.g. :
    * data_miles.csv : all patients  
    * data_without_SD.csv : all patients without stable disease  
    * data_miles_mel.csv : all melanoma patients  
    * data_miles_mel_wo_SD : all melanoma patients without stable disease  
*  Example scripts to run miles can be found here [here](data_for_publication/running_MILES).
    * The paths need to be adjusted for the input of interest (e.g. from [data_for_publication/MILES_input](data_for_publication/MILES_input))  
    * These are slurm scripts to run MILES with different hyperparameter settings on a hpc cluster with slurm  
    * [MILES_cross_validation.py](running_MILES\MILES_cross_validation.py) is the python script that is called by the slurm scripts and runs MILES with nested cross validation  
*  Raw results from MILES can be found  [here](data_for_publication/MILES_results)
*  Downstream analysis scripts to generate the downstream data used for the figures in the manuscript can be found [here](data_analysis/Downstream_analysis)
    * Summarize_MILES_x.Rmd : summarises the performance results from MILES.
    * Summarize_neoantigen_candidate_load_x.RmD : Analysis of the neoantigen candidate load
    * ROC_analysis_neoantigen_candidate_load.RmD : ROC-curve analysis of predicting ICB efficacy based on neoantigen candidate load
