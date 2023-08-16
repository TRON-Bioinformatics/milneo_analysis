
# About

This repository provides code and data related to the manuscript "Multiple instance learning to predict immune checkpoint blockade efficacy using neoantigen candidates" (https://doi.org/10.1101/2022.05.06.490587)

*  Example scripts to run miles can be found here [here](running_MILES).
    * These are slurm scripts to run MILES with different hyperparameter settings on a hpc cluster with slurm  
    * [MILES_cross_validation.py](running_MILES/MILES_cross_validation.py) is the python script that is called by the slurm scripts and runs MILES with nested cross validation  
*  Raw results from MILES can be found  [here](data_for_publication/MILES_results)
