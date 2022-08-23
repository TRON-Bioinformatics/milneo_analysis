
# this is an example SLURM script  how feature importance analysis was run with best hyperparameter set
# here: SNVs + all ICB cohorts

# THIS NEEDS TO BE ADJUSTED BY THE USER
acc=usr # account name
partition=partition #partition name
id=FEATURE_IMPORTANCE_SNV_SUMMIT # job name
pathSVM="./miles_publication" # main path
path_infile=$pathSVM/data_for_publication/MILES_input/SNV # path to input data
path_data=$path_infile/20220310feature_importance/ # path to input data -->  files with feature of interest permutated
path_results=$pathSVM/data_for_publication/MILES_results/SNV/ # path to output
outpath=$path_results/feature_importance_all/ # path to output

mkdir -p $outpath

# choose best hyperparameter set
# e.g.
sigmas2=10000
c_values=0.1


# train models with removed featurese

for infile in $path_data/*.csv
do
  echo $sigmas2"+"$c_values
  echo $infile
  name="${infile/$path_data"/"/}"
  name=${name/".csv"/}
  #name="cohorts_full"
  outfile=$outpath$name".csv"
  error=$outpath/"$name".err
  output=$outpath/"$name".out
  echo $outfile
  echo $infile
  echo "python $pathSVM/data_for_publication/running_MILES/MILES_cross_validation.py $sigmas2 $c_values $infile $outfile"
  sbatch --account $acc -p $partition --job-name $id  --output=$output --error=$error --mem=4G -n 5 -N 1 --wrap="
  python $pathSVM/data_for_publication/running_MILES/MILES_cross_validation.py $sigmas2 $c_values $infile $outfile"

done
