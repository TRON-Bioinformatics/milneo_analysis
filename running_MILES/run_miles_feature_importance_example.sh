
#best conditions for combined: all predictd neoantigen candidates and filtered features
# use data 20211006


acc=usr
id=FEATURE_IMPORTANCE_SNV_SUMMIT

pathSVM="./miles_publication"
path_infile=$pathSVM/data_for_publication/MILES_input/SNV
path_data=$path_infile/20220310feature_importance/

path_results=$pathSVM/data_for_publication/MILES_results/SNV/
outpath=$path_results/feature_importance_all/

mkdir -p $outpath

# best hyperparameter set
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
  echo "python $pathSVM/code/training/MILES_cross_validation.py $sigmas2 $c_values $infile $outfile"
  sbatch --account $acc -p Compute --job-name $id --time=120:00:00 --output=$output --error=$error --mem=4G -n 5 -N 1 --wrap="
  python $pathSVM/code/training/MILES_cross_validation.py $sigmas2 $c_values $infile $outfile"

done
