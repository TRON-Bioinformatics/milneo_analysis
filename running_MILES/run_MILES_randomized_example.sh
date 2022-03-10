
# this is an example script how miles was run on a randomized dataset
# needs to be adjusted for the dataset of interest


acc=usr
id=RANDOM_SNV_MILES_SUMMIT
pathSVM="./miles_publication"
path_data=$pathSVM/data_for_publication/MILES_input/SNV
infile=$path_data/"20220310randomised_data.csv"

path_out=$pathSVM/data_for_publication/MILES_results/SNV/
path_out=$path_out/miles_all/


mkdir -p $path_out


# best hyperparameter set
# e.g.
s=10000
c=0.1

echo $s"+"$c
outfile=$path_out"/miles_"$s"_"$c".tsv"
error=$path_out/"$s"_"$c".err
output=$path_out/"$s"_"$c".out
echo $infile
echo $outfile
echo "python $pathSVM/code/training/MILES_cross_validation.py $s $c $infile $outfile"
sbatch --account $acc -p Compute --job-name $id  --output=$output --error=$error --mem=4G -n 5 -N 1 --wrap="
python $pathSVM/code/training/MILES_cross_validation.py $s $c $infile $outfile"
