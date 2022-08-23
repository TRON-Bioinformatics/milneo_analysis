
# this is an example script how miles was run on a randomized dataset
# here: SNVs + all ICB cohorts

# THIS NEEDS TO BE ADJUSTED BY THE USER
acc=usr # account name
partition=partition #partition name
id=RANDOM_SNV_MILES_SUMMIT # job name
pathSVM="./miles_publication" # main path
path_data=$pathSVM/data_for_publication/MILES_input/SNV # path to input data
infile=$path_data/"20220310randomised_data.csv" # path to randomized data
path_out=$pathSVM/data_for_publication/MILES_results/SNV/ # path to output
path_out=$path_out/miles_all/ # path to output


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
echo "python $pathSVM/data_for_publication/running_MILES/MILES_cross_validation.py $s $c $infile $outfile"
sbatch --account $acc -p $partition --job-name $id  --output=$output --error=$error --mem=4G -n 5 -N 1 --wrap="
python $pathSVM/data_for_publication/running_MILES/MILES_cross_validation.py $s $c $infile $outfile"
