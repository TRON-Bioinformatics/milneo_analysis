
# this is an example script how miles was run on a randomized dataset
# here: SNVs + all ICB cohorts

# THIS NEEDS TO BE ADJUSTED BY THE USER
acc=usr # account name
partition=partition #partition name
id=RANDOM_SNV_MILES_SUMMIT # job name
mem=10G # memory ; depends on no. data points in the input

pathSVM="./milneo_analysis" # main path
infile="path_randomized_input" # path to randomized input data
path_out="path_to_out" # path to output

mkdir -p $path_out


# best hyperparameter set
# e.g.
s=10000
c=0.1

outfile=$path_out"/miles_"$s"_"$c".tsv"
error=$path_out/"$s"_"$c".err
output=$path_out/"$s"_"$c".out

sbatch --account $acc -p $partition --job-name $id  --output=$output --error=$error --mem=$mem -n 1 -N 1 --wrap="
python $pathSVM/running_MILES/MILES_cross_validation.py $s $c $infile $outfile"
