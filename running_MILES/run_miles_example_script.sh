
# this is an example SLURM script how miles was run with different hyperparameter settings on a hpc cluster with slurm
# here: SNVs + all ICB cohorts 

# THIS NEEDS TO BE ADJUSTED BY THE USER
acc=usr # account name
mem=10G # memory ; depends on no. data points in the input
partition=partition #partition name
id=SNV_MILES_SUMMIT # job name
pathSVM="./milneo_analysis" # main path
infile="path_to_input_file" # path to input file
path_out="path_to_output_file" # path to output file

mkdir -p $path_out

# hyperparameter sets
sigmas2=(50 100 500 1000 5000 10000 50000 100000 500000 1000000 5000000 10000000)
c_values=($(seq 0.1 0.1 1))


for s in "${sigmas2[@]}"
do
  for c in "${c_values[@]}"
  do
    echo $s"+"$c
    outfile=$path_out"/miles_"$s"_"$c".tsv"
    error=$path_out/"$s"_"$c".err
    output=$path_out/"$s"_"$c".out
    echo $infile
    echo $outfile

    sbatch --account $acc -p $partition --job-name $id --output=$output --error=$error --mem=$mem -n 1 -N 1 --wrap="
    python $pathSVM/running_MILES/MILES_cross_validation.py $s $c $infile $outfile"

  done
done
