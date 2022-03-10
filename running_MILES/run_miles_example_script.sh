
# this is an example script how miles was run with different hyperparameter settings on a hpc cluster with slurm
# needs to be adjusted for the dataset of interest

acc=usr
id=SNV_MILES_SUMMIT
pathSVM="./miles_publication"
path_data=$pathSVM/data_for_publication/MILES_input/SNV
infile=$path_data/"20220309_data_miles.csv"

path_out=$pathSVM/data_for_publication/MILES_results/SNV/
path_out=$path_out/miles_all/

mkdir -p $path_out



sigmas2=(50 100 500 1000 5000 10000 50000 100000 500000 1000000 10000000)
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
    echo "python $pathSVM/data_for_publication/running_MILES/MILES_cross_validation.py $s $c $infile $outfile"
    sbatch --account $acc -p Compute --job-name $id --output=$output --error=$error --mem=4G -n 4 -N 1 --wrap="
    python $pathSVM/data_for_publication/running_MILES/MILES_cross_validation.py $s $c $infile $outfile"

  done
done
