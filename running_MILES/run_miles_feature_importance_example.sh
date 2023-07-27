
# this is an example SLURM script  how feature importance analysis was run with best hyperparameter set
# here: SNVs + all ICB cohorts

# THIS NEEDS TO BE ADJUSTED BY THE USER
acc=usr # account name
partition=partition #partition name
mem=10G # memory ; depends on no. data points in the input
id=FEATURE_IMPORTANCE_SNV_SUMMIT # job name

pathSVM="./milneo_analysis" # main path

path_data="path_to_input" # path to input data -->  contains a folder for each feature. Each feature folder contains n files with permutaed feature of interest

outpath="path_to_output" # output path 
mkdir -p $outpath

# choose best hyperparameter set
# e.g.
sigmas2=10000
c_values=0.1


# train models with randomized features
for inpath in $path_data/*
do
  feature=${inpath/"$path_data"/}

  for infile in $inpath/*.csv
  do
    
    name="${infile/"$inpath"/}"

    outpath_feature=$outpath/$feature
	
	if [[ -e $infile ]]; then
		
		mkdir -p $outpath_feature
		
		outfile=$outpath_feature/$name.tsv
		error=$outpath_feature/$name.err
		
		output=$outpath_feature/$name.out		
		
		sbatch --account $acc -p $partition --job-name $id --output=$output --error=$error --mem=$mem -n 1 -N 1 --wrap="module load anaconda/3/2019;
		python $pathSVM/running_MILES/MILES_cross_validation.py $s $c $infile $outfile"
		
		sleep 0.75
		
	fi 
     
  done
done
