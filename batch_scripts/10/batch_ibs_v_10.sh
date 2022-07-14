#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --array=1-80
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=ibs_vstm
#SBATCH --mail-type=END
#SBATCH --mail-user=xl1005@nyu.edu
#SBATCH --output=ibs2_%j.out

PROJECT_FOLDER="IBS-2-development"

#model=psycho
model=vstm
#model=fourinarow

proc_id=${SLURM_ARRAY_TASK_ID}

#var1='ibs_alloc'
#SPACE=' '
#var2='50'
#$method = "$var1${SPACE}$var2"
#method=ibs_3
#Nsamples=10
alpha=1
method=ibs_10
#method=fixed
#method=fixed 
#method=fixedb
#method=exact


if [ $method = "exact" ]; then
    workdir=$SCRATCH/${PROJECT_FOLDER}/results/${model}/${method}
else
    workdir=$SCRATCH/${PROJECT_FOLDER}/results/${model}/${method}${Nsamples}
fi

module purge. 
module load matlab/2020b
export MATLABPATH=$HOME/${PROJECT_FOLDER}/matlab
mkdir $SCRATCH/${PROJECT_FOLDER}/results
mkdir $SCRATCH/${PROJECT_FOLDER}/results/${model}
mkdir $workdir
cd $workdir

echo $model $method $Nsamples $proc_id

echo "addpath('$SCRATCH/${PROJECT_FOLDER}/matlab/'); recover_theta('${model}','${method}','${alpha}',${proc_id},${Nsamples}); exit;" 
cat<<EOF | matlab -nodisplay
%job_id = str2num(strjoin(regexp('$proc_id','\d','match'), ''))
job_id = str2num('$proc_id')
alpha = str2num('$alpha')
recover_theta('vstm','ibs_10',alpha,job_id)

EOF

