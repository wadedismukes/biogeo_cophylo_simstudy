#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=14-00:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node 
#SBATCH --mem=32G   # maximum memory per node
#SBATCH --job-name="run-rev-all"
#SBATCH --mail-user=waded@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load revbayes
rb=$1
sub=$2

$rb run_epoch_$(echo $sub)_$(echo $SLURM_ARRAY_TASK_ID).Rev 
