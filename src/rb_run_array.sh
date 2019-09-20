#!/bin/bash

# first arg is rb (i.e. revbayes)
rb=$1
# second arg is rev-scripts/ where all the rev-scripts are
dir=$2

cd $dir

all=$(ls)

# go through the sub-directories of rev-scripts
for sub in $all
do
    cd $sub
    sbatch --array=0-9 --job-name=$sub /work/LAS/phylo-lab/wade/dec-cophylo/utility_scripts/array_sub.sh $rb $sub
    cd ../
done
