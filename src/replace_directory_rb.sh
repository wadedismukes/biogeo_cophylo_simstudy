#!/bin/bash

# shellcheck disable=SC2164
cd rev-scripts/

all=$(ls)
for sub in $all
do

    cd "$sub"
    # shellcheck disable=SC2038
    find -- * | xargs -n 50 sed -i '' "s/data\//\/work\/LAS\/phylo-lab\/wade\/biogeo_cophylo_simstudy\/data\//g"
    cd ../
done
