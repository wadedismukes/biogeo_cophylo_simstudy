#!/bin/bash

cd rev-scripts/

all=$(ls)
for sub in $all
do
    cd $sub
    find *.Rev | xargs sed -i "s/data\//\/work\/LAS\/phylo-lab\/wade\/dec-cophylo\/data\//g"
    cd ../
done
