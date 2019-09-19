#!/bin/bash bash

for f in $(basename -s "_settings.txt" settings/*)
do
    prefix=$f
    suffix="_settings.txt"
    input_fn=settings/$prefix$suffix
    echo "$input_fn"
    ~/projects/treeducken/treeducken -i "$input_fn"
    mv *.tre ./data/$prefix/
done
