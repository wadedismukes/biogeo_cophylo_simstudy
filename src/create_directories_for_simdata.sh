#!/usr/bin/env bash


mkdir -p settings/
mv *_settings.txt settings/
find ./settings/ -name "*_settings.txt" | xargs basename -s "_settings.txt" | xargs -I {} mkdir -p ./data/{}

