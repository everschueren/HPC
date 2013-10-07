#!/usr/bin/env bash
mkdir ./data;
mkdir ./data/input;
>./data/input/$1_keys.txt
>./data/input/$1_design.txt
>./data/input/$1_contrasts.txt
cp ~/Projects/HPCKrogan/Scripts/MSPipeline/tests/SILAC_TEMPLATE.cfg ./data/$1.cfg
