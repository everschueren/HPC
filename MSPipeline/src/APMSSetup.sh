#!/usr/bin/env bash
mkdir ./summary;
mkdir ./networks;
mkdir ./scripts;
mkdir ./data;
mkdir ./data/input;
>./data/input/$1_keys.txt
>./data/input/$1_remove.txt
>./data/input/$1_exclusions.txt
>./data/input/$1_collapse.txt
cp ~/Projects/HPCKrogan/Scripts/MSPipeline/tests/APMS_TEMPLATE.cfg ./data/$1.cfg
