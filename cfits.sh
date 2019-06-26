#!/bin/bash
mkdir -p cfit_outputs
mkdir -p cfit_cfgs
rm ./cfit_outputs/*
rm ./cfit_cfgs/*
module load ROOT
python script/gencfgs.py

for f in `ls ./cfit_cfgs/`
do
	root -b -q script/combinedFit.C < ./cfit_cfgs/$f 2>&1 | tee ./cfit_outputs/$f
done

root -b -q script/parse_cfit_outputs.C
echo cfit_summary.root | root -b -q script/plot_cfit_summary.C
