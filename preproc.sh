#!/bin/bash 

module load ROOT

root -b -q ./script/preprocess.C
root -b -q ./script/drawHists.C
root -b -q ./script/drawHists_e_slices.C
root -b -q ./script/plotIA.C
root -b -q ./script/plotCorr.C
