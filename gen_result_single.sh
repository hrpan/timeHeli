#!/bin/bash
./compile.sh
rm plots/cfits/*.png
rm plots/cfits_slice/*.png
echo $h $analysis
dir="$1_fix"
analysis=$2
echo $dir
mkdir -p $dir
./combinedFit fix $analysis | tee $dir/fit_result 

cp -r plots/cfits $dir/
cp -r plots/cfits_slice $dir/

rm plots/cfits/*.png
rm plots/cfits_slice/*.png

dir="$1_free"
echo $dir
mkdir -p $dir
./combinedFit free $analysis | tee $dir/fit_result 

cp -r plots/cfits $dir/
cp -r plots/cfits_slice $dir/       

