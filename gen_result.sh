#!/bin/bash
version=release_noeff
for d in hists_archive/release_candidate/*/; do
    if [[ $d == *nGd* ]]; then
        continue
    elif [[ $d == *nH* ]]; then
        continue
    elif [[ $d == *unified* ]]; then
        analysis="unified"
    fi

    for h in ${d}*/; do
        suffix=$(basename $h)
        echo 'Processing' $suffix
        rm -r ./hists
        cp -r $h ./hists
        ./compile.sh
        rm plots/cfits/*.png
        rm plots/cfits_slice/*.png
        if [[ $d == *nGd* ]]; then
            analysis="nGd"
        elif [[ $d == *nH* ]]; then
            analysis="nH"
        elif [[ $d == *unified* ]]; then
            analysis="unified"
        fi
        echo $h $analysis
        dir="./results/$version/$analysis/${suffix}_fix"
        echo $dir
        mkdir -p $dir
        #./combinedFit fix $analysis | tee $dir/fit_result 
        ./combinedFit fix $analysis &> $dir/fit_result 

        cp -r plots/cfits $dir/
        cp -r plots/cfits_slice $dir/

        rm plots/cfits/*.png
        rm plots/cfits_slice/*.png

        dir="./results/$version/$analysis/${suffix}_free"
        echo $dir
        mkdir -p $dir
        #./combinedFit free $analysis | tee $dir/fit_result 
        ./combinedFit free $analysis &> $dir/fit_result

        cp -r plots/cfits $dir/
        cp -r plots/cfits_slice $dir/       
    done
done

