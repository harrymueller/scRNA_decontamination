#!/usr/bin/bash
CONFIG_PROFILE="stability_testing"
DIR="/data/Perkins/stability_testing"

if [ -z $1 ]; then
    echo "Number of subsets not supplied"
    exit 1
else
    N_SUBSETS=$1
fi    

# loop through n_subsets
for (( i=1; i<=$N_SUBSETS; i++ )); do
    # Rscript main.R $CONFIG_PROFILE $i

done

# rearrange file structure