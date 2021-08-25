#!/bin/bash
input=$1
output=$2
num_cells=$3

cellbender remove-background --input $input --output $output --expected-cells $num_cells


