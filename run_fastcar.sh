#!/bin/bash
configs=("fastcar" "fastcar_no_mt" "fastcar_reclus")
for i in ${configs[@]}; do
	Rscript main.R $i
done
