#!/bin/bash
#methods=("no_decontamination" "soupx:autoEstCont" "soupx:background_genes" "soupx:top_background_genes" "decontx:no_cell_types" "decontx:with_cell_types" "cellbender")
methods=("no_decontamination" "soupx:autoEstCont" "decontx:no_cell_types" "cellbender")

for m in ${methods[@]}; do
	mkdir -p $1/$m/Rda/decontaminated_samples
	mkdir $1/$m/matrices
	mkdir $1/$m/plots
done
