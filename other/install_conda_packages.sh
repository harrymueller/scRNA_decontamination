#!/bin/bash
# Getting input
[ -z "$1" ] && opt=-h || opt=$1

############################################################
# R packages
############################################################
if [ $opt = -r ] || [ $opt = --rPackages ]; then
	# Basic packages
	echo "r-base & r-essentials"
	conda install -c r r-base r-essentials -y
	echo "varhandle & config"
	conda install -c conda-forge r-varhandle r-config r-spatstat=1.64_01 -y

	# Seurat
	echo "seurat"
	conda install -c conda-forge r-seurat=3.2.3 -y

	# soupx
	echo "soupx"
	conda install -c conda-forge r-soupx -y

	# decontx via celda & biocmanager
	echo "biocmanager"
	conda install -c conda-forge r-biocmanager -y
	echo "celda"
	Rscript -e 'BiocManager::install("celda")'

	# packages required for analysis
	echo "ggplot2, xlsx, reshape2"
	conda install -c r r-ggplot2 r-xlsx r-reshape2 -y
	echo "plotly, cowplot, patchwork"
	conda install -c conda-forge r-plotly r-cowplot r-patchwork -y
	echo "orca"
	conda install -c plotly plotly-orca -y

	# others
	conda install -c conda-forge r-htmlwidgets r-vctrs -y
############################################################
# CellBender (CPU)
############################################################
elif [ $opt = -c ] || [ $opt = --cellbenderCPU ]; then
	echo "Ensure that python 3.7 is installed"
	echo ""
	echo "pytables, pytorch (CPU)"
	conda install -c anaconda pytables -y
	conda install pytorch torchvision -c pytorch -y

	echo "cellbender"
	cd $2
	git clone https://github.com/broadinstitute/CellBender.git
	pip install -e CellBender
############################################################
# CellBender (GPU)
############################################################
elif [ $opt = -g ] || [ $opt = --cellbenderGPU ]; then
	echo "Ensure that python3.7 installed"
	echo ""
	echo "pytables, pytorch (CPU)"
	conda install -c anaconda pytables -y
	conda install pytorch torchvision torchaudio cudatoolkit=10.2 -c pytorch -y

	echo "cellbender"
	cd $2
	git clone https://github.com/broadinstitute/CellBender.git
	pip install -e CellBender
############################################################
# Options
############################################################
else
	echo "Options:"
	echo "  -r || --rPackages     > install R packages in current conda environment"
	echo "  -c || --cellbenderCPU > install cellbender (CPU version) in the current conda environment"
	echo "  -g || --cellbenderGPU > install cellbender (GPU version) in the current conda environment"
	echo ""
	echo "  <dir>                 > directory to install Celbender in"
fi
