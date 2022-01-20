#!/bin/bash
#set -e

echo "This script will setup the pipeline used to create results from the manuscript xxx"
echo 
echo "The analysis consists of three parts and this script will help you to set everything up to get going"
echo "Part 1: Phylogenomic analysis"
echo "Part 2: Gene calling and functional annotation"
echo "Part 3: Comparative genomic analysis"
echo
echo "The analyses performed in each part have to be run individually, look at the instructions in the README file."
echo "You may see git detached HEAD warning during this setup, these warnings can be ignored."
read -p "Would you like to setup the analysis now? (y/n)" proceed

if [[ $proceed == "y"  || $proceed == "Y" ]]
then
	echo "------------------------------------------------------------"
	echo "Will setup Part1 (phylogenomics)..."
	echo "Will clone phylociraptor into the ./phylogeny directory"
	git clone --recursive https://github.com/reslp/phylociraptor.git phylogeny/phylociraptor
	cd phylogeny/phylociraptor
	echo "Will checkout version used in the paper..."
	git checkout a93b4c8
	# add genome download here once we have NCBI accession numbers
	cd ../..
	echo "Will copy over config files for phylociraptor from input-data/phylociraptor."
	echo "These files contain settings and taxon selection used in the study."
	echo "HINT: Currently these files do not contain the final NCBI accession number as we are still waiting for them..."
	cp input-data/phylciraptor/config.yaml phylogeny/phylociraptor/data
	cp input-data/phylociraptor/LFS-phylogenomics-taxon-selection.csv phylogeny/phylociraptor/data	
	echo "Files copied successfully. You may now look at the README file to see how to recreate phylogenomic results."
	
	echo "------------------------------------------------------------"
	echo "Will setup Part2 (genome annotation)..."
	echo "Will clone smsi-funannotate into the genome-annotation directory"
	git clone --recursive https://github.com/reslp/smsi-funannotate.git genome-annotation/smsi-funannotate 
	cd genome-annotation/smsi-funannotate
	echo "Will checkout version used in the paper..."
	git checkout 1ed30e3
	cd ../..
	
else
	echo "OK, will abort setup."
 	exit 0
fi

