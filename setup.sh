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
echo
echo "IMPORTANT:"
echo "	- Installation time will depend on your internet connection speed because several large databases need to be downloaded"
echo "	- Make sure to have 150-200GB of free disk space available to download and store all software and databases."
echo "  - Make sure that you have followed the instructions in the README file prior to running this script!"
echo "  - Make sure that you have acquired a licensed copy of GenMark-ES 4.62 and Signal-P 4.1. Refer to the README file to see how this is done"
echo
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
	echo "Will download additional required software to input-data/genome-annotation/data/external"
	echo "Downloading Interproscan 5.48-83.0, this can take some time... please be patient."
	wget -P input-data/genome-annotation/data/external http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.48-83.0/interproscan-5.48-83.0-64-bit.tar.gz
	echo "Interproscan successfully downloaded. Unpacking..."
	tar -xvfz input-data/genome-annotation/data/external/interproscan-5.48-83.0-64-bit.tar.gz
	echo "Interproscan located in input-data/genome-annotation/data/external/interproscan-5.48-83.0"
	echo
	echo "Downloading eggnog-mapper databases, this can take some time... please be patient."
	echo "This command requires singularity to be installed!"
	mkdir genome-annotation/smsi-funannotate/data/eggnogdb
	singularity run docker://reslp/eggnog-mapper:1.0.3 download_eggnog_data.py NOG -y --data_dir genome-annotation/smsi-funannotate/data/eggnogdb
	echo "Eggnog database successfully downloaded"
	echo
	echo "Will setup GeneMark-ES. Make sure you have placed the files necessary into input-data."
	echo "Extracting tar archive to genome-annotation/smsi-funannotate/data/external/gm_et_linux_64"
	if [[ ! -f input-data/gmes_linux_64.tar.gz ]]
	then
		echo "input-data/gmes_linux_64.tar.gz file was not found, please make sure you have this file."
		echo "Check the README file for instructions on how to get it."
		exit 1
	fi	
	tar -xf input-data/gmes_linux_64.tar.gz -C genome-annotation/smsi-funannotate/data/external/
	mv genome-annotation/smsi-funannotate/data/external/gmes_linux_64/ genome-annotation/smsi-funannotate/data/external/gm_et_linux_64
	if [[ ! -f input-data/gm_key_64.gz ]]
	then
		echo "input-data/gm_key_64.gz file was not found, please make sure you have this file."
		echo "Check the README file for instructions on how to get it."
		exit 1
	fi
	gunzip input-data/gm_key_64.gz
	cp input-data/gm_key_64 genome-annotation/smsi-funannotate/.gm_key
	echo "GeneMark setup complete"
	echo
	echo "Will setup Signal-P now. Make sure you have placed the files necessary into input-data."
	echo "Extracting tar archive to genome-annotation/smsi-funannotate/data/external/signalp-4.1"
	if [[ ! -f input-data/signalp-4.1g.Linux.tar.gz ]]
	then
		echo "input-data/signalp-4.1g.Linux.tar.gz was not found, please make sure you have this file."
		echo "Check the README file for instructions on how to get it."
		exit 1
	fi	
	tar -xf input-data/signalp-4.1g.Linux.tar.gz -C genome-annotation/smsi-funannotate/data/external/
	echo "Will change signalp script now so that paths are correct..."
	sed -i 's$/usr/opt/www/pub/CBS/services/SignalP-4.1/signalp-4.1$/data/external/signalp-4.1$' genome-annotation/smsi-funannotate/data/external/signalp-4.1/signalp
	sed -i 's$/var/tmp$/tmp$' genome-annotation/smsi-funannotate/data/external/signalp-4.1/signalp
	sed -i 's$MAX_ALLOWED_ENTRIES=20000$MAX_ALLOWED_ENTRIES=100000$' genome-annotation/smsi-funannotate/data/external/signalp-4.1/signalp
	echo "Signal-P setup complete"
	echo

else
	echo "OK, will abort setup."
 	exit 0
fi

