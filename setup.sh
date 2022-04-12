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
	echo "Phylociraptor has changed and improved quite a bit since it was used to analyse the data for this paper."
	echo "There is no harm to use the most current version of it. However it is also possible to use the same version as in the paper:"
	read -p "Would you like to use the exact version of phylociraptor as was used in the paper? (y/n)" phyv
	if [[ $phyv == "y" || $phyv == "Y" ]]
	then
		echo "Will checkout version (git commit a93b4c8) used in the paper now..."
		git checkout a93b4c8
	else
		echo "Will keep updated version of phylociraptor"
	fi
	cd ../..
	echo "Will copy over config files for phylociraptor from input-data/phylociraptor."
	echo "These files contain settings and taxon selection used in the study."
	echo "Copy config file"
	cp input-data/phylociraptor/config.yaml phylogeny/phylociraptor/data/
	echo "Copy samples file"
	cp input-data/phylociraptor/LFS-phylogenomics-taxon-selection.csv phylogeny/phylociraptor/data	
	echo "Copy genome_download script"
	cp input-data/phylociraptor/download_assemblies.sh phylogeny/phylociraptor/data
	echo "Copy Dibaes and Graphis assemblies"
	cp input-data/phylociraptor/assemblies/* phylogeny/phylociraptor/data/assemblies/
	echo "Now I will download the remaining assemblies from NCBI"
	cd phylogeny/phylociraptor/data/
	bash download_assemblies.sh
	cd ../../../
	echo "$(pwd)"
	echo " Setup for part 1 (phylogenomics) is complete. You may now look at the README file to see how to recreate phylogenomic results."
	
	echo "------------------------------------------------------------"
	echo "Will setup Part2 (genome annotation)..."
	echo "Will clone smsi-funannotate into the genome-annotation directory"
	git clone --recursive https://github.com/reslp/smsi-funannotate.git genome-annotation/smsi-funannotate 
	cd genome-annotation/smsi-funannotate
	echo "Will checkout version used in the paper..."
	git checkout 1ed30e3
	cd ../..
	echo "Will download additional required software to input-data/genome-annotation/data/external"
	echo "Downloading Interproscan 5.48-83.0, this can take some time... please be patient the."
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
	echo "Now checkiong for RepeatMasker libraries in input-files/Repeatmasker_Libraries.tar.gz"
	echo "If this file is not present, the annotation pipeline will be configured to used TanTan for repeatmasking instead!"
	echo "Please refer to the README file if you would like to know more about this..."
	if [[ -f input-files/Repeatmasker_Libraries.tar.gz ]]
	then
		echo "Repeatmasker Library found, will unpack now..."
		tar -xf input-files/Repeatmasker_Libraries.tar.gz -C genome-annotation/smsi-funannotate/data/
		echo "Repeatmasker Library was unpacked to genome-annotation/smsi-funannotate/data/."
		echo "IMPORTANT: Make sure that the folder of the library is called RepeatMaskerLibraries!"		
	else
		echo "Repeatmasker Library not found, will change Repeatmasking to tantan instead."
		sed -i 's/repeatmasker/tantan/' genome-annotation/smsi-funannotate/data/config.yaml
	fi
        echo "Copy smsi-funannotate config files to the correct directory:"
        cp input-data/genome-annotation/config.yaml genome-annotation/smsi-funannotate/data/
        cp input-data/genome-annotation/data.csv genome-annotation/smsi-funannotate/data/
        echo "Copy assemblies to correct directory:"
        cp phylogeny/phylociraptor/data/assemblies/*.f* genome-annotation/smsi-funannotate/data/assemblies/
        echo "Setup part 2 is done"	

	echo "------------------------------------------------------------"
	echo "Will setup Part3 (comparative genomics)..."
	echo "Downloading data from Dryad (https://datadryad.org/stash/dataset/doi:10.5061/dryad.3xsj3txjb)"
	cd comparative-analysis/data
	wget https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.3xsj3txjb/download -O dryad_dataset.zip
	echo "Download finished. Unpacking data."
	unzip dryad_dataset.zip
	mv README ../README-Dryad
	unzip LFS-Cazy-comparative-phylogeny-annotations-settings.zip
	mv for_dryad/genome_annotations/* 83_genomes/
	mv for_dryad/other/all_ascomycota_peroxidases_redoxibase.fas all_ascomycota_peroxidases_redoxibase.fas
	mv for_dryad/other/sugar_transporters_characterized_reference_seqs.fa sugar_transporters_characterized_reference_seqs.fa
	mv for_dryad/other/stats_genomes.csv 83_genomes/
	mv for_dryad/other/lifestyle 83_genomes/
	mv for_dryad/other/r8s_template.nex 83_genomes/
	mv for_dryad/other/character_information.csv 83_genomes/
#	mkdir settings_and_parameters
#	mv for_dryad/other/apriori_cazyme_sets.txt settings_and_parameters/
	mv for_dryad/other/cloned_putative_cellulases.txt settings_and_parameters/
	mv for_dryad/other/color_information.csv settings_and_parameters/
	mv for_dryad/other/taxonomy_information.csv settings_and_parameters/
	mv for_dryad/phylogeny 83_genomes/
	echo "Setup is finished. Please refer to the README files for further instructions on how to perform analyses."
	
else
	echo "OK, will abort setup."
 	exit 0
fi	

