rule scrape_cazy:
	input:
		cazy = config["funannotate_input"]["cazy"]
	output:
		checkpoint = "results/checkpoints/prepare_scrape_cazy.done",
		all_cazy_info = "results/cazy_characterization/cazy_information/all_cazy_data.csv",
		interesting_cazy = "results/cazy_characterization/cazy_information/all_interesting_cazy.txt"
	singularity:
		"docker://reslp/scrape_cazy:1"
	params:
		prefix = config["prefix"],
		wd = os.getcwd(),
		search_terms = config["scrape_cazy"]["terms"]
	shadow: "shallow"
	shell:
		"""
		#download all information on characterized cazymes from cazy.org
		scrape_cazy.py -f $(cat {params.wd}/{input.cazy} | awk -F "," 'NR > 1 {{print $1;}}' | awk -F "_" '{{print $1}}' | uniq | tr "\\n" "," | sed '$s/,$//')
		echo "Combining files..."
		tail -n +2 -q *_characterized.txt > {output.all_cazy_info}
		# strangely this only works if the terms specified in the config file exist. It works outside of snakemake. Maybe it some problem with the different shell 
		# invoked by snakemake...
		names="{params.search_terms}"
		echo "Searching families: "$names
		if [[ -f {output.interesting_cazy} ]]; then
			rm {output.interesting_cazy}
			touch {output.interesting_cazy}
		else
			touch {output.interesting_cazy}
		fi
		for name in $names; do
			echo "Working on: "$name
			if grep -q "^$name" all_cazy_data.csv; then # this is needed because grep returns non-zero exicode when no match is found
				cat all_cazy_data.csv | grep "^$name" | awk '{{print $1}}' | uniq >> {output.interesting_cazy}
			else 
				echo "Not found: "$name
			fi	
		done	

		touch {output.checkpoint}
		"""

rule extract_cazy_proteins:
	input:
		combined_proteins = config["funannotate_input"]["protein_folder"],
		combined_gff = config["funannotate_input"]["annotations_gff_folder"]
	output:
		cazy_proteins = "results/cazy_characterization/saccharis/all_cazy_proteins.fas",
		checkpoint = "results/checkpoints/extract_cazy_proteins.done"
	singularity:
		"docker://reslp/biopython_plus:1.77"
	shell:
		"""
		#cat {input.combined_proteins}/*.fa >> {output.cazy_proteins}
		bin/select_cazy_from_gff.py -gff {input.combined_gff} -fasta {input.combined_proteins} -item CAZy >> {output.cazy_proteins}
		touch {output.checkpoint}
		"""

# all rules using the saccharis container need the following bindpoints: "-B $(pwd):/data -B $(pwd)/bin:/usr/local/external"
rule prepare_saccharis:
	input:
		combined_proteins = rules.extract_cazy_proteins.output.cazy_proteins,
		cazy = config["funannotate_input"]["cazy"]	
	output:
		renamed_sequences = "results/cazy_characterization/saccharis/renamed_proteins.fas",
		checkpoint = "results/checkpoints/prepare_saccharis.done"
	singularity:
		"docker://reslp/saccharis:1"
	shell:
		"""
		perl -pe 's/\>/$& . U . sprintf("%08d", ++$n) . " "/ge' {input.combined_proteins} > {output.renamed_sequences}
		touch {output.checkpoint}
		"""	
		
rule saccharis:
	input:
		seqs = rules.prepare_saccharis.output.renamed_sequences,
		cazy = rules.scrape_cazy.output.checkpoint,
		apriori_sets = config["other_input"]["apriori_cazy_sets"]
	output:
		checkpoint = "results/checkpoints/saccharis.done"
	singularity:
		"docker://reslp/saccharis:1"
	params:
		saccharis_threads = config["saccharis"]["threads"],
		parallel_jobs = config["saccharis"]["parallel_jobs"],
		prefix = config["prefix"],
		wd = os.getcwd()
	shell:
		"""
		export TMPDIR={params.wd}/tmp
		# Saccharis sometimes does not finsish (reasons unkown, maybe memory), but it still returns a valid exit code 0, like it finished normally.
		# Therefore parallel does not catch the failed jobs. 
		# The only way around this seems to rerun the saccharis rule
		# However it is only necessary to run it for the failed familes.
		# The for loop will check for which families saccharis finished correctly and then rerun the failed ones.
		# It will run for all cazymes if nothing has been run before.
		rm -f {params.wd}/results/cazy_characterization/saccharis/interesting_cazy_families.txt
		for cazy in $(echo $(cat {input.apriori_sets} | tail -n +2 | awk -F "," '{{print $1}}' | awk -F "_" '{{print $1}}') $(cat results/cazy_characterization/cazy_information/all_interesting_cazy_families.txt | sort | uniq | tr '\\n' ' ') | sort | uniq); 
		do 
			if ! grep -Fxq "Cazy Pipeline Finished" {params.wd}/results/cazy_characterization/saccharis/$cazy.log;
			then 
				echo $cazy >> {params.wd}/results/cazy_characterization/saccharis/interesting_cazy_families.txt
				Saccharis.pl -d /data/results/cazy_characterization/saccharis -g characterized -s /data/{input.seqs} -t {params.saccharis_threads} -f $cazy &> {params.wd}/results/cazy_characterization/saccharis/$cazy.log
				rm -rf {params.wd}/results/cazy_characterization/saccharis/$cazy
			else
				continue
			fi
		done
		echo $(cat /data/results/cazy_characterization/saccharis/interesting_cazy_families.txt | tr '\\n' ' ')
		#parallel --joblog {params.wd}/results/cazy_characterization/saccharis/parallel_logfile --results $TMPDIR  -j {params.parallel_jobs} Saccharis.pl -d /data/results/cazy_characterization/saccharis -g characterized -s /data/{input.seqs} -t {params.saccharis_threads} -f {{}} "&>" {params.wd}/results/cazy_characterization/saccharis/{{}}.log ::: $(cat /data/results/cazy_characterization/saccharis/interesting_cazy_families.txt | tr '\\n' ' ')
		touch {output.checkpoint}
		"""

rule get_saccharis_mapping_data:
	input:
		seqs = rules.prepare_saccharis.output.renamed_sequences,
		genome_info = config["stats_genomes"],
		discrete_data = config["other_input"]["character_information"]
	output:
		saccharis_mapping_data = "results/cazy_characterization/saccharis_plotting/saccharis_mapping_information.tsv",
		checkpoint = "results/checkpoints/get_saccharis_mapping_data.done"
	conda:
		"../envs/pyutils.yml"
	shell:
		"""
		bin/get_saccharis_mapping_data.py -fasta {input.seqs} -info {input.genome_info} -map {input.discrete_data} > {output.saccharis_mapping_data}
		touch {output.checkpoint}
		"""

# deeploc needs the following bindpoints: "-B $(pwd)/external/deeploc-1.0/bin/:/external -B $(pwd)/external/deeploc-1.0/DeepLoc:/usr/lib/python3/dist-packages/DeepLoc -B $(pwd):/data"
#
rule deeploc:
	input:	
		checkpoint = rules.saccharis.output.checkpoint
	output:
		checkpoint = "results/checkpoints/deeploc.done",
		outdir = directory("results/cazy_characterization/deeploc/")
	params:
		prefix = config["prefix"],
		wd = os.getcwd()
	singularity:
		"docker://reslp/deeploc:1.0"
	shell:
		"""
		cd /data/{output.outdir}
		for seqfile in $(find /data/results/cazy_characterization/saccharis/*/characterized/dbcan/*extracted.fasta -not -empty); do
			echo $seqfile
			outname=$(basename "$seqfile" .fasta)
			echo $outname
			deeploc -f $seqfile -o $outname 
		done;
		cd /data
		touch {output.checkpoint}
		"""


rule plot_saccharis_trees:
	input:
		saccharis_mapping_data = rules.get_saccharis_mapping_data.output.saccharis_mapping_data,
		checkpoint = rules.scrape_cazy.output.checkpoint,
		checkpoint2 = rules.deeploc.output.checkpoint
	output:
		checkpoint = "results/checkpoints/plot_saccharis_trees.done"
	params:
		wd = os.getcwd(),
		prefix = config["prefix"]
	conda:
		"../envs/rreroot.yml"
	shell:
		"""
		for treename in $(ls results/cazy_characterization/saccharis/*/characterized/*.tree); do	
			echo $treename
			#IFS="_"
			#read -ar fields <<< $(basename $treename)
			famname=$(basename $treename _characterized.tree)
			#famname="${{fields[0]}}"
			echo $famname 
		        Rscript bin/annotate_saccharis_trees.R {params.wd} $treename results/cazy_characterization/cazy_information/"$famname"_characterized.txt {input.saccharis_mapping_data} "results/cazy_characterization/saccharis_plotting" $famname results/cazy_characterization/deeploc/"$famname".extracted.txt
		done
		touch {output.checkpoint}
			#RScript bin/annotate_saccharis_trees.R {params.wd} "results/cazy_characterization/saccharis/"$treename 
		"""

