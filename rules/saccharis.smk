rule prepare_scrape_cazy:
	input:
		cazy = expand("data/{pre}/CAZyme.all.results.csv", pre=config["prefix"]),
	output:
		checkpoint = expand("results/{pre}/checkpoints/prepare_scrape_cazy.done", pre=config["prefix"]),	
		all_cazy_info = expand("results/{pre}/cazy_information/all_cazy_data.csv", pre=config["prefix"])
	singularity:
		"docker://reslp/scrape_cazy:1"
	params:
		prefix = config["prefix"],
		wd = os.getcwd(),
		search_terms = config["scrape_cazy"]["terms"]
	shell:
		"""
		if [[ ! -d results/{params.prefix}/cazy_information ]]
		then
			mkdir results/{params.prefix}/cazy_information
		fi
		cd results/{params.prefix}/cazy_information
		rm -f all_interesting_cazy_families.txt
		scrape_cazy.py -f $(cat {params.wd}/{input.cazy} | awk -F "," 'NR > 1 {{print $1;}}' | awk -F "_" '{{print $1}}' | uniq | tr "\\n" "," | sed '$s/,$//')
		tail -n +2 -q *_characterized.txt > all_cazy_data.csv
		# strangely this only works if the terms specified in the config file exist. It works outside of snakemake. Maybe it some problem with the different shell 
		# invoked by snakemake...
		names="{params.search_terms}"
		for name in $names; do
			touch "$name"_families.txt
			cat all_cazy_data.csv | grep "$name" | awk '{{print $1}}' | uniq > "$name"_families.txt
		done	
		cat *_families.txt > all_interesting_cazy_families.txt
		cd {params.wd}
		touch {output.checkpoint}
		"""

rule extract_cazy_proteins:
	input:
		combined_proteins = expand("data/{pre}/protein_files", pre=config["prefix"]),
		combined_gff = expand("data/{pre}/gff_files", pre=config["prefix"]),
	output:
		cazy_proteins = expand("results/{pre}/saccharis/all_cazy_proteins.fas", pre=config["prefix"]),
		checkpoint = expand("results/{pre}/checkpoints/extract_cazy_proteins.done", pre=config["prefix"])
	singularity:
		"docker://reslp/biopython_plus:1.77"
	shell:
		"""
		bin/select_cazy_from_gff.py -gff {input.combined_gff} -fasta {input.combined_proteins} -item CAZy >> {output.cazy_proteins}
		touch {output.checkpoint}
		"""

# all rules using the saccharis container need the following bindpoints: "-B $(pwd):/data -B $(pwd)/bin:/usr/local/external"
rule prepare_saccharis:
	input:
		combined_proteins = rules.extract_cazy_proteins.output.cazy_proteins,
		cazy = expand("data/{pre}/CAZyme.all.results.csv", pre=config["prefix"])	
	output:
		renamed_sequences = expand("results/{pre}/saccharis/renamed_proteins.fas", pre=config["prefix"]),
		checkpoint = expand("results/{pre}/checkpoints/prepare_saccharis.done", pre=config["prefix"])
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
		cazy = rules.prepare_scrape_cazy.output.checkpoint
	output:
		checkpoint = expand("results/{pre}/checkpoints/saccharis.done", pre=config["prefix"])
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
		rm -f {params.wd}/results/{params.prefix}/saccharis/interesting_cazy_families.txt
		for cazy in $(cat {params.wd}/results/{params.prefix}/cazy_information/all_interesting_cazy_families.txt | sort | uniq | tr '\\n' ' '); 
		do 
			if ! grep -Fxq "Cazy Pipeline Finished" {params.wd}/results/{params.prefix}/saccharis/$cazy.log;
			then 
				echo $cazy >> {params.wd}/results/{params.prefix}/saccharis/interesting_cazy_families.txt
				rm -rf {params.wd}/results/{params.prefix}/saccharis/$cazy
			else
				continue
			fi
		done
		echo $(cat /data/results/{params.prefix}/saccharis/interesting_cazy_families.txt | tr '\\n' ' ')
		parallel --joblog {params.wd}/results/{params.prefix}/saccharis/parallel_logfile --results $TMPDIR  -j {params.parallel_jobs} Saccharis.pl -d /data/results/{params.prefix}/saccharis -g characterized -s /data/{input.seqs} -t {params.saccharis_threads} -f {{}} "&>" {params.wd}/results/{params.prefix}/saccharis/{{}}.log ::: $(cat /data/results/{params.prefix}/saccharis/interesting_cazy_families.txt | tr '\\n' ' ')
		touch {output.checkpoint}
		"""

rule get_saccharis_mapping_data:
	input:
		seqs = rules.prepare_saccharis.output.renamed_sequences,
		genome_info = expand("data/{pre}/stats_genomes.csv", pre=config["prefix"]),
		discrete_data = expand("data/{pre}/character_information.csv",pre = config["prefix"])
	output:
		saccharis_mapping_data = expand("results/{pre}/saccharis_plotting/saccharis_mapping_information.tsv",pre = config["prefix"]),
		checkpoint = expand("results/{pre}/checkpoints/get_saccharis_mapping_data.done", pre=config["prefix"])
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
		checkpoint = expand("results/{pre}/checkpoints/deeploc.done", pre=config["prefix"]),
		outdir = directory(expand("results/{pre}/deeploc/", pre=config["prefix"]))
	params:
		prefix = config["prefix"],
		wd = os.getcwd()
	singularity:
		"docker://reslp/deeploc:1.0"
	shell:
		"""
		cd /data/{output.outdir}
		for seqfile in $(find /data/results/{params.prefix}/saccharis/*/characterized/dbcan/*extracted.fasta -not -empty); do
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
		checkpoint = rules.prepare_scrape_cazy.output.checkpoint,
		checkpoint2 = rules.deeploc.output.checkpoint
	output:
		checkpoint = expand("results/{pre}/checkpoints/plot_saccharis_trees.done", pre=config["prefix"])
	params:
		wd = os.getcwd(),
		prefix = config["prefix"]
	conda:
		"../envs/rreroot.yml"
	shell:
		"""
		for treename in $(ls results/{params.prefix}/saccharis/*/characterized/*.tree); do	
			echo $treename
			#IFS="_"
			#read -ar fields <<< $(basename $treename)
			famname=$(basename $treename _characterized.tree)
			#famname="${{fields[0]}}"
			echo $famname 
		        Rscript bin/annotate_saccharis_trees.R {params.wd} $treename results/{params.prefix}/cazy_information/"$famname"_characterized.txt {input.saccharis_mapping_data} "results/{params.prefix}/saccharis_plotting" $famname results/{params.prefix}/deeploc/"$famname".extracted.txt
		done
		touch {output.checkpoint}
			#RScript bin/annotate_saccharis_trees.R {params.wd} "results/{params.prefix}/saccharis/"$treename 
		"""

