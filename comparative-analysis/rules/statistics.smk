rule genome_statistics:
	input:
		config["stats_genomes"]
	output:
		checkpoint = "results/checkpoints/genome_statistics.done"
	params:
		prefix=config["prefix"],
		wd = os.getcwd(),
		dir = "results/statistics/"
	conda:
		"../envs/rplotting.yml"
	shell:
		"""
		if [[ ! -d results/statistics ]]
		then
			mkdir results/statistics
		fi
		Rscript bin/stats.R {input} {params.wd}/{params.dir}
		touch {output.checkpoint}
		"""

rule summarize_secreted_and_cazy:
	input: 
		gff_dir = expand("data/{pre}/gff_files", pre=config["prefix"])
	output:
		sec_cazy_summary = expand("results/{pre}/secreted_and_cazy/sec_cazy_sum.tsv", pre=config["prefix"]),
		checkpoint = expand("results/{pre}/checkpoints/summarize_secreted_and_cazy.done", pre=config["prefix"]),
	shell:
		"""
		python bin/summarize_secreted_and_cazy_genes.py -gff {input.gff_dir} > {output.sec_cazy_summary}
		touch {output.checkpoint}
		"""

rule plot_genome_overview:
	input:
		cogs_file = expand("data/{pre}/COGS.all.results.csv", pre=config["prefix"]),
		cazy_file = expand("data/{pre}/CAZyme.summary.results.csv", pre=config["prefix"]),
		secmet_file = expand("data/{pre}/SM.summary.results.csv", pre=config["prefix"]),
		stats_file = expand("data/{pre}/genome.stats.summary.csv", pre=config["prefix"]),
		seccazy_file = rules.summarize_secreted_and_cazy.output.sec_cazy_summary
	output:
		genomes_overview = expand("results/{pre}/genome_overviews/genomes_overview.pdf", pre=config["prefix"]),
		checkpoint = expand("results/{pre}/checkpoints/plot_genome_overview.done", pre=config["prefix"])
	conda:
		"../envs/genome_overview.yml"
	params:
		wd = os.getcwd()
	shell:
		"""
		Rscript bin/plot_overview.R {params.wd} {input.cogs_file} {input.cazy_file} {input.secmet_file} {input.stats_file} {input.seccazy_file} {output.genomes_overview}
		touch {output.checkpoint}
		"""	
if config["phylogeny"]["precalculated"] == "no":
	rule orthology_statistics:
	    input:
	        directory = expand("results/{pre}/orthofinder", pre=config["prefix"]),
	        tree = rules.extract_tree.output.ultra_tree,
	        checkpoint1 = rules.infer_orthology.output.checkpoint,
	        checkpoint2 = rules.extract_tree.output.checkpoint
	    output:
	        directory = directory(expand("results/{pre}/orthofinder_statistics/", pre = config["prefix"])),
	        checkpoint = expand("results/{pre}/checkpoints/orthology_statistics.done", pre=config["prefix"])
	    params:
	        prefix = config["prefix"],
	        orthofiles = config["orthofiles"]
	    conda:
	        "../envs/rorthologystatistics.yml"
	    shell:
	        """
	        for i in {params.orthofiles}
	        do
	            echo $i
	            Rscript bin/ortholog_heatmap.R {output.directory}/{params.prefix}_$i.RData {input.tree} {input.directory}/orthofinder/Results_ortho/Comparative_Genomics_Statistics/$i {output.directory}/{params.prefix} $i
	        done
	        touch {output.checkpoint}
	        """
else:
	rule orthology_statistics:
		output:
			checkpoint = expand("results/{pre}/checkpoints/orthology_statistics.done", pre=config["prefix"])
		shell:
			"""
			touch {output.checkpoint}
			"""

rule gene2gene_distance:
	input:
		anno_folder = config["funannotate_input"]["annotations_tsv_folder"]
	output:
		gene_lengths = "results/statistics/gene_length_species_medians.txt",
		gene2gene_distance = "results/statistics/gene2gene_species_medians.txt"
	singularity: "docker://reslp/biopython_plus:1.77"
	shell:
		"""
		files=$(ls {input.anno_folder}/*)
		bin/gene2gene_distance.py -i $files 
		mv gene_length_species_medians.txt {output.gene_lengths}
		mv gene2gene_species_medians.txt {output.gene2gene_distance}
		"""
rule get_cellulase_orthologs:
	input:
		cellulase_names = config["other_input"]["putative_cellulases"],
		ids = "data/ids.txt"
	output:
		cellulase_orthologs = "results/statistics/cellulase_orthologs/cellulase_orthologs.tsv",
		cellulase_orthologs_renamed = "results/statistics/cellulase_orthologs/cellulase_orthologs_renamed.tsv"
	shell:
		"""
		searchstring=$(cat data/settings_and_parameters/cloned_putative_cellulases.txt | tr "\\n" "-" | sed 's/.$//' | sed -e 's#-#i\\\|#')
		echo $searchstring
		filename=$(grep -l "$searchstring" results/orthology/orthofinder/Results_ortho/Orthogroup_Sequences/*.fa)
	 	echo $filename
		cat $filename | grep ">" | cut -d"_" -f 1 | cut -d ">" -f 2 | sort | uniq -c | sed 's/pred//' | awk '{{print $2 "\t" $1}}' > {output.cellulase_orthologs}
		awk 'NR==FNR{{a[$1]=$2; next}}{{$1=a[$1]; print}}' {input.ids} {output.cellulase_orthologs} > {output.cellulase_orthologs_renamed}
		cp $filename results/statistics/cellulase_orthologs/
		""" 

rule genome_overview:
	input:
		cogs_file = config["funannotate_input"]["cogs"], 
		cazy_file = config["funannotate_input"]["cazy_summary"],
		secmet_file = config["funannotate_input"]["secmet"],
		stats_file = config["funannotate_input"]["genome_stats"],
		gene_length = rules.gene2gene_distance.output.gene_lengths,
		color_file = config["other_input"]["color_information"],
		gene2gene_distance = rules.gene2gene_distance.output.gene2gene_distance,
		tree_file = rules.extract_tree.output.ultra_tree,
		taxonomy_file = config["other_input"]["taxonomy_information"]
	output:
		checkpoint = "results/checkpoints/genome_overview.done",
		p1 = "results/statistics/genomes_overview.pdf",
		p2 = "results/statistics/Rplots.pdf",
		p3 = "results/statistics/cazymes_overview_plus_transporters.pdf",
		stats_file = "results/statistics/median_values.txt"
	singularity: "docker://reslp/rphylogenetics:4.0.3_test"
	shadow: "shallow"
	shell:
		"""
		WD=$(pwd)
		Rscript bin/plot_overview_new_compressed.R $WD {input.cogs_file} {input.cazy_file} {input.secmet_file} {input.stats_file} {input.gene2gene_distance} {input.gene_length} {input.tree_file} {input.color_file} {input.taxonomy_file}
		cp genomes_overview.pdf {output.p1}
		cp Rplots.pdf {output.p2}
		cp cazymes_overview_plus_transporters.pdf {output.p3}
		cp median_values_taxonomy.csv {output.stats_file}	
		touch {output.checkpoint}
		"""
