rule genome_statistics:
	input:
		expand("data/{pre}/stats_genomes.csv", pre=config["prefix"])
	output:
		checkpoint = expand("results/{pre}/checkpoints/genome_statistics.done", pre=config["prefix"])
	params:
		prefix=config["prefix"],
		wd = os.getcwd(),
		dir = expand("results/{pre}/statistics/", pre=config["prefix"])
	conda:
		"../envs/rplotting.yml"
	shell:
		"""
		if [[ ! -d results/{params.prefix}/statistics ]]
		then
			mkdir results/{params.prefix}/statistics
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
		anno_folder = directory(expand("data/{pre}/annotation_tsv", pre=config["prefix"]))
	output:
		gene_lengths = expand("results/{pre}/statistics/gene_length_species_medians.txt",pre=config["prefix"]),
		gene2gene_distance = expand("results/{pre}/statistics/gene2gene_species_medians.txt", pre=config["prefix"])
	singularity: "docker://reslp/biopython_plus:1.77"
	shell:
		"""
		files=$(ls {input.anno_folder}/*)
		bin/gene2gene_distance.py -i $files 
		mv gene_length_species_medians.txt {output.gene_lengths}
		mv gene2gene_species_medians.txt {output.gene2gene_distance}
		"""

rule genome_overview:
	input:
		cogs_file = expand("data/{pre}/COGS.all.results.csv", pre=config["prefix"]),
		cazy_file = expand("data/{pre}/CAZyme.summary.results.csv", pre=config["prefix"]),
		secmet_file = expand("data/{pre}/SM.summary.results.csv", pre=config["prefix"]),
		stats_file = expand("data/{pre}/genome.stats.summary.csv", pre=config["prefix"]),
		gene_length = rules.gene2gene_distance.output.gene_lengths,
		lifestyle_file = expand("data/{pre}/lifestyle", pre=config["prefix"]),
		gene2gene_distance = rules.gene2gene_distance.output.gene2gene_distance,
		tree_file = rules.extract_tree.output.ultra_tree
	output:
		checkpoint = expand("results/{pre}/checkpoints/genome_overview.done", pre=config["prefix"]),
		p1 = expand("results/{pre}/statistics/genomes.pdf", pre=config["prefix"]),
		p2 = expand("results/{pre}/statistics/Rplots.pdf", pre=config["prefix"]),
		p3 = expand("results/{pre}/statistics/cazymes_overview_plus_transporters.pdf", pre=config["prefix"])
	singularity: "docker://reslp/rphylogenetics:4.0.3"
	shadow: "shallow"
	shell:
		"""
		WD=$(pwd)
		Rscript bin/plot_overview_new_compressed.R $WD {input.cogs_file} {input.cazy_file} {input.secmet_file} {input.stats_file} {input.gene2gene_distance} {input.gene_length} {input.tree_file} {input.lifestyle_file}
		cp genomes.pdf {output.p1}
		cp Rplots.pdf {output.p2}
		cp cazymes_overview_plus_transporters.pdf {output.p3}
		touch {output.checkpoint}
		"""
