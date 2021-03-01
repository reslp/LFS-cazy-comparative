rule statistics:
	input:
		expand("data/{pre}/stats_genomes.csv", pre=config["prefix"])
	output:
		checkpoint = expand("results/{pre}/checkpoints/statistics.done", pre=config["prefix"]),
		dir = directory(expand("results/{pre}/statistics/", pre=config["prefix"]))
	params:
		prefix=config["prefix"],
		wd = os.getcwd()
	conda:
		"../envs/rplotting.yml"
	shell:
		"""
		if [[ ! -d results/{params.prefix}/statistics ]]
		then
			mkdir results/{params.prefix}/statistics
		fi
		Rscript bin/stats.R {input} {params.wd}/{output.dir}
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
