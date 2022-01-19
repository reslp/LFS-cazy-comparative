rule pca:
	input:
		cazy_data = rules.ancestral_states_all_cazy.output.rdata,
		check1 = rules.ancestral_states_all_cazy.output.checkpoint,
		apriori_sets = "data/settings_and_parameters/apriori_cazyme_sets.txt",
		genome_stats = config["other_input"]["taxonomy_information"],
		colors_file = "data/settings_and_parameters/color_information.csv",
		peroxi_gene_counts = rules.parse_orthofinder_peroxi.output.gene_counts
	output:
		checkpoint = "results/geneset_similarity/pca.done",
	params:
		wd = os.getcwd(),
		prefix = config["prefix"],
		dir = directory("results/geneset_similarity/")

	#conda:
	#	"../envs/rreroot.yml"
	singularity:
		"docker://reslp/rphylogenetics:4.0.3"
	shell:
		"""
		Rscript bin/phyl_pca.R {input.cazy_data} {input.colors_file} {input.apriori_sets} {input.genome_stats} {input.peroxi_gene_counts} {params.dir} 
		touch {output.checkpoint}
		"""

rule cazy_statistics:
	input:
		taxonomy = config["other_input"]["taxonomy_information"],
		cazy_info = config["funannotate_input"]["cazy_summary"]
	output:
		"results/geneset_similarity/overview_cazyme_counts.csv"
	params:
		wd = os.getcwd()
	singularity:
		"docker://reslp/rphylogenetics:4.0.3"
	shell:
		"""
		Rscript bin/cazy_statistics.R {params.wd}/results/geneset_similarity/ {params.wd}/{input.taxonomy} {params.wd}/{input.cazy_info}
		"""

