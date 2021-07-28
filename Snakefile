singularity: "docker://continuumio/miniconda3:4.7.10"
#singularity: "docker://conda/miniconda3-centos7"
configfile: "data/config.yaml"

#print(config)

#rule all:
#	input:
#		#expand("results/{pre}/checkpoints/get_upstream_inputfiles.done", pre=config["prefix"]),
#		expand("results/{pre}/checkpoints/statistics.done", pre=config["prefix"]),
#		#expand("results/{pre}/checkpoints/infer_orthology.done", pre=config["prefix"]),
#		#expand("results/{pre}/checkpoints/rename_ortholog_sequences.done", pre=config["prefix"]),
#		#expand("results/{pre}/checkpoints/align_aa.done", pre=config["prefix"]),
#		#expand("results/{pre}/checkpoints/trim.done", pre=config["prefix"]),
#		expand("results/{pre}/checkpoints/iqtree_concat.done", pre=config["prefix"]),
#		expand("results/{pre}/checkpoints/iqtree_gene_trees.done", pre=config["prefix"]),
#		expand("results/{pre}/checkpoints/iqtree_gene_concordance.done", pre=config["prefix"]),
#		expand("results/{pre}/checkpoints/reroot_tree.done", pre=config["prefix"]),
#		expand("results/{pre}/checkpoints/create_r8s_controlfile.done", pre=config["prefix"]),
#		expand("results/{pre}/checkpoints/run_r8s.done", pre=config["prefix"]),
#		expand("results/{pre}/checkpoints/extract_tree.done", pre=config["prefix"]),
#		expand("results/{pre}/checkpoints/cazy_anc_summary.done", pre=config["prefix"]),
#		#expand("results/{pre}/checkpoints/ancestral_states_cazy.done", pre=config["prefix"]),
#		#expand("results/{pre}/checkpoints/similarity_clustering.done", pre=config["prefix"]),
#		expand("results/{pre}/checkpoints/orthology_statistics.done", pre=config["prefix"]),
#		expand("results/{pre}/checkpoints/plot_phylogeny.done", pre=config["prefix"]),
#		expand("results/{pre}/checkpoints/create_gene_family_table.done", pre=config["prefix"]),
#		# next lines are CAFE related parts
#		#expand("results/{pre}/checkpoints/filter_cafe_family_table.done", pre=config["prefix"]),
#		#expand("results/{pre}/checkpoints/create_cafe_style_tree.done", pre=config["prefix"]),
#		#expand("results/{pre}/checkpoints/create_cafe_commands.done", pre=config["prefix"]),
#		#expand("results/{pre}/checkpoints/run_cafe.done", pre=config["prefix"]),
#		#expand("results/{pre}/checkpoints/parse_cafe_output.done", pre=config["prefix"]),
#		#expand("results/{pre}/checkpoints/visualize_cafe_output.done", pre=config["prefix"]),
#		#expand("results/{pre}/checkpoints/extract_functional_from_cafe.done", pre=config["prefix"]),
#		#expand("results/{pre}/checkpoints/visualize_function_annotations_cafe.done", pre=config["prefix"]),
#		expand("results/{pre}/checkpoints/astral_species_tree.done", pre=config["prefix"]),
#		#expand("results/{pre}/checkpoints/extract_gh5_genes.done", pre=config["prefix"]),
#		#expand("results/{pre}/checkpoints/combine_gh5_genes_and_reference.done", pre=config["prefix"]),
#		#expand("results/{pre}/checkpoints/align_gh5.done", pre=config["prefix"]),
#		#expand("results/{pre}/checkpoints/trim_gh5.done", pre=config["prefix"]),
#		#expand("results/{pre}/checkpoints/iqtree_gh5.done", pre=config["prefix"]),
#		expand("results/{pre}/checkpoints/ancestral_states_cazy_all.done", pre=config["prefix"]),
#		expand("results/{pre}/checkpoints/plot_genome_overview.done", pre=config["prefix"]),
#		expand("results/{pre}/checkpoints/summarize_secreted_and_cazy.done", pre=config["prefix"]),
#		#expand("results/{pre}/checkpoints/plot_ancestral_states_cazy_all.done", pre=config["prefix"]),
#		expand("results/{pre}/checkpoints/character_correlation.done", pre=config["prefix"]),
#		expand("results/{pre}/checkpoints/extract_cazy_proteins.done", pre=config["prefix"]),
#		expand("results/{pre}/checkpoints/saccharis.done", pre=config["prefix"]),
#		expand("results/{pre}/checkpoints/plot_saccharis_trees.done", pre=config["prefix"]),
#		expand("results/{pre}/checkpoints/phylosig.done", pre=config["prefix"]),
#		expand("results/{pre}/checkpoints/pca.done", pre = config["prefix"]),
#		#expand("results/{pre}/checkpoints/create_codon_alignments.done", pre=config["prefix"]),
#		#expand("results/{pre}/checkpoints/run_codeml.done", pre=config["prefix"])

rule all:
	input:
		"results/checkpoints/ancestral_states.done",
		"results/checkpoints/cazy_characterization.done",
		expand("results/{pre}/checkpoints/characterize_transporters.done", pre=config["prefix"]),
		"results/checkpoints/statistics.done",
		"results/checkpoints/orthology.done"

rule phylogeny:
	input:
		#expand("results/{pre}/checkpoints/iqtree_concat.done", pre=config["prefix"]),
		#expand("results/{pre}/checkpoints/iqtree_gene_trees.done", pre=config["prefix"]),
		#expand("results/{pre}/checkpoints/iqtree_gene_concordance.done", pre=config["prefix"]),
		#expand("results/{pre}/checkpoints/reroot_tree.done", pre=config["prefix"]),
		#expand("results/{pre}/checkpoints/create_r8s_controlfile.done", pre=config["prefix"]),
		#expand("results/{pre}/checkpoints/run_r8s.done", pre=config["prefix"]),
		#expand("results/{pre}/checkpoints/extract_tree.done", pre=config["prefix"]),
		"results/checkpoints/plot_phylogeny.done"
	output:
		"results/checkpoints/all_phylogeny.done"
	shell:
		"""
		touch {output}
		"""

rule ancestral_states:
	input:
		rules.phylogeny.output,
		#"results/checkpoints/plot_ancestral_states_cazy_all.done",
		"results/checkpoints/plot_ancestral_states.done"
	output:
		"results/checkpoints/ancestral_states.done"
	shell:
		"""
		touch {output}
		"""

rule geneset_similarity:
	input:
		"results/geneset_similarity/pca.done",
		"results/geneset_similarity/overview_cazyme_counts.csv"

terms = config["cazy_characterization"]["terms"].split(" ")
rule cazy_characterization:
	input:
		expand("results/cazy_characterization/saccharis/saccharis_{term}.done", term=terms),
		expand("results/cazy_characterization/deeploc/checkpoints/deeploc_{term}.done", term=terms),
		expand("results/cazy_characterization/saccharis_plotting/checkpoints/plot_saccharis_tree_{term}.done", term=terms)
	output:
		"results/checkpoints/cazy_characterization.done"
	shell:
		"""
		touch {output}
		"""
rule statistics:
	input:
		"results/checkpoints/genome_overview.done",
		"results/checkpoints/genome_statistics.done",
		"results/statistics/cellulase_orthologs/cellulase_orthologs.tsv"
	output:
		"results/checkpoints/statistics.done"
	shell:
		"""
		touch {output}
		"""
rule orthology:
	input: "results/orthology/checkpoints/infer_orthology.done"
	output: "results/checkpoints/orthology.done"
	shell:
		"""
		touch {output}
		"""
			
#rule all_transporter_tree:
#	input:
#                #expand("results/{pre}/checkpoints/extract_transporter_genes.done", pre=config["prefix"]),
#                #expand("results/{pre}/checkpoints/combine_transporter_genes_and_reference.done", pre=config["prefix"]),
#                #expand("results/{pre}/checkpoints/align_transporter.done", pre=config["prefix"]),
#                #expand("results/{pre}/checkpoints/trim_transporter.done", pre=config["prefix"]),
#                expand("results/{pre}/checkpoints/iqtree_transporter.done", pre=config["prefix"])
#	output:
#		expand("results/{pre}/checkpoints/characterize_transporters.done", pre=config["prefix"])
#	shell:
#		"""
#		touch {output}
#		"""

rule characterize_transporters:
	input:
		"results/transporter_characterization/orthofinder.done"

rule gene_family_evolution:
	input:
		"results/gene_family_evolution/cafe/raw_cafe_input_file.tab",
		"results/gene_family_evolution/cafe/cafe_tree.tre",
		"results/gene_family_evolution/checkpoints/run_cafe.done",
		"results/gene_family_evolution/cafe_significant_families_per_model.txt",
		"results/gene_family_evolution/cafe_summary_tree.pdf"
	output:
		"results/gene_family_evolution/gene_family_evolution.done"
	shell:
		"""
		touch {output}
		"""


def remove_donefile(files):
	new_files = []
	for file in files:
		if ".done" in file:
			continue
		else:
			new_files.append(file)
	return new_files

include: "rules/functions.smk"
include: "rules/orthofinder.smk"
include: "rules/phylogeny.smk"
include: "rules/statistics.smk"
include: "rules/peroxidase_characterization.smk"
include: "rules/ancestral_state_reconstruction.smk"
include: "rules/geneset_similarity.smk"
include: "rules/saccharis.smk"
include: "rules/character_correlation.smk"
include: "rules/cafe.smk"
include: "rules/transporter_phylogeny.smk"
