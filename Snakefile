singularity: "docker://continuumio/miniconda3:4.7.10"
#singularity: "docker://conda/miniconda3-centos7"
configfile: "data/config.yaml"

#print(config)

rule all:
	input:
		expand("results/{pre}/checkpoints/get_upstream_inputfiles.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/statistics.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/infer_orthology.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/rename_ortholog_sequences.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/align_aa.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/trim.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/iqtree_concat.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/iqtree_gene_trees.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/iqtree_gene_concordance.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/reroot_tree.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/create_r8s_controlfile.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/run_r8s.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/extract_tree.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/cazy_anc_summary.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/ancestral_states_cazy.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/similarity_clustering.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/orthology_statistics.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/plot_phylogeny.done", pre=config["prefix"]),
		#expand("results/{pre}/checkpoints/create_gene_family_table.done", pre=config["prefix"]),
		#expand("results/{pre}/checkpoints/filter_cafe_family_table.done", pre=config["prefix"]),
		#expand("results/{pre}/checkpoints/create_cafe_style_tree.done", pre=config["prefix"]),
		#expand("results/{pre}/checkpoints/create_cafe_commands.done", pre=config["prefix"]),
		#expand("results/{pre}/checkpoints/run_cafe.done", pre=config["prefix"]),
		#expand("results/{pre}/checkpoints/parse_cafe_output.done", pre=config["prefix"]),
		#expand("results/{pre}/checkpoints/visualize_cafe_output.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/astral_species_tree.done", pre=config["prefix"]),
		#expand("results/{pre}/checkpoints/extract_functional_from_cafe.done", pre=config["prefix"]),
		#expand("results/{pre}/checkpoints/visualize_function_annotations_cafe.done", pre=config["prefix"]),
		#expand("results/{pre}/checkpoints/extract_gh5_genes.done", pre=config["prefix"]),
		#expand("results/{pre}/checkpoints/combine_gh5_genes_and_reference.done", pre=config["prefix"]),
		#expand("results/{pre}/checkpoints/align_gh5.done", pre=config["prefix"]),
		#expand("results/{pre}/checkpoints/trim_gh5.done", pre=config["prefix"]),
		#expand("results/{pre}/checkpoints/iqtree_gh5.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/ancestral_states_cazy_all.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/plot_genome_overview.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/summarize_secreted_and_cazy.done", pre=config["prefix"]),
		#expand("results/{pre}/checkpoints/plot_ancestral_states_cazy_all.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/character_correlation.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/extract_cazy_proteins.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/saccharis.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/plot_saccharis_trees.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/phylosig.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/pca.done", pre = config["prefix"]),
		#expand("results/{pre}/checkpoints/create_codon_alignments.done", pre=config["prefix"]),
		#expand("results/{pre}/checkpoints/run_codeml.done", pre=config["prefix"])

def remove_donefile(files):
	new_files = []
	for file in files:
		if ".done" in file:
			continue
		else:
			new_files.append(file)
	return new_files
include: "rules/orthofinder.smk"
include: "rules/phylogeny.smk"
include: "rules/statistics.smk"
include: "rules/ancestral_state_reconstruction.smk"
include: "rules/saccharis.smk"
include: "rules/character_correlation.smk"
