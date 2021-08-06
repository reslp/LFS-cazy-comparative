singularity: "docker://continuumio/miniconda3:4.7.10"
configfile: "data/config.yaml"


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

rule reconstruct_ancestral_states:
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
rule characterize_cazymes:
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

rule characterize_transporters:
	input:
		"results/transporter_characterization/orthofinder.done"

rule infer_gene_family_evolution:
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

rule characterize_peroxidases:
	input:
		"results/peroxidase_characterization/pod_heatmap.pdf"


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
include: "rules/transporter_characterization.smk"
