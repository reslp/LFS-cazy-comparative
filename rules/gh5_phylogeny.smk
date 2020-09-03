# the following rules are currently not included in the analysis. There is now a better way using Saccharis and Deeploc
rule extract_gh5_genes:
	input:
		gff_files = expand("data/{pre}/gff_files/", pre=config["prefix"]),
		prot_files = expand("data/{pre}/protein_files/", pre=config["prefix"]),
	output:
		gh5_genes = expand("results/{pre}/gh5_tree/extracted_gh5_genes_{pre}.fas", pre = config["prefix"]),
		checkpoint = expand("results/{pre}/checkpoints/extract_gh5_genes.done", pre=config["prefix"])
	conda:
		"envs/pyutils.yml"
	shell:
		"""
		bin/select_genes_for_annotation.py -gff {input.gff_files} -fasta {input.prot_files} -item CAZy:GH5 > {output.gh5_genes}
		touch {output.checkpoint}
		"""

rule combine_gh5_genes_and_reference:
	input:
		known_gh5_genes = "data/cazy_gh5_characterized_reference_seqs.fa",
		extracted_gh5_genes = rules.extract_gh5_genes.output.gh5_genes,
	output:
		all_gh5_genes = expand("results/{pre}/gh5_tree/all_gh5_genes_{pre}.fas", pre = config["prefix"]),
		checkpoint = expand("results/{pre}/checkpoints/combine_gh5_genes_and_reference.done", pre=config["prefix"])
	shell:
		"""
		cat {input.known_gh5_genes} {input.extracted_gh5_genes} > {output.all_gh5_genes}
		touch {output.checkpoint}
		"""

rule align_gh5:
	input:
		all_gh5_genes = rules.combine_gh5_genes_and_reference.output.all_gh5_genes
	log:
		"log/align_gh5.log"
	output:
		checkpoint = expand("results/{pre}/checkpoints/align_gh5.done", pre=config["prefix"]),
		gh5_alignment =  expand("results/{pre}/gh5_tree/all_gh5_genes_alignment_{pre}.fas", pre = config["prefix"]),
	conda:
		"envs/phylogenomics.yml"
	shell:
		"""
		mafft --quiet --auto {input.all_gh5_genes}  > {output.gh5_alignment}
		touch {output.checkpoint}
		"""

rule trim_gh5:
	input:
		rules.align_gh5.output.gh5_alignment
	output:
		checkpoint = expand("results/{pre}/checkpoints/trim_gh5.done", pre=config["prefix"]),
		trimmed_alignment = expand("results/{pre}/gh5_tree/all_gh5_genes_alignment_trimmed_{pre}.fas", pre = config["prefix"])
	conda:
		"envs/phylogenomics.yml"
	shell:
		"""
		trimal -gappyout -in {input} -out {output.trimmed_alignment}
		touch {output.checkpoint}
		"""

rule iqtree_gh5:
	input:
		rules.trim_gh5.output.trimmed_alignment
	output:
		checkpoint = expand("results/{pre}/checkpoints/iqtree_gh5.done", pre=config["prefix"]),
	params:
		wd = os.getcwd(),
		prefix = config["prefix"]
	singularity:
		"docker://reslp/iqtree:2.0rc2"
	threads:
		config["iqtree"]["threads"]
	shell:
		"""
		cd results/{params.prefix}/gh5_tree
		iqtree -s {params.wd}/{input} --prefix gh5_tree -bb 1000 -nt {threads} -m TEST -redo --no-terrace
		cd {params.wd}
		touch {output.checkpoint}
		"""
