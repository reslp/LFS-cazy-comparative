ruleorder: extract_transporter_genes > orthofinder_transporter

species_list = [sp.split("/")[-1].split(".gff3")[0] for sp in glob.glob(config["funannotate_input"]["annotations_gff_folder"] + "/*.gff3")]

rule extract_transporter_genes:
	input: 
		gff_file = config["funannotate_input"]["annotations_gff_folder"] + "/{sp}.gff3",
		protein_file = config["funannotate_input"]["protein_folder"] + "/{sp}.proteins.fa"
	output:
		transporter_proteins = "results/transporter_characterization/proteins/{sp}.fa"
	singularity:
		"docker://reslp/biopython_plus:1.77"
	shell:
		"""
		bin/select_genes_for_annotation.py -gff {input.gff_file} -fasta {input.protein_file} -item PFAM:PF00083 > {output.transporter_proteins}
		"""

rule orthofinder_transporter:
	input:
		prot_files = expand("results/transporter_characterization/proteins/{sp}.fa", sp=species_list),
		characterized_genes = config["other_input"]["characterized_transporters"]
	output:
		"results/transporter_characterization/orthofinder.done"
	singularity:
		"docker://reslp/orthofinder:2.5.2"
	threads: config["threads"]["infer_orthology"]
	shell:
		"""
		cp {input.characterized_genes} results/transporter_characterization/proteins/Characterized_transporters.fa
		ulimit -n 3000	
		orthofinder -f results/transporter_characterization/proteins -o results/transporter_characterization/orthofinder -n ortho -S diamond -t {threads}
		touch {output}
		"""
		
		
		

#rule extract_transporter_genes:
#	input:
#		gff_files = expand("data/{pre}/gff_files", pre=config["prefix"]),
#		prot_files = expand("data/{pre}/protein_files", pre=config["prefix"]),
#	output:
#		transporter_genes = expand("results/{pre}/PF00083_tree/extracted_transporter_genes_{pre}.fas", pre = config["prefix"]),
#		checkpoint = expand("results/{pre}/checkpoints/extract_transporter_genes.done", pre=config["prefix"])
#	conda:
#		"../envs/pyutils.yml"
#	shell:
#		"""
#		bin/select_genes_for_annotation.py -gff {input.gff_files} -fasta {input.prot_files} -item PFAM:PF00083 > {output.transporter_genes}
#		touch {output.checkpoint}
#		"""
#
#rule combine_transporter_genes_and_reference:
#	input:
#		known_transporter_genes = "data/sugar_transporters_characterized_reference_seqs.fa",
#		extracted_transporter_genes = rules.extract_transporter_genes.output.transporter_genes,
#	output:
#		all_transporter_genes = expand("results/{pre}/PF00083_tree/all_transporter_genes_{pre}.fas", pre = config["prefix"]),
#		checkpoint = expand("results/{pre}/checkpoints/combine_transporter_genes_and_reference.done", pre=config["prefix"])
#	shell:
#		"""
#		cat {input.known_transporter_genes} {input.extracted_transporter_genes} > {output.all_transporter_genes}
#		touch {output.checkpoint}
#		"""
#
#rule align_transporter:
#	input:
#		all_transporter_genes = rules.combine_transporter_genes_and_reference.output.all_transporter_genes
#	log:
#		"log/align_transporter.log"
#	output:
#		checkpoint = expand("results/{pre}/checkpoints/align_transporter.done", pre=config["prefix"]),
#		transporter_alignment =  expand("results/{pre}/PF00083_tree/all_transporter_genes_alignment_{pre}.fas", pre = config["prefix"]),
#	conda:
#		"../envs/phylogenomics.yml"
#	shell:
#		"""
#		mafft --quiet --auto {input.all_transporter_genes}  > {output.transporter_alignment}
#		touch {output.checkpoint}
#		"""
#
#rule trim_transporter:
#	input:
#		rules.align_transporter.output.transporter_alignment
#	output:
#		checkpoint = expand("results/{pre}/checkpoints/trim_transporter.done", pre=config["prefix"]),
#		trimmed_alignment = expand("results/{pre}/PF00083_tree/all_transporter_genes_alignment_trimmed_{pre}.fas", pre = config["prefix"])
#	conda:
#		"../envs/phylogenomics.yml"
#	shell:
#		"""
#		trimal -gappyout -in {input} -out {output.trimmed_alignment}
#		touch {output.checkpoint}
#		"""
#
#rule iqtree_transporter:
#	input:
#		rules.trim_transporter.output.trimmed_alignment
#	output:
#		checkpoint = expand("results/{pre}/checkpoints/iqtree_transporter.done", pre=config["prefix"]),
#	params:
#		wd = os.getcwd(),
#		prefix = config["prefix"]
#	singularity:
#		"docker://reslp/iqtree:2.0rc2"
#	threads:
#		config["threads"]["iqtree_transporter"]
#	shell:
#		"""
#		cd results/{params.prefix}/PF00083_tree
#		iqtree -s {params.wd}/{input} --prefix PF00083_tree -bb 1000 -nt {threads} -m TEST -redo --no-terrace
#		cd {params.wd}
#		touch {output.checkpoint}
#		"""
