species_list = [sp.split("/")[-1].split(".gff3")[0] for sp in glob.glob(config["funannotate_input"]["annotations_gff_folder"] + "/*.gff3")]

rule create_diamond_db:
	input:
		protein_file = config["funannotate_input"]["protein_folder"] + "/{sp}.proteins.fa"
	output:
		diamond_db = "results/peroxidase_characterization/diamond_dbs/{sp}.dmnd"
	params:
		species = "{sp}",
		wd = os.getcwd()
	singularity:
		"docker://reslp/diamond:0.9.22"
	shell:
		"""
		diamond makedb --in {input.protein_file} -d {params.wd}/results/peroxidase_characterization/diamond_dbs/{params.species}
		"""

rule peroxidase_diamond:
	input:
		protein_file = config["funannotate_input"]["protein_folder"] + "/{sp}.proteins.fa",
		diamond_db = rules.create_diamond_db.output.diamond_db,
		peroxi_seqs = config["other_input"]["peroxidase_sequences"]
	output:
		diamond_results = "results/peroxidase_characterization/diamond_results/{sp}_diamond_results.txt"
	singularity:
		"docker://reslp/diamond:0.9.22"
	threads: 4
	shell:
		"""
		diamond blastp --threads {threads} -d {input.diamond_db} -q {input.peroxi_seqs} -o {output.diamond_results}
		if [[ ! -f {output.diamond_results} ]]; then
			touch {output.diamond_results}
		fi
		"""
rule filter_peroxidase_sequences:
	input:
		diamond_results = rules.peroxidase_diamond.output.diamond_results,
		protein_file = config["funannotate_input"]["protein_folder"] + "/{sp}.proteins.fa"
	output:
		peroxidases = "results/peroxidase_characterization/peroxidase_sequences/{sp}.fa"

	params:
		species = "{sp}"
	singularity:
		"docker://reslp/biopython_plus:1.77"
	shell:
		"""
		cat {input.diamond_results} | awk '{{print $2}}' | sort | uniq > results/peroxidase_characterization/diamond_results/{params.species}_ids.txt 
		bin/select_transcripts.py -i results/peroxidase_characterization/diamond_results/{params.species}_ids.txt -f {input.protein_file} > {output.peroxidases}
		"""


rule orthofinder_peroxi:
	input: 
		protein_files = expand("results/peroxidase_characterization/peroxidase_sequences/{sp}.fa", sp=species_list),
		peroxi_seqs = config["other_input"]["peroxidase_sequences"]
	output:
		"results/peroxidase_characterization/orthofinder.done"
	singularity:
		"docker://reslp/orthofinder:2.5.2"
	threads: config["threads"]["infer_orthology"]
	shell:
		"""
		#cat {input.peroxi_seqs} | sed 's/|/_/g' | sed 's/[][ ]//g' > results/peroxidase_characterization/sequences/Characterized_peroxidases.fa
		cat {input.peroxi_seqs} | awk -F"|" '{{if ($1 ~ />/ ) {{print $1"_"$2}} else {{print}}}}' > results/peroxidase_characterization/peroxidase_sequences/Characterized_peroxidases.fa
		echo "Ulimit hard: "$(ulimit -Hn)
		echo "Ulimit soft: "$(ulimit -Sn)	
		orthofinder -f results/peroxidase_characterization/peroxidase_sequences -o results/peroxidase_characterization/orthofinder -n ortho -S diamond -t {threads}
		touch {output}
		"""

rule parse_orthofinder_peroxi:
	input:
		rules.orthofinder_peroxi.output,
		peroxi_seqs = config["other_input"]["peroxidase_sequences"]
	output:
		pod_mapping = "results/peroxidase_characterization/POD_name_mapping.txt",
		gene_counts = "results/peroxidase_characterization/Orthogroups.GeneCount.renamed.tsv"
	singularity:
		"docker://reslp/biopython_plus:1.77"
	shell:
		"""
		cat {input.peroxi_seqs} | grep ">" | sed 's/>//' | awk -F"|" '{{printf$1"_"$2"\\t"$5"\\n"}}' > {output.pod_mapping}
		bin/rename_orthofinder_counts.py --orthogroups-file results/peroxidase_characterization/orthofinder/Results_ortho/Orthogroups/Orthogroups.txt --mapping-file {output.pod_mapping} --counts-file results/peroxidase_characterization/orthofinder/Results_ortho/Orthogroups/Orthogroups.GeneCount.tsv > {output.gene_counts}
		"""


rule plot_peroxidases:
	input:
		countsfile = rules.parse_orthofinder_peroxi.output.gene_counts,
		treefile = rules.extract_tree.output.ultra_tree		
	output:
		"results/peroxidase_characterization/pod_heatmap.pdf"
	singularity:
		"docker://reslp/rphylogenetics:4.0.3"
	params:
		wd = os.getcwd()
	shell:
		"""
		Rscript bin/POD_heatmap.R {params.wd}/results/peroxidase_characterization {input.treefile} {input.countsfile}	
		"""
