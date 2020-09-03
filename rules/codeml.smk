#these rules are currently not used
rule create_codon_alignments:
	input:
		nucl_file = expand("data/{pre}/{pre}_transcripts.fas", pre=config["prefix"]),
		aa_alignments = rules.align_aa.output.dir,
		checkpoint = rules.align_aa.output.checkpoint,
		orth_dir = expand("results/{pre}/orthofinder/", pre=config["prefix"]),
		checkpoint2 = rules.infer_orthology.output.checkpoint,
		id = expand("data/{pre}/ids.txt", pre=config["prefix"])
	output:
		dir = directory(expand("results/{pre}/codon_alignments/", pre=config["prefix"])),
		checkpoint = expand("results/{pre}/checkpoints/create_codon_alignments.done", pre=config["prefix"])
	conda:
		"envs/pyutils.yml"
	shell:
		"""
		for file in $(find {input.aa_alignments} -type f -name "*.fa_renamed_aligned");
		do
			echo "Codon alignment for file: "$file
			outname=$(basename $file)
			echo $outname
			python -W ignore bin/create_codon_alignment.py -prot {input.aa_alignments}/"$outname" -nucl {input.nucl_file} -orth {input.orth_dir}/orthofinder/Results_ortho/Orthogroups/Orthogroups.tsv -id {input.id}> {output.dir}/"$outname"_codon
		done
		touch {output.checkpoint}
		"""
rule run_codeml:
# this rule needs xvfb-run which creates a fake X Server to plot the created tree. sudo apt-get install xvfb
	input:
		dir = rules.create_codon_alignments.output.dir,
		checkpoint1 = rules.create_codon_alignments.output.checkpoint,
		tree = rules.extract_tree.output,
		checkpoint2 = rules.extract_tree.output.checkpoint
	output:
		dir = directory(expand("results/{pre}/codeml/", pre=config["prefix"])),
		checkpoint = expand("results/{pre}/checkpoints/run_codeml.done", pre=config["prefix"])
	singularity:
		"docker://reslp/ete3:3.1.1"
	threads: 8
	shell:
		"""
		for file in $(find {input.dir} -type f -size +0c -name "*_codon");
		do
			date
			echo "Running evolution analysis on file: "$file
			outname=$(basename $file | cut -d. -f1)
			#ete3 evol -v 3 -t {input.tree} --alg $file --models M0 M2 bsA bsA1 --mark Xylographa_parallela,,,Xylographa_trunciseda -o {output.dir}"$outname"_results --clear_all --cpu {threads} --tests bsA,bsA1 --noimg --codeml_param cleandata,1 &> {output.dir}"$outname".log
			#ete3 evol -v 3 -t {input.tree} --alg $file --models M0 -o {output.dir}"$outname"_results --clear_all --cpu {threads} --noimg --codeml_param runmode,0 cleandata,1 fix_blength,1 verbose,1 noisy,1 &> {output.dir}"$outname".log
			ete3 evol -v 3 -t {input.tree} --alg $file --models M0 -o {output.dir}"$outname"_results --clear_all --cpu {threads} --codeml_param method,1 cleandata,1 fix_blength,-1 verbose,1 noisy,2 --noimg &> {output.dir}"$outname".log
		done
		"""
