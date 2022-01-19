# attention: this uses a trick. Singularity does not resule symlinks to the folder that are above the folder from where singularity is run.
# # However snakemake needs the files in the correct places to find them as inputs.
# # I had to create the symblink normally: data/upstream => linking to the upstream funannotate analysis
# # next I need to bind the directory in singularity with -B (in submit.sh) to overwrite the symlink to be able to move into subfolders of this folder inside the running container
#
# #idsfile = "data/upstream/results/downstream/ids.txt",
# #combined_gff = "data/upstream/results/downstream/combined.gff",
# #combined_proteins = "data/upstream/results/downstream/combined_proteins.fa",
# #combined_transcripts = "data/upstream/results/downstream/combined_proteins.fa",
# #cazyme_results = "data/upstream/results/downstream/CAZyme.all.results.csv", 
# #cazyme_summary_results = "data/upstream/results/downstream/CAZyme.summary.results.csv",
# #interproscan_results = "data/upstream/results/downstream/interproscan.results.csv",
# #pfam_results = "data/upstream/results/downstream/pfam.results.csv",
# #protein_files = "data/upstream/results/downstream/protein_files"
# #cp {input.idsfile} {output.idsfile}
# #cp {input.combined_gff} {output.combined_gff}
# #cp {input.combined_proteins} {output.combined_proteins}
# #cp {input.combined_transcripts} {output.combined_transcripts}
# #cp {input.cazyme_results} {output.cazyme_results}
# #cp {input.cazyme_summary_results} {output.cazyme_summary_results}
# # cp {input.interproscan_results} {output.interproscan_results}
# # cp {input.pfam_results} {output.pfam_results}
# # mkdir -p {output.protein_files}
# # cp -r {input.protein_files}/* {output.protein_files}
# #touch {output.checkpoint}

# currently this rules does not check for the gff and protein file directories!
rule get_upstream_inputfiles:
	input:
		idsfile = expand("data/{pre}/ids.txt", pre=config["prefix"]),
		combined_gff = expand("data/{pre}/{pre}_combined.gff", pre=config["prefix"]),
		combined_proteins = expand("data/{pre}/{pre}_proteins.fas", pre=config["prefix"]),
		combined_transcripts = expand("data/{pre}/{pre}_transcripts.fas", pre=config["prefix"]),
		cazyme_results = expand("data/{pre}/CAZyme.all.results.csv", pre=config["prefix"]),
		cazyme_summary_results = expand("data/{pre}/CAZyme.summary.results.csv", pre=config["prefix"]),
		interproscan_results = expand("data/{pre}/interproscan.results.csv", pre=config["prefix"]),
		pfam_results = expand("data/{pre}/pfam.results.csv", pre=config["prefix"])
	output:	
		checkpoint = expand("results/{pre}/checkpoints/get_upstream_inputfiles.done", pre=config["prefix"]), 
	shell:
		"""
		touch {output.checkpoint}	
		"""

rule infer_orthology:
	input:
		prot_files = get_protein_files
	output:
		checkpoint = "results/orthology/checkpoints/infer_orthology.done"
	params:
		proteins = config["funannotate_input"]["protein_folder"],
		dir = "results/orthology"
	singularity: "docker://reslp/orthofinder:2.5.2"
	threads: config["threads"]["infer_orthology"]
	shell:
		"""
		orthofinder -f {params.proteins} -o {params.dir}/orthofinder -n ortho -S diamond -t {threads}
		touch {output.checkpoint}
		"""

# this rules renames the orthology results from orthofinder so that they can be used in downstream analysis
# the ids file for renaming is specified here, which is probably not optimal and should be moved to the config file
rule rename_ortholog_sequences:
    input:
        checkpoint = rules.infer_orthology.output.checkpoint,
        dir = expand("results/{pre}/orthofinder", pre=config["prefix"]),
        ids_file = expand("data/{pre}/ids.txt", pre=config["prefix"])
    output:
        checkpoint = expand("results/{pre}/checkpoints/rename_ortholog_sequences.done", pre=config["prefix"]),
        dir = directory(expand("results/{pre}/renamed_sc_seq_files", pre=config["prefix"]))
    conda:
        "../envs/pyutils.yml"
    params:
        wd = os.getcwd()
    shell:
        """
        mkdir {output.dir}/aa
        mkdir {output.dir}/nu
        echo "Renaming aa sequences..."
        for file in $(find {input.dir}/orthofinder/Results_ortho/Single_Copy_Orthologue_Sequences/ -type f -name "*.fa");
        do
	    echo $file
            outname=$(basename $file)
            python {params.wd}/bin/rename_sequences.py {input.ids_file} $file > {output.dir}/aa/$outname"_renamed"
        done
        echo "Renaming nu sequences..."
        for file in $(find {input.dir}/orthofinder/Results_ortho/Single_Copy_Orthologue_Sequences/ -type f -name "*.fa");
        do
        	outname=$(basename $file)
        	python {params.wd}/bin/rename_sequences.py {input.ids_file} $file > {output.dir}/nu/$outname"_renamed"
        done
        touch {output.checkpoint}
        """
