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
		#expand("results/{pre}/checkpoints/similarity_clustering.done", pre=config["prefix"]),
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
		expand("results/{pre}/checkpoints/extract_gh5_genes.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/combine_gh5_genes_and_reference.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/align_gh5.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/trim_gh5.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/iqtree_gh5.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/ancestral_states_cazy_all.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/plot_genome_overview.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/summarize_secreted_and_cazy.done", pre=config["prefix"]),
		expand("results/{pre}/checkpoints/plot_ancestral_states_cazy_all.done", pre=config["prefix"])
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
# attention: this uses a trick. Singularity does not resule symlinks to the folder that are above the folder from where singularity is run.
# However snakemake needs the files in the correct places to find them as inputs.
# I had to create the symblink normally: data/upstream => linking to the upstream funannotate analysis
# next I need to bind the directory in singularity with -B (in submit.sh) to overwrite the symlink to be able to move into subfolders of this folder inside the running container

#idsfile = "data/upstream/results/downstream/ids.txt",
#combined_gff = "data/upstream/results/downstream/combined.gff",
#combined_proteins = "data/upstream/results/downstream/combined_proteins.fa",
#combined_transcripts = "data/upstream/results/downstream/combined_proteins.fa",
#cazyme_results = "data/upstream/results/downstream/CAZyme.all.results.csv", 
#cazyme_summary_results = "data/upstream/results/downstream/CAZyme.summary.results.csv",
#interproscan_results = "data/upstream/results/downstream/interproscan.results.csv",
#pfam_results = "data/upstream/results/downstream/pfam.results.csv",
#protein_files = "data/upstream/results/downstream/protein_files"
#cp {input.idsfile} {output.idsfile}
#cp {input.combined_gff} {output.combined_gff}
#cp {input.combined_proteins} {output.combined_proteins}
#cp {input.combined_transcripts} {output.combined_transcripts}
#cp {input.cazyme_results} {output.cazyme_results}
#cp {input.cazyme_summary_results} {output.cazyme_summary_results}
# cp {input.interproscan_results} {output.interproscan_results}
# cp {input.pfam_results} {output.pfam_results}
# mkdir -p {output.protein_files}
# cp -r {input.protein_files}/* {output.protein_files}
#touch {output.checkpoint}

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
		"envs/rplotting.yml"
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
		gff_dir = expand("data/{pre}/gff_files/", pre=config["prefix"])
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
		"envs/genome_overview.yml"
	params:
		wd = os.getcwd()
	shell:
		"""
		Rscript bin/plot_overview.R {params.wd} {input.cogs_file} {input.cazy_file} {input.secmet_file} {input.stats_file} {input.seccazy_file} {output.genomes_overview}
		touch {output.checkpoint}
		"""	

rule infer_orthology:
	input:
		prot_files = expand("data/{pre}/protein_files", pre=config["prefix"]),
		check = rules.get_upstream_inputfiles.output.checkpoint 
	output:
		dir = directory(expand("results/{pre}/orthofinder/", pre=config["prefix"])),
		checkpoint = expand("results/{pre}/checkpoints/infer_orthology.done", pre=config["prefix"])
	conda: "envs/orthofinder.yml"
	threads: config["orthofinder"]["threads"]
	shell:
		"""
		orthofinder -f {input.prot_files} -o {output.dir}/orthofinder -n ortho -S diamond -t {threads}
		touch {output.checkpoint}
		"""


# this rules renames the orthology results from orthofinder so that they can be used in downstream analysis
# the ids file for renaming is specified here, which is probably not optimal and should be moved to the config file
rule rename_ortholog_sequences:
    input:
        checkpoint = rules.infer_orthology.output.checkpoint,
        dir = rules.infer_orthology.output.dir,
        ids_file = expand("data/{pre}/ids.txt", pre=config["prefix"])
    output:
        checkpoint = expand("results/{pre}/checkpoints/rename_ortholog_sequences.done", pre=config["prefix"]),
        dir = directory(expand("results/{pre}/renamed_sc_seq_files/", pre=config["prefix"]))
    conda:
        "envs/pyutils.yml"
    params:
        wd = os.getcwd()
    shell:
        """
        mkdir {output.dir}/aa
        mkdir {output.dir}/nu
        echo "Renaming aa sequences..."
        for file in $(find {input.dir}orthofinder/Results_ortho/Single_Copy_Orthologue_Sequences/ -type f -name "*.fa");
        do
            outname=$(basename $file)
            python {params.wd}/bin/rename_sequences.py {input.ids_file} $file > {output.dir}/aa/$outname"_renamed"
        done
        echo "Renaming nu sequences..."
        for file in $(find {input.dir}orthofinder/Results_ortho/Single_Copy_Orthologue_Sequences/ -type f -name "*.fa");
        do
        	outname=$(basename $file)
        	python {params.wd}/bin/rename_sequences.py {input.ids_file} $file > {output.dir}/nu/$outname"_renamed"
        done
        touch {output.checkpoint}
        """

rule align_aa:
    input:
        dir = rules.rename_ortholog_sequences.output.dir,
        checkpoint = rules.rename_ortholog_sequences.output.checkpoint
    log:
        "log/align_aa.log"
    output:
        checkpoint = expand("results/{pre}/checkpoints/align_aa.done", pre=config["prefix"]),
        dir = directory(expand("results/{pre}/aa_alignments/", pre=config["prefix"]))
    conda:
        "envs/phylogenomics.yml"
    shell:
        """
        for file in $(find {input.dir}/aa -type f -name "*_renamed");
        do
            outname=$(basename $file)
            echo "Aligning sequences: "$outname
            mafft --quiet --auto {input.dir}/aa/$outname  > {output.dir}/$outname"_aligned"
        done
        touch {output.checkpoint}
        """

rule trim:
    input:
        dir = rules.align_aa.output.dir,
        checkpoint = rules.align_aa.output.checkpoint
    output:
        checkpoint = expand("results/{pre}/checkpoints/trim.done", pre=config["prefix"]),
        dir = directory(expand("results/{pre}/aa_trimmed/", pre=config["prefix"]))
    conda:
        "envs/phylogenomics.yml"
    shell:
        """
        for file in $(find {input.dir} -type f -name "*_aligned");
        do
            echo "Trimming file: "$file
            outname=$(basename $file)
            trimal -gappyout -in $file -out {output.dir}$outname"_trimmed"
        done
        touch {output.checkpoint}
        """

rule iqtree_concat:
    input:
        trims = rules.trim.output.dir,
        checkpoint = rules.trim.output.checkpoint,
        ids_file = expand("data/{pre}/ids.txt", pre=config["prefix"])
    output:
        checkpoint = expand("results/{pre}/checkpoints/iqtree_concat.done", pre=config["prefix"]),
        tree = expand("results/{pre}/phylogeny/concatenated/concat.treefile", pre=config["prefix"])
    params:
        wd = os.getcwd(),
        prefix = config["prefix"]
    singularity:
        "docker://reslp/iqtree:2.0rc2"
    threads:
        config["iqtree"]["threads"]
    shell:
        """
        rm -rf results/{params.prefix}/phylogeny/concatenated/algn
        cd results/{params.prefix}/phylogeny/concatenated
        mkdir algn
        cp {params.wd}/{input.trims}*trimmed algn
        iqtree -p algn/ --prefix concat -bb 1000 -nt AUTO -m WAG -redo -T {threads} --no-terrace
        rm -r algn
        cd {params.wd}
        touch {output.checkpoint}
        """

rule iqtree_gene_trees:
    input:
        trims = rules.trim.output.dir,
        checkpoint = rules.trim.output.checkpoint,
        ids_file = expand("data/{pre}/ids.txt", pre=config["prefix"])
    output:
        checkpoint = expand("results/{pre}/checkpoints/iqtree_gene_trees.done", pre=config["prefix"]),
        trees = expand("results/{pre}/phylogeny/trees/loci.treefile", pre=config["prefix"])
    params:
        wd = os.getcwd(),
        prefix = config["prefix"]
    threads:
        config["iqtree"]["threads"]
    singularity:
        "docker://reslp/iqtree:2.0rc2"
    shell:
        """
        rm -rf results/{params.prefix}/phylogeny/trees/algn
        mkdir results/{params.prefix}/phylogeny/trees/algn
        cd results/{params.prefix}/phylogeny/trees
        cp {params.wd}/{input.trims}*trimmed algn
        iqtree -S algn --prefix loci -nt AUTO -m WAG -redo -T {threads}
        cd {params.wd}
        touch {output.checkpoint}
        """

rule iqtree_gene_concordance:
    input:
        loci = rules.iqtree_gene_trees.output.trees,
        checkpoint1 = rules.iqtree_gene_trees.output.checkpoint,
        heckpoint2 = rules.iqtree_concat.output.checkpoint,
        concat = rules.iqtree_concat.output.tree,
        alignments = rules.trim.output.dir
    output:
        checkpoint = expand("results/{pre}/checkpoints/iqtree_gene_concordance.done", pre=config["prefix"]),
        treefile = expand("results/{pre}/phylogeny/{pre}.cf.tree", pre=config["prefix"])
    params:
        wd = os.getcwd(),
        prefix = config["prefix"]
    singularity:
        "docker://reslp/iqtree:2.0rc2"
    threads:
        config["iqtree"]["threads"]
    shell:
        """
        cd results/{params.prefix}/phylogeny
        mkdir -p algn
        cp {params.wd}/{input.alignments}*trimmed algn
        iqtree -t {params.wd}/{input.concat} --no-terrace --gcf {params.wd}/{input.loci} -p algn --scf 100 --prefix {params.prefix} -T {threads} -redo
        rm -r algn
        cd {params.wd}
        touch {output.checkpoint}
        """

rule reroot_tree:
    input:
        tree = rules.iqtree_concat.output.tree,
        checkpoint = rules.iqtree_concat.output.checkpoint
    output:
        treefile = expand("results/{pre}/phylogeny/{pre}_concat_rerooted.tre", pre=config["prefix"]),
        checkpoint = expand("results/{pre}/checkpoints/reroot_tree.done", pre=config["prefix"])
    conda:
        "envs/rreroot.yml"
    params:
        wd = os.getcwd(),
        root = config["root"]
    shell:
        """
        Rscript bin/reroot_tree.R {params.wd}/{input.tree} {params.root} {output.treefile}
        touch {output.checkpoint}
        """

rule create_r8s_controlfile:
    input:
        template = expand("data/{pre}/r8s_template.nex", pre=config["prefix"]),
        tree = rules.reroot_tree.output.treefile,
        checkpoint = rules.reroot_tree.output.checkpoint
    output:
        r8s_controlfile = expand("results/{pre}/phylogeny/{pre}_r8s_controlfile.nex", pre=config["prefix"]),
        checkpoint = expand("results/{pre}/checkpoints/create_r8s_controlfile.done", pre=config["prefix"])
    params:
        prefix = config["prefix"]
    conda:
        "envs/phylogenomics.yml"
    shell:
        """
        python bin/create_r8s.py -t {input.template} -tr {input.tree} -l results/{params.prefix}/phylogeny/concatenated/concat.log > {output.r8s_controlfile}
        touch {output.checkpoint}
        """

rule run_r8s:
    input:
        r8s_controlfile = rules.create_r8s_controlfile.output.r8s_controlfile,
        checkpoint = rules.create_r8s_controlfile.output.checkpoint
    output:
        r8s = expand("results/{pre}/phylogeny/{pre}_r8s_output.txt", pre=config["prefix"]),
        checkpoint = expand("results/{pre}/checkpoints/run_r8s.done", pre=config["prefix"])
    singularity:
        "docker://reslp/r8s:1.81"
    shell:
        """
        # this is a bit strange but r8s does not go well with snakemake, so it has to be done like this
        set +euo pipefail
        out=$(r8s -f {input.r8s_controlfile})
        echo $out > {output.r8s}
        touch {output.checkpoint}
        """

rule extract_tree:
    input:
        r8s = rules.run_r8s.output.r8s,
        checkpoint = rules.run_r8s.output.checkpoint
    output:
        ultra_tree = expand("results/{pre}/phylogeny/{pre}_ultra.tre", pre=config["prefix"]),
        checkpoint = expand("results/{pre}/checkpoints/extract_tree.done", pre=config["prefix"])
    shell:
        """
        python bin/extract_tree_from_r8s.py -r {input.r8s} > {output.ultra_tree}
        touch {output.checkpoint}
        """

rule cazy_anc_summary:
    input:
        cazy_summary = expand("data/{pre}/CAZyme.summary.results.csv", pre=config["prefix"]),
        ultra_tree = rules.extract_tree.output.ultra_tree,
        checkpoint = rules.extract_tree.output.checkpoint
    output:
        all_out = expand("results/{pre}/cazy_ancestral_states_summary/all_out.pdf", pre = config["prefix"]),
        overview_heatmap = expand("results/{pre}/cazy_ancestral_states_summary/overview_heatmap.pdf", pre = config["prefix"]),
        checkpoint = expand("results/{pre}/checkpoints/cazy_anc_summary.done", pre=config["prefix"])
    params:
        wd = os.getcwd(),
        prefix = config["prefix"]
    conda:
        "envs/rreroot.yml"
    log:
        "log/output_cazy_anc_summary.log"
    shell:
        """
        if [[ ! -d results/{params.prefix}/cazy_ancestral_states_summary ]]
        then
            mkdir results/{params.prefix}/cazy_ancestral_states_summary
        fi
        cd results/{params.prefix}/cazy_ancestral_states_summary
        Rscript {params.wd}/bin/snake_anc_cazy_xylo_summary.R {params.wd}/results/{params.prefix}/cazy_ancestral_states_summary {params.wd}/{input.cazy_summary} {params.wd}/{input.ultra_tree} {params.prefix} #2>&1 > {params.wd}/{log}
        cd {params.wd}
        touch {output.checkpoint}
        """

rule ancestral_states_cazy:
    input:
        csv = expand("data/{pre}/CAZyme.all.results.csv",pre = config["prefix"]),
        tree = rules.extract_tree.output.ultra_tree
    output:
        dir = directory(expand("results/{pre}/cazy_ancestral_states_genesets/", pre = config["prefix"])),
        checkpoint = expand("results/{pre}/checkpoints/ancestral_states_cazy.done", pre=config["prefix"])
    params:
        wd = os.getcwd(),
        prefix = config["prefix"],
        set = config["set"],
	outprefix = "cazy_ancestral_states_genesets"
    conda:
        "envs/rreroot.yml"
    shell:
        """
        for i in {params.set}
        do
            echo $i
            Rscript bin/snake_anc_cazy_xylo_all.R {params.wd} {input.csv} {input.tree} {params.prefix} $i {params.outprefix}
        done
        touch {output.checkpoint}
        """

rule ancestral_states_all_cazy:
    input:
        csv = expand("data/{pre}/CAZyme.all.results.csv",pre = config["prefix"]),
        tree = rules.extract_tree.output.ultra_tree
    output:
        dir = directory(expand("results/{pre}/cazy_ancestral_states_all_cazy/", pre = config["prefix"])),
        checkpoint = expand("results/{pre}/checkpoints/ancestral_states_cazy_all.done", pre=config["prefix"])
    params:
        wd = os.getcwd(),
        prefix = config["prefix"],
	outprefix = "cazy_ancestral_states_all_cazy"
    conda:
        "envs/rreroot.yml"
    shell:
        """
        Rscript bin/snake_anc_cazy_xylo_all.R {params.wd} {input.csv} {input.tree} {params.prefix} all {params.outprefix}
        touch {output.checkpoint}
        """
rule plot_ancestral_states_cazy_all:
	input:
		rdata = expand("results/{pre}/cazy_ancestral_states_all_cazy/{pre}_all_anc_cazy_all.RData", pre = config["prefix"]),
		checkpoint = rules.ancestral_states_all_cazy.output.checkpoint
	output:
		plot = expand("results/{pre}/plot_ancestral_states_all_cazy/{pre}_all_cazymes.pdf", pre = config["prefix"]),
		checkpoint = expand("results/{pre}/checkpoints/plot_ancestral_states_cazy_all.done", pre=config["prefix"])
	params:
		wd = os.getcwd(),
		prefix = config["prefix"]
	conda:
		"envs/rreroot.yml"
	shell:
		"""
		Rscript bin/plot_all_cazy_anc.R {params.wd} {input.rdata} {params.prefix}
		touch {output.checkpoint}
		""" 
rule similarity_clustering:
    input:
        cazy_data = expand("data/{pre}/CAZyme.all.results.csv",pre = config["prefix"]),
        iprscan_data = expand("data/{pre}/interproscan.results.csv",pre = config["prefix"]),
        pfam_data = expand("data/{pre}/pfam.results.csv",pre = config["prefix"])
    output:
        cazy_clustering = expand("results/{pre}/similarity_clustering/cazy_clustering_all.pdf", pre=config["prefix"]),
        pfam_clustering = expand("results/{pre}/similarity_clustering/pfam_clustering_all.pdf", pre=config["prefix"]),
        iprscan_clustering = expand("results/{pre}/similarity_clustering/interpro_clustering_all.pdf", pre=config["prefix"]),
        checkpoint = expand("results/{pre}/checkpoints/similarity_clustering.done", pre=config["prefix"])
    params:
        wd = os.getcwd(),
        prefix = config["prefix"]
    log:
	"log/simclust.log"
    conda:
        "envs/rreroot.yml"
    shell:
        """
        Rscript bin/snake_clustering.R {params.wd} {input} {params.prefix}
        touch {output.checkpoint}
        """

rule orthology_statistics:
    input:
        directory = rules.infer_orthology.output.dir,
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
        "envs/rorthologystatistics.yml"
    shell:
        """
        for i in {params.orthofiles}
        do
            echo $i
            Rscript bin/ortholog_heatmap.R {output.directory}/{params.prefix}_$i.RData {input.tree} {input.directory}orthofinder/Results_ortho/Comparative_Genomics_Statistics/$i {output.directory}/{params.prefix} $i
        done
        touch {output.checkpoint}
        """

rule astral_species_tree:
        input:
                trees = rules.iqtree_gene_trees.output.trees,
                checkpoint = rules.iqtree_gene_trees.output.checkpoint
        output:
                species_tree = expand("results/{pre}/astral/species_tree.tre", pre=config["prefix"]),
                checkpoint = expand("results/{pre}/checkpoints/astral_species_tree.done", pre=config["prefix"])
        params:
                wd = os.getcwd()
        singularity:
                "docker://reslp/astral:5.7.1"
        log:
                expand("results/{pre}/astral/astral.log", pre=config["prefix"])
        shell:
                """
                java -jar /ASTRAL-5.7.1/Astral/astral.5.7.1.jar -i {input.trees} -o {output.species_tree} 2> {log}
                touch {output.checkpoint}
                """

rule plot_phylogeny:
    input:
        tree = rules.iqtree_gene_concordance.output.treefile,
	astral = rules.astral_species_tree.output.species_tree,
        checkpoint = rules.iqtree_gene_concordance.output.checkpoint
    output:
        phylogeny = expand("results/{pre}/phylogeny/{pre}_phylogeny.pdf", pre = config["prefix"]),
	phylogeny2 = expand("results/{pre}/phylogeny/{pre}_phylogeny2.pdf", pre = config["prefix"]),
        checkpoint = expand("results/{pre}/checkpoints/plot_phylogeny.done", pre=config["prefix"])
    params:
        wd = os.getcwd()
    conda:
        "envs/rorthologystatistics.yml"
    shell:
        """
        Rscript bin/visualize_tree.R {params.wd} {input.tree} {input.astral} {output.phylogeny} {output.phylogeny2}
        touch {output.checkpoint}
        """

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

rule create_gene_family_table:
    input:
        orthodir = rules.infer_orthology.output.dir,
        checkpoint = rules.infer_orthology.output.checkpoint
    output:
        raw_file = expand("results/{pre}/cafe/{pre}_raw_cafe_input_file.tab", pre=config["prefix"]),
        checkpoint = expand("results/{pre}/checkpoints/create_gene_family_table.done", pre=config["prefix"])
    conda:
        "envs/pyutils.yml"
    shell:
        """
        python bin/orthofinder_to_cafe.py -i {input.orthodir}/orthofinder/Results_ortho/Orthogroups/Orthogroups.GeneCount.tsv -o {output.raw_file}
        touch {output.checkpoint}
        """

if config["cafe"]["filter"] == "yes":
	rule filter_cafe_family_table:
		input:
			raw_file = rules.create_gene_family_table.output.raw_file,
			checkpoint = rules.create_gene_family_table.output.checkpoint
		output:
        		filtered_file = expand("results/{pre}/cafe/{pre}_filtered_cafe_input_file.tab", pre=config["prefix"]),
        		checkpoint = expand("results/{pre}/checkpoints/filter_cafe_family_table.done", pre=config["prefix"])
		params:
			wd = os.getcwd(),
			prefix = config["prefix"]
		conda:
			"envs/pyutils.yml"
		shell:
			"""
			cd results/{params.prefix}/cafe/
			name=$(basename "{params.wd}/{input.raw_file}")
			python {params.wd}/bin/clade_and_size_filter.py -i "$name" -o {params.wd}/{output.filtered_file} -s
			cd {params.wd}
			touch {output.checkpoint}
			"""
else:
	rule filter_cafe_family_table:
		input:
			raw_file = rules.create_gene_family_table.output.raw_file,
                	checkpoint = rules.create_gene_family_table.output.checkpoint
		output:
        		filtered_file = expand("results/{pre}/cafe/{pre}_filtered_cafe_input_file.tab", pre=config["prefix"]),
        		checkpoint = expand("results/{pre}/checkpoints/filter_cafe_family_table.done", pre=config["prefix"])
		shell:
			"""
			cp {input.raw_file} {output.filtered_file}
			touch {output.checkpoint}
			"""

rule create_cafe_style_tree:
    input:
        tree = rules.extract_tree.output.ultra_tree,
        checkpoint = rules.extract_tree.output.checkpoint
    output:
        tree = expand("results/{pre}/cafe/{pre}_tree.tre", pre=config["prefix"]),
        cafe_tree = expand("results/{pre}/cafe/{pre}_cafe_tree_single.tre", pre=config["prefix"]),
        checkpoint = expand("results/{pre}/checkpoints/create_cafe_style_tree.done", pre=config["prefix"])
    params:
        wd = os.getcwd()
    conda:
        "envs/rreroot.yml"
    shell:
        """
        Rscript bin/get_cafe_tree.R {params.wd} {input.tree} {output.tree} {output.cafe_tree}
        touch {output.checkpoint}
        """

rule create_cafe_commands:
    input:
        data = rules.filter_cafe_family_table.output.filtered_file,
        template = expand("data/{pre}/cafe_commands_template.sh",pre = config["prefix"]),
        tree = rules.create_cafe_style_tree.output.tree,
        cafe_tree = rules.create_cafe_style_tree.output.cafe_tree
    output:
        cafe_commands = expand("results/{pre}/cafe/{pre}_cafe_commands.sh", pre=config["prefix"]),
        checkpoint = expand("results/{pre}/checkpoints/create_cafe_commands.done", pre=config["prefix"])
    params:
        prefix = config["prefix"]
    conda:
        "envs/pyutils.yml"
    shell:
        """
        python bin/create_cafe.py -t {input.tree} -ct {input.cafe_tree} -template {input.template} -data {input.data} -pre results/{params.prefix}/cafe/{params.prefix} > {output.cafe_commands}
        touch {output.checkpoint}
        """

rule run_cafe:
    input:
        cafe_commands = rules.create_cafe_commands.output.cafe_commands,
        checkpoint = expand("results/{pre}/checkpoints/create_cafe_commands.done", pre=config["prefix"])
    output:
        cafe_results = expand("results/{pre}/cafe/{pre}_cafe_results.cafe", pre=config["prefix"]),
        checkpoint = expand("results/{pre}/checkpoints/run_cafe.done", pre=config["prefix"])
    singularity:
        "docker://reslp/cafe:4.2.1"
    shell:
        """
        cafe {input.cafe_commands}
        touch {output.checkpoint}
        """

rule parse_cafe_output:
    input:
        cafe_results = rules.run_cafe.output.cafe_results,
        checkpoint = rules.run_cafe.output.checkpoint
    output:
        anc = expand("results/{pre}/cafe/{pre}_anc.txt", pre=config["prefix"]),
        fams = expand("results/{pre}/cafe/{pre}_fams.txt", pre=config["prefix"]),
        node = expand("results/{pre}/cafe/{pre}_node.txt", pre=config["prefix"]),
        pub = expand("results/{pre}/cafe/{pre}_pub.txt", pre=config["prefix"]),
        table = expand("results/{pre}/cafe/{pre}_table.txt", pre=config["prefix"]),
        tree = expand("results/{pre}/cafe/{pre}_cafe_labeled_tree.tre", pre=config["prefix"]),
        checkpoint = expand("results/{pre}/checkpoints/parse_cafe_output.done", pre=config["prefix"])
    params:
        prefix = config["prefix"],
        wd = os.getcwd()
    conda:
        "envs/pyutils.yml"
    shell:
        """
        cd results/{params.prefix}/cafe
        python {params.wd}/bin/cafetutorial_report_analysis.py -i {params.wd}/{input.cafe_results} -o {params.prefix}
        cd {params.wd}
        touch {output.checkpoint}
        """

rule visualize_cafe_output:
    input:
        table = rules.parse_cafe_output.output.table,
        tree =  rules.parse_cafe_output.output.tree,
        plottree = rules.extract_tree.output.ultra_tree,
        checkpoint1 = rules.extract_tree.output.checkpoint,
        checkpoint2 = rules.parse_cafe_output.output.checkpoint
    output:
        node_table = expand("results/{pre}/cafe/{pre}_node_number_table.csv", pre=config["prefix"]),
        mapping_file = expand("results/{pre}/cafe/{pre}_cafe_tree_node_mapping.txt", pre=config["prefix"]),
        checkpoint = expand("results/{pre}/checkpoints/visualize_cafe_output.done", pre=config["prefix"])
    params:
        prefix = config["prefix"],
        wd = os.getcwd()
    conda:
        "envs/rreroot.yml"
    shell:
        """
        Rscript bin/visualize_cafe.R {params.wd}/results/{params.prefix}/cafe {input.tree} {input.table} {input.plottree} {params.prefix}
        touch {output.checkpoint}
        """

rule extract_functional_from_cafe:
	input:
		fams = rules.parse_cafe_output.output.fams,
		checkpoint1 = rules.parse_cafe_output.output.checkpoint,
		mapping_file = rules.visualize_cafe_output.output.mapping_file,
		checkpoint2 = rules.visualize_cafe_output.output.checkpoint,
		orthgroups = rules.infer_orthology.output.dir,
		checkpoint3 = rules.infer_orthology.output.checkpoint,
		gff = expand("data/{pre}/{pre}_combined.gff", pre=config["prefix"])
	output:
		exp_cazy = expand("results/{pre}/gene_family_evolution/{pre}_expanded_cazy.txt", pre=config["prefix"]),
		exp_egg = expand("results/{pre}/gene_family_evolution/{pre}_expanded_egg.txt", pre=config["prefix"]),
		exp_go = expand("results/{pre}/gene_family_evolution/{pre}_expanded_go.txt", pre=config["prefix"]),
		exp_interpro = expand("results/{pre}/gene_family_evolution/{pre}_expanded_interpro.txt", pre=config["prefix"]),
		exp_pfam = expand("results/{pre}/gene_family_evolution/{pre}_expanded_pfam.txt", pre=config["prefix"]),
		exp_prod = expand("results/{pre}/gene_family_evolution/{pre}_expanded_prod.txt", pre=config["prefix"]),
		exp_stats = expand("results/{pre}/gene_family_evolution/{pre}_expanded_statistics.txt", pre=config["prefix"]),
		con_cazy = expand("results/{pre}/gene_family_evolution/{pre}_contracted_cazy.txt", pre=config["prefix"]),
		con_egg = expand("results/{pre}/gene_family_evolution/{pre}_contracted_egg.txt", pre=config["prefix"]),
		con_go = expand("results/{pre}/gene_family_evolution/{pre}_contracted_go.txt", pre=config["prefix"]),
		con_interpro = expand("results/{pre}/gene_family_evolution/{pre}_contracted_interpro.txt", pre=config["prefix"]),
		con_pfam = expand("results/{pre}/gene_family_evolution/{pre}_contracted_pfam.txt", pre=config["prefix"]),
		con_prod = expand("results/{pre}/gene_family_evolution/{pre}_contracted_prod.txt", pre=config["prefix"]),
		con_stats = expand("results/{pre}/gene_family_evolution/{pre}_contracted_statistics.txt", pre=config["prefix"]),
		checkpoint = expand("results/{pre}/checkpoints/extract_functional_from_cafe.done", pre=config["prefix"])
	params:
		prefix = config["prefix"]
	conda:
		"envs/pyutils.yml"
	shell:
		"""
		python bin/get_functional_for_cafe_families.py -fam {input.fams} -cm {input.mapping_file} -ortho {input.orthgroups}/orthofinder/Results_ortho/Orthogroups/Orthogroups.tsv -gff {input.gff} -node internal -which + -pre results/{params.prefix}/gene_family_evolution/{params.prefix}
		python bin/get_functional_for_cafe_families.py -fam {input.fams} -cm {input.mapping_file} -ortho {input.orthgroups}/orthofinder/Results_ortho/Orthogroups/Orthogroups.tsv -gff {input.gff} -node internal -which - -pre results/{params.prefix}/gene_family_evolution/{params.prefix}
		touch {output.checkpoint}
		"""

rule visualize_function_annotations_cafe:
	input:
		tree = 	rules.extract_tree.output.ultra_tree,
		checkpoint = rules.extract_tree.output.checkpoint,
		iprmapping = expand("data/{pre}/entry.list", pre=config["prefix"]),
		pfammapping = expand("data/{pre}/pfam_all_description.txt", pre=config["prefix"]),
		files = remove_donefile(expand(rules.extract_functional_from_cafe.output))
	output:
		dir = directory(expand("results/{pre}/cafe_functional_plots/", pre = config["prefix"])),
		checkpoint = expand("results/{pre}/checkpoints/visualize_function_annotations_cafe.done", pre=config["prefix"])
	params:
		wd = os.getcwd(),
		prefix = config["prefix"]
	conda:
		"envs/rreroot.yml"
	shell:
		"""
		for file in {input.files}
		do
		Rscript bin/visualize_cafe_annotation.R {params.wd} {input.tree} $file {params.prefix} {input.iprmapping} {input.pfammapping} {output.dir}
		done
		touch {output.checkpoint}
		"""
rule create_codon_alignments:
	input:
		nucl_file = expand("data/{pre}/{pre}_transcripts.fas", pre=config["prefix"]),
		aa_alignments = rules.align_aa.output.dir,
		checkpoint = rules.align_aa.output.checkpoint,
		orth_dir = rules.infer_orthology.output.dir,
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
