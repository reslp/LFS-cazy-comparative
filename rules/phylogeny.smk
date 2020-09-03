rule align_aa:
    input:
        dir = rules.rename_ortholog_sequences.output.dir,
        checkpoint = rules.rename_ortholog_sequences.output.checkpoint
    log:
        "log/align_aa.log"
    output:
        checkpoint = expand("results/{pre}/checkpoints/align_aa.done", pre=config["prefix"]),
        dir = directory(expand("results/{pre}/aa_alignments", pre=config["prefix"]))
    conda:
        "../envs/phylogenomics.yml"
    shell:
        """
	export TMPDIR="$(pwd)/tmp"
	if [[ ! -d {output.dir} ]]; then
		mkdir {output.dir}
	fi
        for file in $(find {input.dir}/aa -type f -name "*_renamed");
        do
            outname=$(basename $file)
            echo "Aligning sequences: "$outname
            mafft --quiet --auto {input.dir}/aa/$outname > {output.dir}/$outname"_aligned"
        done
        touch {output.checkpoint}
        """

rule trim:
    input:
        dir = rules.align_aa.output.dir,
        checkpoint = rules.align_aa.output.checkpoint
    output:
        checkpoint = expand("results/{pre}/checkpoints/trim.done", pre=config["prefix"]),
        dir = directory(expand("results/{pre}/aa_trimmed", pre=config["prefix"]))
    conda:
        "../envs/phylogenomics.yml"
    shell:
        """
	export TMPDIR="$(pwd)/tmp"
        if [[ ! -d {output.dir} ]]; then
                mkdir {output.dir}
        fi
        
	for file in $(find {input.dir}/ -type f -name "*_aligned");
        do
            echo "Trimming file: "$file
            outname=$(basename $file)
            trimal -gappyout -in $file -out {output.dir}/$outname"_trimmed"
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
        cp {params.wd}/{input.trims}/*trimmed algn
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
        cp {params.wd}/{input.trims}/*trimmed algn
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
        cp {params.wd}/{input.alignments}/*trimmed algn
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
        "../envs/rreroot.yml"
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
        "../envs/phylogenomics.yml"
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
        "../envs/rorthologystatistics.yml"
    shell:
        """
        Rscript bin/visualize_tree.R {params.wd} {input.tree} {input.astral} {output.phylogeny} {output.phylogeny2}
        touch {output.checkpoint}
        """
