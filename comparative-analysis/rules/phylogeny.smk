if config["phylogeny"]["precalculated"] == "no":
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
            config["threads"]["iqtree_concat"]
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
else:
	rule iqtree_concat:
		input:
			config["phylogeny"]["concatenated"]
		output:
			tree = "results/phylogeny/concatenated/concat.treefile",
			checkpoint = "results/checkpoints/iqtree_concat.done"
		shell:
			"""
			cp $(dirname {input})/* $(dirname {output.tree})
			touch {output.checkpoint}
			"""

	rule iqtree_gene_trees:
		input:
			config["phylogeny"]["gene_trees"]
		output:
			trees = "results/phylogeny/trees/loci.treefile",
			checkpoint = "results/checkpoints/iqtree_gene_trees.done"
		shell:
			"""
			cp {input} {output.trees}
			touch {output.checkpoint}
			"""
	rule iqtree_gene_concordance:
		input:
			loci = rules.iqtree_gene_trees.output.trees,
			checkpoint1 = rules.iqtree_gene_trees.output.checkpoint,
			heckpoint2 = rules.iqtree_concat.output.checkpoint,
			concat = rules.iqtree_concat.output.tree
		output:
			checkpoint = "results/checkpoints/iqtree_gene_concordance.done"
		params:
			wd = os.getcwd(),
			prefix = config["prefix"],
			alignments = config["phylogeny"]["alignments"]
		singularity:
			"docker://reslp/iqtree:2.0rc2"
		threads:
			config["threads"]["iqtree_concat"]
		shell:
			"""
			cd results/phylogeny
			mkdir -p algn
			cp {params.wd}/{params.alignments}*trimmed.fas algn
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
        treefile = "results/phylogeny/concat_rerooted.tre",
        checkpoint = "results/checkpoints/reroot_tree.done"
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
        tree = rules.reroot_tree.output.treefile,
        checkpoint = rules.reroot_tree.output.checkpoint
    output:
        r8s_controlfile = "results/phylogeny/r8s_controlfile.nex",
        checkpoint = "results/checkpoints/create_r8s_controlfile.done"
    params:
        prefix = config["prefix"],
	template = config["r8s_controlfile"]
    conda:
        "../envs/phylogenomics.yml"
    shell:
        """
        python bin/create_r8s.py -t {params.template} -tr {input.tree} -l results/{params.prefix}/phylogeny/concatenated/concat.log > {output.r8s_controlfile}
        touch {output.checkpoint}
        """

rule run_r8s:
    input:
        r8s_controlfile = rules.create_r8s_controlfile.output.r8s_controlfile,
        checkpoint = rules.create_r8s_controlfile.output.checkpoint
    output:
        r8s = "results/phylogeny/r8s_output.txt",
        checkpoint = "results/checkpoints/run_r8s.done"
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
        ultra_tree = "results/phylogeny/phylogeny_r8s_ultra.tre",
        checkpoint = "results/checkpoints/extract_tree.done"
    shell:
        """
        python bin/extract_tree_from_r8s.py -r {input.r8s} > {output.ultra_tree}
        touch {output.checkpoint}
        """
if config["phylogeny"]["precalculated"] == "no":
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
else:
        rule astral_species_tree:
                input:
                        config["phylogeny"]["species_tree"]
                output:
                        species_tree = "results/phylogeny/astral/species_tree.tre",
                        checkpoint = "results/checkpoints/astral_species_tree.done"
                shell:
                        """
                        cp {input} {output.species_tree}
                        touch {output.checkpoint}
                        """

rule plot_phylogeny:
    input:
        checkpoint_tree = rules.iqtree_gene_concordance.output.checkpoint,
	astral = rules.astral_species_tree.output.species_tree,
        checkpoint = rules.iqtree_gene_concordance.output.checkpoint
    output:
        phylogeny = "results/phylogeny/phylogeny.pdf",
	phylogeny2 = "results/phylogeny/phylogeny2.pdf",
        checkpoint = "results/checkpoints/plot_phylogeny.done"
    params:
        wd = os.getcwd(),
	prefix = config["prefix"]
    conda:
        "../envs/rorthologystatistics.yml"
    shell:
        """
        Rscript bin/visualize_tree.R {params.wd} results/phylogeny/{params.prefix}.cf.tree {input.astral} {output.phylogeny} {output.phylogeny2}
        touch {output.checkpoint}
        """
