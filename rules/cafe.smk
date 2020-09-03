# these rules are currently not used
rule create_gene_family_table:
    input:
        orthodir = expand("results/{pre}/orthofinder/", pre=config["prefix"]),
        checkpoint = rules.infer_orthology.output.checkpoint
    output:
        raw_file = expand("results/{pre}/cafe/{pre}_raw_cafe_input_file.tab", pre=config["prefix"]),
        checkpoint = expand("results/{pre}/checkpoints/create_gene_family_table.done", pre=config["prefix"])
    conda:
        "../envs/pyutils.yml"
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
			"../envs/pyutils.yml"
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
        "../envs/rreroot.yml"
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
        "../envs/pyutils.yml"
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
        "../envs/pyutils.yml"
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
        "../envs/rreroot.yml"
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
		orthgroups = expand("results/{pre}/orthofinder/", pre=config["prefix"]),
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
		"../envs/pyutils.yml"
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
		"../envs/rreroot.yml"
	shell:
		"""
		for file in {input.files}
		do
		Rscript bin/visualize_cafe_annotation.R {params.wd} {input.tree} $file {params.prefix} {input.iprmapping} {input.pfammapping} {output.dir}
		done
		touch {output.checkpoint}
		"""
