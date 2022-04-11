
rule create_gene_family_table:
    input:
        checkpoint = rules.infer_orthology.output.checkpoint,
        cazy_file = "data/83_genomes/CAZyme.all.results.csv" 
    output:
        raw_file_cz = "results/gene_family_evolution/cafe/cazy_raw_cafe_input_file.tab"
    singularity:
        "docker://reslp/biopython_plus:1.77"
    params:
        prefix = config["prefix"]
    shell:
        """
        python bin/cazy_counts_to_cafe.py -i {input.cazy_file} -o {output.raw_file_cz}
	"""


#if config["cafe"]["filter"] == "yes":
#	rule filter_cafe_family_table:
#		input:
#			raw_file = rules.create_gene_family_table.output.raw_file_of
#		output:
#			orthofinder_filtered_file = "results/gene_family_evolution/cafe/orthofinder_filtered_cafe_input_file.tab"
#		params:
#			wd = os.getcwd(),
#			prefix = config["prefix"]
#		conda:
#			"../envs/pyutils.yml"
#		shell:
#			"""
#			cd results/gene_family_evolution/cafe/
#			name=$(basename "{params.wd}/{input.raw_file}")
#			python {params.wd}/bin/clade_and_size_filter.py -i "$name" -o {params.wd}/{output.orthofinder_filtered_file} -s
#			cd {params.wd}
#			"""
#else:
#	rule filter_cafe_family_table:
#		input:
#			raw_file = rules.create_gene_family_table.output.raw_file_of
#		output:
#			orthofinder_filtered_file = "results/gene_family_evolution/cafe/orthofinder_filtered_cafe_input_file.tab"
#		shell:
#			"""
#			cp {input.raw_file} {output.filtered_file}
#			"""

rule create_cafe_style_tree:
    input:
        tree = rules.extract_tree.output.ultra_tree,
        checkpoint = rules.extract_tree.output.checkpoint
    output:
        tree = "results/gene_family_evolution/cafe/cafe_tree.tre",
        cafe_tree = "results/gene_family_evolution/cafe/cafe_tree_rates.tre"
    params:
        wd = os.getcwd()
    conda:
        "../envs/rreroot.yml"
    shell:
        """
        Rscript bin/get_cafe_tree.R {params.wd} {input.tree} {output.tree} {output.cafe_tree} Xylographa_vitiligo,Acarospora_spec
        """

#rule create_cafe_commands:
#    input:
#        #data = rules.filter_cafe_family_table.output.filtered_file,
#	data = rules.create_gene_family_table.output.raw_file,
#        template = expand("data/{pre}/cafe_commands_template.sh",pre = config["prefix"]),
#        tree = rules.create_cafe_style_tree.output.tree,
#        cafe_tree = rules.create_cafe_style_tree.output.cafe_tree
#    output:
#        cafe_commands = expand("results/{pre}/cafe/{pre}_cafe_commands.sh", pre=config["prefix"]),
#        checkpoint = expand("results/{pre}/checkpoints/create_cafe_commands.done", pre=config["prefix"])
#    params:
#        prefix = config["prefix"]
#    conda:
#        "../envs/pyutils.yml"
#    shell:
#        """
#        python bin/create_cafe.py -t {input.tree} -ct {input.cafe_tree} -template {input.template} -data {input.data} -pre results/{params.prefix}/cafe/{params.prefix} > {output.cafe_commands}
#        touch {output.checkpoint}
#        """

rule run_cafe:
	input:
		data_cz = rules.create_gene_family_table.output.raw_file_cz,
		tree = rules.create_cafe_style_tree.output.tree,
		cafe_tree = rules.create_cafe_style_tree.output.cafe_tree
	output:
		#cafe_results = "results/gene_family_evolution/cafe/cafe_output/cafe_results.cafe",
		checkpoint = "results/gene_family_evolution/checkpoints/run_cafe.done"
	log:
		"log/cafe.log"
	threads: config["threads"]["run_cafe"]
	singularity: "docker://reslp/cafe:08d27a1"
	shell:
		"""
		echo "Estimate an error model for a single rate":
		cafe -t {input.tree} -i {input.data_cz} -p -o results/gene_family_evolution/cafe/error_model_single_rate -e
		
		echo "Estimate an error model for two rates":
		cafe -t {input.tree} -i {input.data_cz} -y {input.cafe_tree} -p -o results/gene_family_evolution/cafe/error_model_two_rates -e
		
		echo "Run cafe for CAZy data SINGLE RATE:"
		cafe -t {input.tree} -i {input.data_cz} -p -o results/gene_family_evolution/cafe/cafe_output_cazy_single_rate_1 | tee {log}
		cafe -t {input.tree} -i {input.data_cz} -p -o results/gene_family_evolution/cafe/cafe_output_cazy_single_rate_2 | tee {log}
		cafe -t {input.tree} -i {input.data_cz} -p -o results/gene_family_evolution/cafe/cafe_output_cazy_single_rate_3 | tee {log}
		cafe -t {input.tree} -i {input.data_cz} -p -o results/gene_family_evolution/cafe/cafe_output_cazy_single_rate_4 | tee {log}
		cafe -t {input.tree} -i {input.data_cz} -p -o results/gene_family_evolution/cafe/cafe_output_cazy_single_rate_5 | tee {log}
		
		echo "Run cafe for CAZy data SINGLE RATE + error model:"
		cafe -t {input.tree} -i {input.data_cz} -p -o results/gene_family_evolution/cafe/cafe_output_cazy_single_rate_1_errm -eresults/gene_family_evolution/cafe/error_model_single_rate/Base_error_model.txt | tee {log}
		cafe -t {input.tree} -i {input.data_cz} -p -o results/gene_family_evolution/cafe/cafe_output_cazy_single_rate_2_errm -eresults/gene_family_evolution/cafe/error_model_single_rate/Base_error_model.txt | tee {log}
		cafe -t {input.tree} -i {input.data_cz} -p -o results/gene_family_evolution/cafe/cafe_output_cazy_single_rate_3_errm -eresults/gene_family_evolution/cafe/error_model_single_rate/Base_error_model.txt | tee {log}
		cafe -t {input.tree} -i {input.data_cz} -p -o results/gene_family_evolution/cafe/cafe_output_cazy_single_rate_4_errm -eresults/gene_family_evolution/cafe/error_model_single_rate/Base_error_model.txt | tee {log}
		cafe -t {input.tree} -i {input.data_cz} -p -o results/gene_family_evolution/cafe/cafe_output_cazy_single_rate_5_errm -eresults/gene_family_evolution/cafe/error_model_single_rate/Base_error_model.txt | tee {log}

		echo "Run cafe for CAZy data TWO RATES:"
		cafe -t {input.tree} -i {input.data_cz} -p -y {input.cafe_tree} -o results/gene_family_evolution/cafe/cafe_output_cazy_two_rates_1 | tee {log}
		cafe -t {input.tree} -i {input.data_cz} -p -y {input.cafe_tree} -o results/gene_family_evolution/cafe/cafe_output_cazy_two_rates_2 | tee {log}
		cafe -t {input.tree} -i {input.data_cz} -p -y {input.cafe_tree} -o results/gene_family_evolution/cafe/cafe_output_cazy_two_rates_3 | tee {log}
		cafe -t {input.tree} -i {input.data_cz} -p -y {input.cafe_tree} -o results/gene_family_evolution/cafe/cafe_output_cazy_two_rates_4 | tee {log}
		cafe -t {input.tree} -i {input.data_cz} -p -y {input.cafe_tree} -o results/gene_family_evolution/cafe/cafe_output_cazy_two_rates_5 | tee {log}

		echo "Run cafe for CAZy data TWO RATES + error model"
		cafe -t {input.tree} -i {input.data_cz} -p -y {input.cafe_tree} -o results/gene_family_evolution/cafe/cafe_output_cazy_two_rates_1_errm -eresults/gene_family_evolution/cafe/error_model_two_rates/Base_error_model.txt| tee {log}
		cafe -t {input.tree} -i {input.data_cz} -p -y {input.cafe_tree} -o results/gene_family_evolution/cafe/cafe_output_cazy_two_rates_2_errm -eresults/gene_family_evolution/cafe/error_model_two_rates/Base_error_model.txt | tee {log}
		cafe -t {input.tree} -i {input.data_cz} -p -y {input.cafe_tree} -o results/gene_family_evolution/cafe/cafe_output_cazy_two_rates_3_errm -eresults/gene_family_evolution/cafe/error_model_two_rates/Base_error_model.txt | tee {log}
		cafe -t {input.tree} -i {input.data_cz} -p -y {input.cafe_tree} -o results/gene_family_evolution/cafe/cafe_output_cazy_two_rates_4_errm -eresults/gene_family_evolution/cafe/error_model_two_rates/Base_error_model.txt | tee {log}
		cafe -t {input.tree} -i {input.data_cz} -p -y {input.cafe_tree} -o results/gene_family_evolution/cafe/cafe_output_cazy_two_rates_5_errm -eresults/gene_family_evolution/cafe/error_model_two_rates/Base_error_model.txt | tee {log}
		
		touch {output.checkpoint}
		"""

rule parse_cafe_output:
	input:
		rules.run_cafe.output.checkpoint
	output:
		sign_families = "results/gene_family_evolution/cafe_significant_families_per_model.txt",
		sign_trees = "results/gene_family_evolution/significant_trees.tre"
	params:
		wd = os.getcwd()
	shell:
		"""
		cd results/gene_family_evolution/cafe
		for folder in $(ls -d cafe_output*/); do 
			families=$(grep -P "\ty" $folder/Base_family_results.txt | awk '{{printf($1":"$2";")}}'); 
			outf=$(basename $folder); 
			printf "$outf\t$families\n"; 
			done > {params.wd}/{output.sign_families}
		echo $'#nexus\nbegin trees;' > {params.wd}/{output.sign_trees}
		for folder in $(ls -d cafe_output*/); do
			fd=$(basename $folder)
			export fd
			grep "*" $folder\Base_asr.tre | perl -pe 's/(TREE )/TREE $ENV{{fd}}_/' >> {params.wd}/{output.sign_trees}
			done
		echo "end;">> {params.wd}/{output.sign_trees}
		"""

rule plot_cafe_summary_tree:
	input:
		sign_trees = rules.parse_cafe_output.output.sign_trees
	output:
		cafe_summary_tree = "results/gene_family_evolution/cafe_summary_tree.pdf"
	params:
		wd = os.getcwd()
	singularity:
		"docker://reslp/rphylogenetics:4.0.3"
	shell:
		"""
		Rscript bin/parse_cafe_trees.R {params.wd} {input.sign_trees} {output.cafe_summary_tree}
		"""



#rule parse_cafe_output:
#    input:
#        cafe_results = rules.run_cafe.output.cafe_results,
#        checkpoint = rules.run_cafe.output.checkpoint
#    output:
#        anc = expand("results/{pre}/cafe/{pre}_anc.txt", pre=config["prefix"]),
#        fams = expand("results/{pre}/cafe/{pre}_fams.txt", pre=config["prefix"]),
#        node = expand("results/{pre}/cafe/{pre}_node.txt", pre=config["prefix"]),
#        pub = expand("results/{pre}/cafe/{pre}_pub.txt", pre=config["prefix"]),
#        table = expand("results/{pre}/cafe/{pre}_table.txt", pre=config["prefix"]),
#        tree = expand("results/{pre}/cafe/{pre}_cafe_labeled_tree.tre", pre=config["prefix"]),
#        checkpoint = expand("results/{pre}/checkpoints/parse_cafe_output.done", pre=config["prefix"])
#    params:
#        prefix = config["prefix"],
#        wd = os.getcwd()
#    conda:
#        "../envs/pyutils.yml"
#    shell:
#        """
#        cd results/{params.prefix}/cafe
#        python {params.wd}/bin/cafetutorial_report_analysis.py -i {params.wd}/{input.cafe_results} -o {params.prefix}
#        cd {params.wd}
#        touch {output.checkpoint}
#        """

#rule visualize_cafe_output:
#    input:
#        table = rules.parse_cafe_output.output.table,
#        tree =  rules.parse_cafe_output.output.tree,
#        plottree = rules.extract_tree.output.ultra_tree,
#        checkpoint1 = rules.extract_tree.output.checkpoint,
#        checkpoint2 = rules.parse_cafe_output.output.checkpoint
#    output:
#        node_table = expand("results/{pre}/cafe/{pre}_node_number_table.csv", pre=config["prefix"]),
#        mapping_file = expand("results/{pre}/cafe/{pre}_cafe_tree_node_mapping.txt", pre=config["prefix"]),
#        checkpoint = expand("results/{pre}/checkpoints/visualize_cafe_output.done", pre=config["prefix"])
#    params:
#        prefix = config["prefix"],
#        wd = os.getcwd()
#    conda:
#        "../envs/rreroot.yml"
#    shell:
#        """
#        Rscript bin/visualize_cafe.R {params.wd}/results/{params.prefix}/cafe {input.tree} {input.table} {input.plottree} {params.prefix}
#        touch {output.checkpoint}
#        """
#
#rule extract_functional_from_cafe:
#	input:
#		fams = rules.parse_cafe_output.output.fams,
#		checkpoint1 = rules.parse_cafe_output.output.checkpoint,
#		mapping_file = rules.visualize_cafe_output.output.mapping_file,
#		checkpoint2 = rules.visualize_cafe_output.output.checkpoint,
#		#orthgroups = expand("results/{pre}/orthofinder/", pre=config["prefix"]),
#		checkpoint3 = rules.infer_orthology.output.checkpoint,
#		gff = expand("data/{pre}/{pre}_combined.gff", pre=config["prefix"])
#	output:
#		exp_cazy = expand("results/{pre}/gene_family_evolution/{pre}_expanded_cazy.txt", pre=config["prefix"]),
#		exp_egg = expand("results/{pre}/gene_family_evolution/{pre}_expanded_egg.txt", pre=config["prefix"]),
#		exp_go = expand("results/{pre}/gene_family_evolution/{pre}_expanded_go.txt", pre=config["prefix"]),
#		exp_interpro = expand("results/{pre}/gene_family_evolution/{pre}_expanded_interpro.txt", pre=config["prefix"]),
#		exp_pfam = expand("results/{pre}/gene_family_evolution/{pre}_expanded_pfam.txt", pre=config["prefix"]),
#		exp_prod = expand("results/{pre}/gene_family_evolution/{pre}_expanded_prod.txt", pre=config["prefix"]),
#		exp_stats = expand("results/{pre}/gene_family_evolution/{pre}_expanded_statistics.txt", pre=config["prefix"]),
#		con_cazy = expand("results/{pre}/gene_family_evolution/{pre}_contracted_cazy.txt", pre=config["prefix"]),
#		con_egg = expand("results/{pre}/gene_family_evolution/{pre}_contracted_egg.txt", pre=config["prefix"]),
#		con_go = expand("results/{pre}/gene_family_evolution/{pre}_contracted_go.txt", pre=config["prefix"]),
#		con_interpro = expand("results/{pre}/gene_family_evolution/{pre}_contracted_interpro.txt", pre=config["prefix"]),
#		con_pfam = expand("results/{pre}/gene_family_evolution/{pre}_contracted_pfam.txt", pre=config["prefix"]),
#		con_prod = expand("results/{pre}/gene_family_evolution/{pre}_contracted_prod.txt", pre=config["prefix"]),
#		con_stats = expand("results/{pre}/gene_family_evolution/{pre}_contracted_statistics.txt", pre=config["prefix"]),
#		checkpoint = expand("results/{pre}/checkpoints/extract_functional_from_cafe.done", pre=config["prefix"])
#	params:
#		prefix = config["prefix"]
#	conda:
#		"../envs/pyutils.yml"
#	shell:
#		"""
#		python bin/get_functional_for_cafe_families.py -fam {input.fams} -cm {input.mapping_file} -ortho results/{params.prefix}/orthofinder//orthofinder/Results_ortho/Orthogroups/Orthogroups.tsv -gff {input.gff} -node internal -which + -pre results/{params.prefix}/gene_family_evolution/{params.prefix}
#		python bin/get_functional_for_cafe_families.py -fam {input.fams} -cm {input.mapping_file} -ortho results/{params.prefix}/orthofinder/orthofinder/Results_ortho/Orthogroups/Orthogroups.tsv -gff {input.gff} -node internal -which - -pre results/{params.prefix}/gene_family_evolution/{params.prefix}
#		touch {output.checkpoint}
#		"""
#
#rule visualize_function_annotations_cafe:
#	input:
#		tree = 	rules.extract_tree.output.ultra_tree,
#		checkpoint = rules.extract_tree.output.checkpoint,
#		iprmapping = expand("data/{pre}/entry.list", pre=config["prefix"]),
#		pfammapping = expand("data/{pre}/pfam_all_description.txt", pre=config["prefix"]),
#		files = remove_donefile(expand(rules.extract_functional_from_cafe.output))
#	output:
#		dir = directory(expand("results/{pre}/cafe_functional_plots/", pre = config["prefix"])),
#		checkpoint = expand("results/{pre}/checkpoints/visualize_function_annotations_cafe.done", pre=config["prefix"])
#	params:
#		wd = os.getcwd(),
#		prefix = config["prefix"]
#	conda:
#		"../envs/rreroot.yml"
#	shell:
#		"""
#		for file in {input.files}
#		do
#		Rscript bin/visualize_cafe_annotation.R {params.wd} {input.tree} $file {params.prefix} {input.iprmapping} {input.pfammapping} {output.dir}
#		done
#		touch {output.checkpoint}
#		"""
