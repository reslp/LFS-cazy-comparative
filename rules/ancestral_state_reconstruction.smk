rule cazy_anc_summary:
    input:
        cazy_summary = config["funannotate_input"]["cazy_summary"],
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
        "../envs/rreroot.yml"
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

rule phylosig:
	input:
		tree = rules.extract_tree.output.ultra_tree,
		cazy_data = config["funannotate_input"]["cazy"]
	output:
		phylosig_families = "results/ancestral_states_cazy/phylosig/phylosig_sign_cazymes.txt",
		checkpoint = "results/checkpoints/phylosig.done"
	params:
		wd = os.getcwd(),
		prefix = config["prefix"]
	log:
		"log/phylosig.log"
	conda:
		"../envs/rreroot.yml"
	shell:
		"""
		Rscript bin/phylosig.R {input.tree} {input.cazy_data} {params.wd}/results/ancestral_states_cazy/phylosig/ &> {log}
		touch {output.checkpoint}
		"""


rule ancestral_states_all_cazy:
    input:
        csv = config["funannotate_input"]["cazy"],
        tree = rules.extract_tree.output.ultra_tree
    output:
        checkpoint = "results/checkpoints/ancestral_states_cazy_all.done",
	rdata = "results/ancestral_states_cazy/anc_cazy_all.RData"
    params:
        wd = os.getcwd(),
        prefix = config["prefix"],
	outprefix = "ancestral_states_cazy"
    conda:
        "../envs/rreroot.yml"
    shell:
        """
        Rscript bin/ancestral_state_reconstruction.R {params.wd} {input.csv} {input.tree} {params.prefix} {params.outprefix}
        touch {output.checkpoint}
        """

rule plot_ancestral_states:
	input:
		rdata = rules.ancestral_states_all_cazy.output.rdata,
		phylosig = rules.phylosig.output.checkpoint
	output:
		checkpoint = "results/checkpoints/plot_ancestral_states.done",
		plot_all = "results/ancestral_states_cazy/plots/all_cazymes.pdf",
		rdata_gs = "results/ancestral_states_cazy/plots/ancestral_states.rData",
		rdata_all = "results/ancestral_states_cazy/plots/ancestral_states_all.rData"
	params:
		wd = os.getcwd(),
		prefix = config["prefix"]
	conda:
		"../envs/rreroot.yml"
	shell:
		"""
		# plot the heatmaps and tree for different genesets:
		Rscript bin/plot_ancestral_states.R {params.wd} {params.prefix} results/ancestral_states_cazy/plots/ {input.rdata} {output.rdata_gs}
		# plot heatmap for all cazymes:
		Rscript bin/plot_all_cazy_anc.R {params.wd} {input.rdata} {output.plot_all} {output.rdata_all}
		touch {output.checkpoint}		
		"""
rule pca:
	input:
		cazy_data = rules.ancestral_states_all_cazy.output.rdata,
		check1 = rules.ancestral_states_all_cazy.output.checkpoint,
		apriori_sets = "data/settings_and_parameters/apriori_cazyme_sets.txt",
		genome_stats = config["other_input"]["taxonomy_information"],
		colors_file = "data/settings_and_parameters/color_information.csv"
	output:
		checkpoint = "results/checkpoints/pca.done",
		dir = directory("results/ancestral_states_cazy/pca/")
	params:
		wd = os.getcwd(),
		prefix = config["prefix"]
	conda:
		"../envs/rreroot.yml"
	shell:
		"""
		Rscript bin/phyl_pca.R {input.cazy_data} {input.colors_file} {input.apriori_sets} {input.genome_stats} {output.dir} 
		touch {output.checkpoint}
		"""

# this rule should be incorporated into the other rule, which should plot the heatmap for all cazymes
rule plot_ancestral_states_cazy_all:
	input:
		checkpoint = rules.ancestral_states_all_cazy.output.checkpoint,
		rdata = rules.ancestral_states_all_cazy.output.rdata
	output:
		plot = "results/ancestral_states_cazy/plots/all_cazymes.pdf",
		checkpoint = "results/checkpoints/plot_ancestral_states_cazy_all.done"
	params:
		wd = os.getcwd(),
		prefix = config["prefix"]
	conda:
		"../envs/rreroot.yml"
	shell:
		"""
		#rdata={params.wd}/results/ancestral_states_cazy/plots/{params.prefix}_all_anc_cazy_all.RData
		Rscript bin/plot_all_cazy_anc.R {params.wd} {input.rdata} {output.plot} {params.wd}/results/ancestral_states_cazy/all_anc_cazy_all.RData
		touch {output.checkpoint}
		""" 
