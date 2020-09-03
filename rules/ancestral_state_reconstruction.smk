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
		cazy_data = expand("data/{pre}/CAZyme.all.results.csv",pre = config["prefix"])
	output:
		phylosig_families = expand("results/{pre}/phylosig/phylosig_sign_cazymes.txt", pre=config["prefix"]),
		checkpoint = expand("results/{pre}/checkpoints/phylosig.done", pre=config["prefix"])
	params:
		wd = os.getcwd(),
		prefix = config["prefix"]
	log:
		"log/phylosig.log"
	conda:
		"../envs/rreroot.yml"
	shell:
		"""
		Rscript bin/phylosig.R {input.tree} {input.cazy_data} {params.wd}/results/{params.prefix}/phylosig/ &> {log}
		touch {output.checkpoint}
		"""


rule ancestral_states_all_cazy:
    input:
        csv = expand("data/{pre}/CAZyme.all.results.csv",pre = config["prefix"]),
        tree = rules.extract_tree.output.ultra_tree
    output:
        #dir = directory(expand("results/{pre}/cazy_ancestral_states_all_cazy/", pre = config["prefix"])),
        checkpoint = expand("results/{pre}/checkpoints/ancestral_states_cazy_all.done", pre=config["prefix"]),
	rdata = expand("results/{pre}/cazy_ancestral_states_all_cazy/{pre}_anc_cazy_all.RData", pre=config["prefix"])
    params:
        wd = os.getcwd(),
        prefix = config["prefix"],
	outprefix = "cazy_ancestral_states_all_cazy"
    conda:
        "../envs/rreroot.yml"
    shell:
        """
        Rscript bin/ancestral_state_reconstruction.R {params.wd} {input.csv} {input.tree} {params.prefix} {params.outprefix}
        touch {output.checkpoint}
        """

rule plot_ancestral_states:
	input:
		rdata = rules.ancestral_states_all_cazy.output.rdata
	output:
		checkpoint = expand("results/{pre}/checkpoints/plot_ancestral_states.done", pre=config["prefix"]),
		rdata = expand("results/{pre}/plots_ancestral_states/ancestral_states.rData", pre=config["prefix"]),
	params:
		wd = os.getcwd(),
		prefix = config["prefix"],
	conda:
		"../envs/rreroot.yml"
	shell:
		"""
		Rscript bin/plot_ancestral_states.R {params.wd} {params.prefix} results/{params.prefix}/plots_ancestral_states/ {input.rdata}
		touch {output.checkpoint}		
		"""
rule pca:
	input:
		cazy_data = rules.ancestral_states_all_cazy.output.rdata,
		check1 = rules.ancestral_states_all_cazy.output.checkpoint,
		phylosig_cazy = rules.phylosig.output.phylosig_families,
		genome_stats = expand("data/{pre}/stats_genomes.csv",pre = config["prefix"])
	output:
		checkpoint = expand("results/{pre}/checkpoints/pca.done", pre = config["prefix"]),
		dir = directory(expand("results/{pre}/pca/", pre = config["prefix"]))
	params:
		wd = os.getcwd(),
		prefix = config["prefix"]
	conda:
		"../envs/rreroot.yml"
	shell:
		"""
		Rscript bin/phyl_pca.R {input.cazy_data} {input.phylosig_cazy} {input.genome_stats} {output.dir}
		touch {output.checkpoint}
		"""

# this rule should be incorporated into the other rule, which should plot the heatmap for all cazymes
rule plot_ancestral_states_cazy_all:
	input:
		checkpoint = rules.ancestral_states_all_cazy.output.checkpoint,
		rdata = rules.ancestral_states_all_cazy.output.rdata
	output:
		plot = expand("results/{pre}/plot_ancestral_states_all_cazy/{pre}_all_cazymes.pdf", pre = config["prefix"]),
		checkpoint = expand("results/{pre}/checkpoints/plot_ancestral_states_cazy_all.done", pre=config["prefix"])
	params:
		wd = os.getcwd(),
		prefix = config["prefix"]
	conda:
		"../envs/rreroot.yml"
	shell:
		"""
		#rdata={params.wd}/results/{params.prefix}/cazy_ancestral_states_all_cazy/{params.prefix}_all_anc_cazy_all.RData
		Rscript bin/plot_all_cazy_anc.R {params.wd} {input.rdata} {params.prefix}
		touch {output.checkpoint}
		""" 
