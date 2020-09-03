# this rule is currently not used
rule character_correlation:
	input:
		ultra_tree = rules.extract_tree.output.ultra_tree,
		cazy_data = expand("data/{pre}/CAZyme.all.results.csv",pre = config["prefix"]),
		discrete_data = expand("data/{pre}/character_information.csv",pre = config["prefix"]),
	output:
		checkpoint = expand("results/{pre}/checkpoints/character_correlation.done", pre=config["prefix"]),
		outdir = directory(expand("results/{pre}/character_correlation", pre=config["prefix"]))
	params:
		wd = os.getcwd()
	conda:
		"envs/rreroot.yml"
	threads: 8
	shell:
		"""
		if [[ ! -d {output.outdir} ]]
        	then
            		mkdir {output.outdir}
        	fi
		Rscript bin/character_correlation.R {params.wd} {input.ultra_tree} {input.discrete_data} {input.cazy_data} {output.outdir} {threads}
		touch {output.checkpoint}
		"""
