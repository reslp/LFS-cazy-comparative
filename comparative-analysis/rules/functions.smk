import glob

def get_protein_files(wildcards):
	searchstring = config["funannotate_input"]["protein_folder"]+"/*.fa"
	file_list = glob.glob(searchstring)
	return file_list
