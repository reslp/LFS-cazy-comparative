# A containerized snakemake comparative genomics pipeline

This snakemake pipeline implements comparative genomics analyses (originally developed and currently heavly focused on) lichen-forming fungi. It performs orthology detection, phylogenomics (gene-trees, species tree and concatenated analyses), gene-family evolution analysis (using cafe), ancestral state estimation (for CAZymes) and more. It heavily uses containerized software (using Singularity) and anaconda enviroments.



## Prerequisites
The pipline was designed in such a way that it can run desktop computers (although this is discouraged), solitary linux servers or large HPC clusters. As a result, needed requirements depend.

Local computer or solitary server:

- Linux or MacOS operating system
- Docker (with the possibility to run in privileged mode)
If Docker is not available:
- globally installed singularity 3.4.1+ 
- installed snakemake 5.10.0+ (eg. in an anaconda environment)

On a cluster:

- installed snakemake 5.10.0+ (eg. in an anaconda environment)
- globally installed singularity 3.4.1+
- SGE or SLURM job scheduling system


## Known issues:

### Conda and Singularity combined don't work on VSC4 (found workaround)

On VSC4 originally, the pipeline failed. It looks like the conda environments could not be created successfully. I created an [issue](https://github.com/snakemake/snakemake/issues/304) about this on the snakemake Github page. The error only seems to occur when --use-conda and --use-singularity is combined. Also on sauron this does not happen. After lots of additional investigation I found a workaround. It works by creating an environment variable which points to a directory to store the downloaded conda packages. Othwerise the packages will be downloaded in the global /tmp directory which seems to result in symlinks being not correctly resolved in `.snakemake/conda`. This is currently only done on SLURM clusters and not on sauron. I don't understand exactly why this happens and I am also not 100% sure if this is a problem of snakemake, singularity or anaconda.

## Rulegraph
*Note:* The PDF of the rulegraph is more up-to-date than this png:
<img src="https://github.com/reslp/smsi-comparative/blob/master/rulegraph.png" height="500">

## Configuration

For the pipeline to run you will have to change a few things:

1. Adjust the `data/cluster-config.yaml` file to fit your own cluster enviroment (if you work on a HPC cluster).
2. Change `data/config.yaml` to change parts of the analyses.
3. Place input files in the correct folder. The prefix of the filenames has to correspond to the `data/config.yaml` file.
For example if `config.yaml`has a specified prefix 80_genomes. The data folder should contain a folder called `80_genomes` with the following files:

```
$ ls -lah
insgesamt 2,2G
drwxr-xr-x 4 reslp p71312    17 29. Apr 10:53 .
drwxr-xr-x 4 reslp p71312     6 29. Apr 09:48 ..
-rw-r--r-- 1 reslp p71312  679M 24. Apr 17:46 82_genomes_combined.gff
-rw-r--r-- 1 reslp p71312  404M 24. Apr 17:46 82_genomes_proteins.fas
-rw-r--r-- 1 reslp p71312  1,2G 24. Apr 17:46 82_genomes_transcripts.fas
-rwxr-xr-x 1 reslp p71312  1,4K 24. Apr 17:47 cafe_commands_template.sh
-rw-r--r-- 1 reslp p71312  111K 24. Apr 17:46 CAZyme.all.results.csv
-rw-r--r-- 1 reslp p71312  4,1K 24. Apr 17:46 CAZyme.summary.results.csv
-rw-r--r-- 1 reslp p71312  9,9K 29. Apr 10:49 COGS.all.results.csv
-rw------- 1 reslp p71312  1,9M 24. Apr 17:47 entry.list
-rw-r--r-- 1 reslp p71312   12K 29. Apr 09:31 genome.stats.summary.csv
drwxr-xr-x 2 reslp p71312    83 24. Apr 17:46 gff_files
-rw-r--r-- 1 reslp p71312  2,3K 28. Apr 17:19 ids.txt
-rw-r--r-- 1 reslp p71312  1,8M 24. Apr 17:46 interproscan.results.csv
-rw-r--r-- 1 reslp p71312  703K 24. Apr 17:47 pfam_all_description.txt
-rw-r--r-- 1 reslp p71312 1019K 24. Apr 17:46 pfam.results.csv
drwxr-xr-x 2 reslp p71312    83 24. Apr 17:46 protein_files
-rw-r--r-- 1 reslp p71312   274 24. Apr 17:47 r8s_template.nex
-rw-r--r-- 1 reslp p71312  2,5K 29. Apr 09:31 SM.summary.results.csv
```
The folder `protein_files` contains protein files for each species included in the analysis.

mafft creates tmp files in the folder `/usertmp`. This is a nonstandard binpoint in singularity. Therefore I pass the `tmp` directory as `/usertmp` to singularity as snakemake is called: `--singularity-args -B /tmp:/usertmp`.

## Run the pipeline:

On a local computer (or server) the pipeline can be run by executing this command in the project folder:

```
snakemake --use-conda --use-singularity 
```

In an HPC environment the pipeline is run as follows. It needs an extra submission script to submit all jobs at once. Again this command is run in the project folder:

```
./submit.sh -c data/cluster_config.yaml -s "--use-conda"
```

Also on a cluster the pipeline could be invoked like on a single server, however this is usually not desired.




 
