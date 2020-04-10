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


## Rulegraph
<img src="https://github.com/reslp/smsi-comparative/blob/master/rulegraph.png" height="500">

## Configuration

For the pipeline to run you will have to change a few things:

1. Adjust the `data/cluster-config.yaml` file to fit your own cluster enviroment (if you work on a HPC cluster).
2. Change `data/config.yaml` to change parts of the analyses.
3. Place input files in the correct folder. The prefix of the filenames has to correspond to the `data/config.yaml` file.
For example if `config.yaml`has a specified prefix 80_genomes. The data folder should contain a folder called `80_genomes` with the following files:

```
$ ls -lah
insgesamt 2,6G
drwxr-xr-x 3 reslph domain_users   16  6. Apr 19:54 .
drwxr-xr-x 4 reslph domain_users    9  7. Apr 09:20 ..
-rw-r--r-- 1 reslph domain_users 664M  1. Apr 16:10 80_genomes_combined.gff
-rw-r--r-- 1 reslph domain_users 395M  1. Apr 16:10 80_genomes_proetins.fas
-rw-r--r-- 1 reslph domain_users 395M  6. Apr 19:43 80_genomes_proteins.fas
-rw-r--r-- 1 reslph domain_users 1,2G  1. Apr 16:10 80_genomes_transcripst.fas
-rwxr-xr-x 1 reslph domain_users 1,4K 30. Mär 08:55 cafe_commands_template.sh
-rw-r--r-- 1 reslph domain_users 109K  1. Apr 16:10 CAZyme.all.results.csv
-rw-r--r-- 1 reslph domain_users 4,0K  1. Apr 16:10 CAZyme.summary.results.csv
-rw------- 1 reslph domain_users 1,9M 30. Mär 17:30 entry.list
-rw-r--r-- 1 reslph domain_users 2,3K  1. Apr 16:10 ids.txt
-rw-r--r-- 1 reslph domain_users 1,8M  1. Apr 16:10 interproscan.results.csv
-rw-r--r-- 1 reslph domain_users 703K 30. Mär 08:55 pfam_all_description.txt
-rw-r--r-- 1 reslph domain_users 999K  1. Apr 16:10 pfam.results.csv
drwxr-xr-x 2 reslph domain_users   82  1. Apr 16:13 protein_files
-rw-r--r-- 1 reslph domain_users  274 30. Mär 08:55 r8s_template.nex
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




 