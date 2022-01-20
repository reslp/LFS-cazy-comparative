# CAZymes in Lichen fungal symbionts

This repository contains the analysis pipeline used to generate the comparative genomic results of the Manuscript: Large differences in carbohydrate degradation and transport potential in the genomes of lichen fungal symbionts. This study is currently available as preprint on BioRxiv: https://www.biorxiv.org/content/10.1101/2021.08.01.454614v1. *If you use parts of this repository please cite us.*

## Contents

[Important information](#Important-information)

[Requirements](#Hardware-and-Software-requirements)

[Installation](#Installation)

[Part 1: Phylogenomics](#Part-1:-Recreate-phylogenomic-analysis)

[Part 2: Genome annotation](#Part-2:-Genome-annotation)

[Part 3: Comparative genomics](#Part-3:-Comparative-genomics)

## Important information

Currently this repository is being worked on to include all files, scripts and streamlined instructions to run all analyses. This will be finished before the manuscript is accepted for publication.

## Hardware and Software requirements
The workflows used to analyse the data are designed in such a way that they can run on desktop computers (although this is strongly discouraged), solitary linux servers or large HPC clusters. We strongly encourage to run this on an HPC system, to be able to benefit from parallelization of different tasks. We support the commonly used job sheduling systems SGE and SLURM.


Local computer or solitary server:

- Linux operating system
- globally installed singularity 3.4.1+ 
- git 2.4+
- installed snakemake 6.0.2+ (eg. installed in an anaconda environment)

On a Linux based HPC cluster:

- installed snakemake 5.19.3+ (eg. installed in an anaconda environment)
- globally installed singularity 3.4.1+
- git 2.4+
- SGE or SLURM job scheduling system
- installed snakemake 6.0.2+ (eg. installed in an anaconda environment)


## Installation

Follow these instructions carefully to install all necessary software and databases used in this study:

1. Get the code in this repository:

```
git clone --recursive https://github.com/reslp/LFS-cazy-comparative.git
```

2. Download external dependecies

Before the setup can proceed, several external dependencies need to be downloaded manually. This cannot be done automatically, since these software packages require you to acquire a personalized license.
Place all downloaded files in `input-data` and make sure that they are correctly named (see below). Otherwise the setup script will complain. 

*GeneMark ES*

GeneMark-ES can be downloaded here: [topaz.gatech.edu/GeneMark/license_download.cgi](http://topaz.gatech.edu/GeneMark/license_download.cgi)

You will get two files: `gmes_linux_64.tar.gz` and `gm_key_64.gz`. Place them into the folder `input-files` so that the setup script is able to find them.

The used version in the paper is 4.62.


*Signal-P*

Download Signal-P from [https://services.healthtech.dtu.dk/service.php?SignalP-4.1](https://services.healthtech.dtu.dk/service.php?SignalP-4.1).

Place the downloaded file into the directory `input-files` and make sure it is named: `signalp-4.1g.Linux.tar.gz` so that the setup script is able to find it.

The used version in the paper is Signal-P 4.1.

*RepeatMasker Libraries*

ADD


3. Run the `setup.sh` script to download remaining workflows for phylogenomics and genome annotation as well as the input data from NCBI.

```
$ ./setup.sh
```

**Currently the setup script is under development and we have to wait for NCBI accession number to finalize it. The script currently downloads all software but not the input data. Downloading input files will be added before the paper gets published**


## Running the workflows to recreate results

### Part 1: Recreate phylogenomic analysis

We used our [phylociraptor](https://github.com/reslp/phylociraptor) pipeline to calculate phylogenomic trees. This part of the analysis can be found in `./phylogeny/phylociraptor`. 
Make sure that you run the `setup.sh` script before.
This script will download the phylociraptor GitHub repository to `./phylogeny` and revert it to the version (based on git commits) that was used in the paper.
Additionally, it will copy the input files necessary to run phylociraptor from `./input-files/phylociraptor` so that the analysis is ready to go. 

To recreate phylogenomic results you should run these commands in this order. Mind you that this can take significant amount of time depending on you computation resources. It is highly recommended to run this on HPC clusters only. The example commands here assume a SLURM submission system.

**IMORTANT: All commands need to be run from within the `phylogeny/phylociraptor` directory. Snakemake and Singularity (see above) need to available** 

1. Setup the phylogenomic analysis:

```
./phylociraptor --setup
```

2. Run orthology detection with BUSCO:

```
./phylociraptor -m busco -t slurm -c data/cluster-config-SLURM.yaml.template
```

4. Create multiple-sequence alignments and filter them:

```
./phylociraptor -m align -t slurm -c data/cluster-config-SLURM.yaml.template
```

5. Perform substitution model testing:

```
./phylociraptor -m model -t slurm -c data/cluster-config-SLURM.yaml.template
```

6. Calculate gene-trees and a species tree:
```
./phylociraptor -m speciestree -t slurm -c data/cluster-config-SLURM.yaml.template
```

7. Calculate a concatenated phylogeny:
```
./phylociraptor -m tree -t slurm -c data/cluster-config-SLURM.yaml.template
```

**IMPORTANT: If you compare these commands to how phylociraptor is currently run you may notice a few differences. This is because phylociraptor has been greatly enhanced since it was used for the LFS-cazy study. You may also want to try using a more recent version**

After you finsihed running step 7 we have recreated the phylogenomic results presented in the LFS-cazy paper.

### Part 2: Genome annotation

For genome annotation we used [funannotate](https://github.com/nextgenusfs/funannotate) using our in-house pipeline [smsi-funannotate](https://github.com/reslp/smsi-funannotate). Running this part of the analysis is independent of Part 1, however you need to make sure to run `setup.sh` and that you have snakemake and singularity installed correctly.

Due to many dependencies and different databases used this part requires additional setup steps. What is described here is largely according to the README files and instructions also given in the [smsi-funannotate](https://github.com/reslp/smsi-funannotate) README file.

**IMPORTANT: All commands need to be run in the `genome-annotation/smsi-funannotate` folder and snakemake and singularity need to be available**

1. Setup funannotate databases:
```
$ singularity shell -B $(pwd)/data/database:/data/database docker://reslp/funannotate:1.7.4
$ funannotate setup -i all
$ funannotate database
```

You should now see a list of all installed funannotate databases. For this study, the database looked like this:

```
  Database          Type        Version      Date         Num_Records   Md5checksum
  merops            diamond     12.0         2017-10-04          5009   a6dd76907896708f3ca5335f58560356
  uniprot           diamond     2020_06      2020-12-02        563972   19435273583d04a417e24aa158f5f541
  dbCAN             hmmer3      9.0          2020-08-04           641   04696dfba1c3bb82ff9b72cfbb3e4a65
  pfam              hmmer3      33.1         2020-04            18259   228db640818af54066a1c23404a3ba38
  repeats           diamond     1.0          2021-01-25         11950   4e8cafc3eea47ec7ba505bb1e3465d21
  go                text        2021-01-01   2021-01-01         47198   24c95428e037fdc4d5899414a8e097e8
  mibig             diamond     1.4          2021-01-25         31023   118f2c11edde36c81bdea030a0228492
  interpro          xml         83.0         2020-12-03         38345   07904253628e12e5aeb2ab1620ba8b99
  busco_outgroups   outgroups   1.0          2021-01-25             8   6795b1d4545850a4226829c7ae8ef058
  gene2product      text        1.65         2020-10-05         33749   0f2d46d3f7d1f08a87c1d4fec1eb05bb
```

**IMPORTANT: These databases are constantly updated, so it may happen that you will see different version numbers, Md5checksum or Dates in your own databases. These differences may lead to slight differences during gene-calling and annotation compared to what is reported in the published version of the paper**

2. Setup additional dependecies of funannotate which require a license

Funannotate uses several additional programs such as eggnog-mapper, Signal-P and GeneMark ES, which require you to get a license. For eggnog-mapper we created a Docker container (which can be used with Singularity as well).

 






1. 


