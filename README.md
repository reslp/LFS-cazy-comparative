# CAZymes in Lichen fungal symbionts

This repository contains the analysis pipeline used to generate the comparative genomic results of the Manuscript: Large differences in carbohydrate degradation and transport potential in the genomes of lichen fungal symbionts. This study is currently available as preprint on BioRxiv: https://www.biorxiv.org/content/10.1101/2021.08.01.454614v1. *If you use parts of this repository please cite us.*

## Contents

[Important information](#Important-information)
[Requirements](#Hardware-and-Software-requirements)
[Installation](#Installation)
[Part 1: Phylogenomics](#Part-1:-Recreate-phylogenomic-analysis)

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

Installation should take only a few minutes. 

1. Get the code in this repository run the following command:

```
git clone --recursive https://github.com/reslp/LFS-cazy-comparative.git
```

2. Run the `setup.sh` script to download remaining workflows for phylogenomics and genome annotation as well as the input data from NCBI.

**Currently the setup script is under development and we have to wait for NCBI accession number to finalize it. The script currently downloads all software but not the input data. Downloading input files will be added before the paper gets published**


## Running the workflows to recreate results

# Part 1: Recreate phylogenomic analysis

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



 
