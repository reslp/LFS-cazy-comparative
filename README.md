# CAZymes in Lichen fungal symbionts

This repository contains the analysis pipeline used to generate the comparative genomic results of the Manuscript: Large differences in carbohydrate degradation and transport potential in the genomes of lichen fungal symbionts. This study is currently available as preprint on BioRxiv: https://www.biorxiv.org/content/10.1101/2021.08.01.454614v1. *If you use parts of this repository please cite us.*

## Important information

Currently this repository is being worked to include all files and streamlined instructions to run all analyses. This will be finished before the manuscript is accepted for publication in a scientific journal.

## Prerequisites
The workflows used to analyse the data are designed in such a way that they can run on desktop computers (although this is strongly discouraged), solitary linux servers or large HPC clusters. We strongly encourage to run this on an HPC system, to be able to benefit from parallelization of different tasks. We support the commonly used job sheduling systems SGE and SLURM.

## Hardware and Software requirements
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
