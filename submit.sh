#!/bin/bash

set -e

usage() { 
	echo "Welcome to the pipeline submission script. A script helps to submit jobs to SLURM and SGE clusters with snakemake and singularity"
	echo
	echo "Usage: $0 [-v] [-c <cluster_config_file>] [-s <snakemke_args>]" 
	echo
	echo "Options:"
	echo "	-c <cluster_config_file> Path to cluster config file in YAML format (mandatory). "
	echo "	-s <snakemake_args> Additional arguments passed on to the snakemake command (optional). snakemake is run with --immediate-submit -pr --notemp --latency-wait 600 --use-singularity --jobs 1001 by default." 
	1>&2; exit 1; }
	
version() {
	echo "$0 v0.1"
	exit 0
}

while getopts ":v:c:s:" option;
	do
		case "${option}"
		in
			v) version;;
			c) CLUSTER_CONFIG=${OPTARG};;
			s) SM_ARGS=${OPTARG};;
			*) echo "Illegal option --$OPTARG\n" >&2; usage;;
			?) echo "Illegal option --$OPTARG\n" >&2 usage;;
		esac
	done
if [ $OPTIND -eq 1 ]; then usage; fi

CLUSTER=""
command -v qsub >/dev/null 2>&1 && { echo >&2 "SGE detected, will use qsub to submit jobs."; CLUSTER="sge"; }
command -v sbatch >/dev/null 2>&1 && { echo >&2 "SLURM detected, will use sbatch to submit jobs."; CLUSTER="slurm"; }

#snakemake --use-singularity --jobs 1001 --cluster-config data/cluster_config-sauron.yaml --cluster "$(pwd)/bin/immediate_submit.py $(echo {dependencies} | sed 's#(##g' | sed 's#)##g' | sed 's#"##g') sge" --immediate-submit -pr --notemp --latency-wait 600 $1
snakemake --use-singularity --jobs 1001 --cluster-config $CLUSTER_CONFIG --cluster "$(pwd)/bin/immediate_submit.py '{dependencies}' $CLUSTER" --immediate-submit -pr --notemp --latency-wait 600 $SM_ARGS
#snakemake --use-singularity --jobs 1001 --cluster-config data/cluster_config-sauron.yaml --immediate-submit -pr --notemp --latency-wait 600 $1
