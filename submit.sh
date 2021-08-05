#!/bin/bash

set -e

usage() { 
	echo "Welcome to the pipeline submission script. A script helps to submit jobs to SLURM and SGE clusters with snakemake and singularity"
	echo
	echo "Usage: $0 [-v] [-c <cluster_config_file>] [-s <snakemke_args>]" 
	echo
	echo "Options:"
	echo "  -t <submission_system> Specify available submission system. Options: sge, slurm, serial (no submission system). Default: Automatic detection."
	echo "	-c <cluster_config_file> Path to cluster config file in YAML format (mandatory). "
	echo "	-s \"<snakemake_args>\" Additional arguments passed on to the snakemake command (optional). snakemake is run with --immediate-submit -pr --notemp --latency-wait 600 --use-singularity --jobs 1001 by default." 
	echo "	-i \"<singularity_args>\" Additional arguments passed on to singularity (optional). Singularity is run with -B /tmp:/usertmp by default."
	1>&2; exit 1; }
	
version() {
	echo "$0 v0.1"
	exit 0
}

CLUSTER=""
while getopts ":v:t:c:s:i:" option;
	do
		case "${option}"
		in
			v) version;;
			t) CLUSTER=${OPTARG};;
			c) CLUSTER_CONFIG=${OPTARG};;
			s) SM_ARGS=${OPTARG};;
			i) SI_ARGS=${OPTARG};;
			*) echo "Illegal option --$OPTARG\n" >&2; usage;;
			?) echo "Illegal option --$OPTARG\n" >&2 usage;;
		esac
	done
if [ $OPTIND -eq 1 ]; then usage; fi

# Determine submission system:
if [[ $CLUSTER == "sge" ]]; then
	echo "SGE (Sun Grid Engine) submission system specified. Will use qsub to submit jobs."
elif [[ $CLUSTER == "slurm" ]]; then
	echo "SLURM submission system specified. Will use sbatch to submit jobs."
elif [[ $CLUSTER == "serial" ]]; then
  echo "Serial execution without job submission specified."
else
	echo "No or unknown submission system specified, will try to detect the system automatically."
	CLUSTER=""
	command -v qsub >/dev/null 2>&1 && { echo >&2 "SGE detected, will use qsub to submit jobs."; CLUSTER="sge"; }
	command -v sbatch >/dev/null 2>&1 && { echo >&2 "SLURM detected, will use sbatch to submit jobs."; CLUSTER="slurm"; }
  if [[ $CLUSTER == "" ]]; then
    echo "Submission system could not be detected. You may be able to run the pipeline without job submission."
    exit 1
  fi
fi

echo "Additional arguments passed on to singularity: $SI_ARGS"
echo "Additional arguments passed on to snakemake: $SM_ARGS"

if [ $CLUSTER = "slurm" ]; then
	export CONDA_PKGS_DIRS="$(pwd)/.conda_pkg_tmp"
	mkdir -p .conda_pkg_tmp
	snakemake --use-conda --use-singularity --singularity-args "-B $(pwd)/external/deeploc-1.0/bin/:/external -B $(pwd)/external/deeploc-1.0/DeepLoc:/usr/lib/python3/dist-packages/DeepLoc -B /tmp:/usrtmp -B $(pwd):/data -B $(pwd)/bin:/usr/local/external $SI_ARGS" --jobs 1001 --cluster-config $CLUSTER_CONFIG --cluster '$(pwd)/bin/immediate-submit/immediate_submit.py {dependencies} slurm' --immediate-submit -pr --notemp --latency-wait 60 $SM_ARGS
	unset CONDA_PKGS_DIRS
elif [ $CLUSTER = "sge" ]; then
	snakemake --use-conda --use-singularity --singularity-args "-B $(pwd)/external/deeploc-1.0/bin/:/external -B $(pwd)/external/deeploc-1.0/DeepLoc:/usr/lib/python3/dist-packages/DeepLoc -B $(pwd):/data -B $(pwd)/bin:/usr/local/external -B /tmp:/usertmp $SI_ARGS" --jobs 1001 --cluster-config $CLUSTER_CONFIG --cluster "$(pwd)/bin/immediate-submit/immediate_submit.py '{dependencies}' sge" --immediate-submit -pr --notemp --latency-wait 60 $SM_ARGS
elif [ $CLUSTER = "serial" ]; then
	snakemake --use-conda --use-singularity --singularity-args "-B $(pwd)/external/deeploc-1.0/bin/:/external -B $(pwd)/external/deeploc-1.0/DeepLoc:/usr/lib/python3/dist-packages/DeepLoc -B $(pwd):/data -B $(pwd)/bin:/usr/local/external -B /tmp:/usertmp $SI_ARGS" -pr --notemp --latency-wait 60 $SM_ARGS
else
	echo "Submission system not recognized"
	exit 1
fi
