#!/usr/bin/env python3

import os
import sys
import shutil

from snakemake.utils import read_job_properties

# last command-line argument is the job script
jobscript = sys.argv[-1]
# last but one argument is the submission system
subs = sys.argv[-2]


print(sys.argv, file=sys.stderr)
# all other command-line arguments are the dependencies
dependencies = ""
if subs == "slurm":
	dependencies = set(sys.argv[1:-2])
elif subs == "sge":
	dependencies = dependencies.join(sys.argv[1:-2])
else:
	print("Cannot get dependencies for submission system")
	sys.exit(1)
#print("Dependencies found")
#print(sys.argv[1:-2])
# parse the job script for the job properties that are encoded by snakemake within
# this also includes the information contained in the cluster-config file as job_properties["cluster"]
job_properties = read_job_properties(jobscript)
print("Job properties: ", job_properties, file=sys.stderr)
print("Rule wildcards: ", job_properties["wildcards"], file=sys.stderr)
print("Submission system: ", subs, file=sys.stderr)
print("Dependencies: ", dependencies, file=sys.stderr)
cmdline=[]
# create list with command line arguments
if subs == "slurm":
	cmdline = ["sbatch"]
	if "sample" in job_properties["wildcards"]:
		job_properties["cluster"]["J"] = job_properties["cluster"]["J"]+"-"+job_properties["wildcards"]["sample"]
		prefix = job_properties["wildcards"]["sample"] + "-" + job_properties["rule"] + "-slurm"
	else:
		prefix=job_properties["cluster"]["J"]
	job_properties["cluster"]["output"] = job_properties["cluster"]["output"].replace("slurm", prefix)
	job_properties["cluster"]["error"] = job_properties["cluster"]["error"].replace("slurm", prefix)
	# create string for slurm submit options for rule
	slurm_args = "--partition={partition} --time={time} --qos={qos} --ntasks={ntasks} --ntasks-per-node={ntasks-per-node} --hint={hint} --output={output} --error={error} -N {N} -J {J}".format(**job_properties["cluster"])
	cmdline.append(slurm_args)

	# now work on dependencies
	if dependencies:
		cmdline.append("--dependency")
		# only keep numbers (which are the jobids) in dependencies list. this is necessary because slurm returns more than the jobid. For other schedulers this could be different!
		dependencies = [x for x in dependencies if x.isdigit()]
		print(dependencies)
		cmdline.append("afterok:" + ",".join(dependencies))
		print(dependencies)
elif subs == "sge":
	name = job_properties["cluster"]["N"]
	job_properties["cluster"]["N"] = job_properties["cluster"]["N"]
	prefix = "comparative-" + job_properties["rule"] + "-sge"
	job_properties["cluster"]["output"] = job_properties["cluster"]["output"].replace("slurm", prefix).replace("%j",name)
	job_properties["cluster"]["error"] = job_properties["cluster"]["error"].replace("slurm", prefix).replace("%j",name)
	cmdline = ["qsub"]
	sge_args = "-cwd -V -q {queue} -l h_vmem={mem} -pe {queue} {ntasks} -o {output} -e {error} -N {N}".format(**job_properties["cluster"])	
	cmdline.append(sge_args)

	#now work on dependencies
	if dependencies:
		cmdline.append("-hold_jid")
		# only keep numbers (which are the jobids) in dependencies list. this is necessary because sge on sauron returns more than the jobid. For other schedulers this could be different!
		dependencies = [x for x in dependencies.split(" ") if x.isdigit()]
		cmdline.append(",".join(dependencies))
		print(dependencies, file=sys.stderr)
else:
	print("Immediate submit error: Unkown submission system!", file=sys.stderr)
	sys.exit(1)


cmdline.append(jobscript)


#now write final commandback to the system
print(" ".join(cmdline), file=sys.stderr)
os.system(" ".join(cmdline))


