#!/bin/bash

#SBATCH --time=16:00:00

#SBATCH --mem=4G

#SBATCH -e snakeSlurm.err.%J

#SBATCH -o snakeSlurm.out.%J

#SBATCH --partition=scu-cpu

## USAGE: runSnakemake.sh can be run either as an independent script or as an sbatch script
## sbatch runSnakemake.sh [snakeFile]
## bash runSnakemake.sh [snakeFile]

snakeFile=$1

#source cut-n-run-pipeline/condaInit.sh
conda activate snakemake

snakemake --version

if [[ ! -d ./Err_Out/ ]]; then
	mkdir ./Err_Out/
fi

if [[ ! -f ${snakeFile} ]];
then
	echo "snakefile does not exist!"
	exit 1
fi

echo "Using Snakefile: ${snakeFile}"

snakemake -s ${snakeFile} \
	--cores 1 \
	--rerun-incomplete \
	--cluster-config Config/slurmConfig.json \
	--use-conda \
	--conda-frontend conda \
	--latency-wait 200 \
	--printshellcmds \
	--cluster \
	"sbatch --partition={cluster.partition} -J {rule} -o Err_Out/slurm-%j.out -e Err_Out/slurm-%j.err -N1 -n {cluster.threads} --time {cluster.time} --mem={cluster.mem}" \
	--jobs 500 

	#--profile cut-n-run-pipeline/Config/slurmConfig.json \

#snakemake --rerun-incomplete --cluster-config cut-n-run-pipeline/Config/slurmConfig.json --printshellcmds -s ${snakeFile} --rerun-incomplete --latency-wait 200 --cluster  \
#	"sbatch --partition={cluster.partition} -J {rule} -o ./Err_Out/slurm-%j.out -e ./Err_Out/slurm-%j.err -N1 -n {cluster.threads} --time {cluster.time} --mem={cluster.mem}" --jobs 500


#snakemake --rerun-incomplete --cluster-config chip_pipeline/Config/slurmConfig.json --printshellcmds -s ${snakeFile} --rerun-incomplete --latency-wait 200 --cluster  \
#	"sbatch --partition=panda -J {rule} -o ./Err_Out/slurm-%j.out -e ./Err_Out/slurm-%j.err -N1 -n {cluster.threads} --time {cluster.time} --mem={cluster.mem}" --jobs 500
#pipepath=cut-n-run-pipeline
#snakemake --rerun-incomplete --printshellcmds -s ${snakeFile} --latency-wait 60 --cluster "sbatch -J {rule} -o ./Err_Out/slurm-%j.out -e ./Err_Out/slurm-%j.err -N1 --time 2:00:00 --mem=16G" --jobs 500
