#!/bin/bash
#SBATCH --job-name=bracken_build
#SBATCH --account=def-someaccount
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=36:00:00

echo "Job started at $(date)"
echo "Running on $(hostname)"
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "Memory allocated: $SLURM_MEM_PER_NODE"

module load kraken2/2.1.3
module load bracken/3.0

#For insturction on how to build a bracken database see: https://ccb.jhu.edu/software/bracken/index.shtml?t=manual and https://github.com/jenniferlu717/Bracken

cd ~/scratch/kraken_dbs/transcriptomes_Species_of_interest_and_Standard

bracken-build -d . -t ${SLURM_CPUS_PER_TASK} -k 35 -l 65

echo "Job finished at $(date)"
