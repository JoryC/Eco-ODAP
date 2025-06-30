#!/bin/bash
#SBATCH --job-name=kraken_build
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

#For instructions on how to download and build custom kraken2 databases, see: https://github.com/DerrickWood/kraken2/wiki/Manual#custom-databases

cd ~/scratch/kraken_dbs/transcriptomes_Species_of_interest_and_Standard

kraken2-build --build --threads ${SLURM_CPUS_PER_TASK} --db .

echo "Job finished at $(date)"
