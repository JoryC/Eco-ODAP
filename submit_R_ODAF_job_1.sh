#!/bin/bash
#SBATCH --job-name=r-odaf_preprocess_and_quant_X
#SBATCH --account=def-someaccount
#SBATCH --cpus-per-task=30
#SBATCH --mem=40G
#SBATCH --time=12:00:00

# Load Apptainer module
module load apptainer/1.2.4

# Get the current working directory (note, this script must be executed from the appropriate root project directory, i.e., the root cloned directory if using cloned git repo)
CURRENT_DIR=$(pwd)

# Define the target script to run
TARGET_SCRIPT="R-ODAF_1_sequencing_DataPreProcess.sh"

# Define the Study id
STUDY_ID="Study_dir"

# Print the STUDY_ID_DIR to verify
echo "Processing directory: $STUDY_ID_DIR"

# Define params
# Define the study ID directory
export STUDY_ID_DIR=${STUDY_ID}
# Specify extention of input files (".fastq" or ".fastq.gz")
export SUFFIX_INPUTFILES='.fastq.gz'
# Define Sequencing type (RNASeq or TempOSeq)
export SEQTYPE="RNASeq"
# Define Sequencing mode ("paired" or "single" end mode)
export SEQMODE="paired"
#Uncomment out suffixes if 'paired' and define suffixes...
export PAIRED_END_SUFFIX_FORWARD="_1"
export PAIRED_END_SUFFIX_REVERSE="_2"
# Define organism genome ID
# NOTE: This ID is a label specific for this script and is made for the user to identify which genome version was used. It can contain any text
export ORGANISM_GENOME_ID=""
# Define the dir to find the genome files like the .fa and .gtf files
export GENOME_FILES_DIR="/mnt/genome"
# Define the genome file names
export GENOME_FILE_NAME="somefile.fa"
export GTF_FILE_NAME="somefile.gtf"
# Specify if your reads are below the length of 50bp on average (Separate indexes need to be made for very small reads - STAR). "Yes" or "No".
export BP50ORLESS="No"
# Define the expected length of your reads
export EXPECTEDLENGTHOFREADS=100
# Whether the genome indexing has already been done. When "Yes" is specified, the indexing will be skipped. If "No" The index will be made
export GENOME_INDEX_ALREADY_AVAILABLE="No" #Specify "Yes" or "No"
export RSEM_INDEX_ALREADY_AVAILABLE="No"   #Specify "Yes" or "No"
#Specify whether you are working with a large (=human) genome. Specify "Yes" when working with human or "No"
export LARGE_GENOME="Yes" #Zebrafish genome index keeps failing due to insufficient memory
export SKIP_TRIMMING="No" #Should the trimming step be skipped, if you have already trimmed the reads? "Yes" or "No"
# Set CPU allocation based on SLURM settings
export CPU_FOR_ALIGNMENT=${SLURM_CPUS_PER_TASK}
export CPU_FOR_OTHER=${SLURM_CPUS_PER_TASK}
export MEMORY_IN_BYTES=$((${SLURM_MEM_PER_NODE} * 1048576)) #Converting memory in MB to bytes

# Debugging messages:
echo "Current Directory: ${CURRENT_DIR}"
echo "Selected Study ID Directory: ${STUDY_ID_DIR}"
echo "Target Script: ${TARGET_SCRIPT}"
echo "Genome Files Directory: ${GENOME_FILES_DIR}"
echo "Genome File Name: ${GENOME_FILE_NAME}"
echo "GTF File Name: ${GTF_FILE_NAME}"
echo "Memory in Bytes: ${MEMORY_IN_BYTES}"
echo "CPU for Alignment: ${CPU_FOR_ALIGNMENT}"
echo "CPU for Other: ${CPU_FOR_OTHER}"

# Run the apptainer command with the current directory bound
# Note: ~/projects/def-someaccount/group_writable/r-odaf_default.sif should point to wherever your apptainer image is stored
# chmod 750 is owner can read, write and execute, group can read and execute... first apply to folders, then also to files
apptainer exec --fakeroot --bind ${CURRENT_DIR}:/mnt ~/projects/def-someaccount/group_writable/r-odaf_default.sif \
bash -c "find /mnt -type d -exec chmod 750 {} \; && \
find /mnt -type f -exec chmod 750 {} \; && \
/mnt/${STUDY_ID_DIR}/scripts/${TARGET_SCRIPT}"