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
STUDY_ID="Study_id_N_Compound_n_day_exposure"

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
# IF you defined "single" end mode, define the following two fragment length distribution variables...
#export FLDMEAN=100 #Integer number
#export FLDSD=20 #Integer number

#Uncomment out suffixes if 'paired' and define suffixes...
export PAIRED_END_SUFFIX_FORWARD="_1"
export PAIRED_END_SUFFIX_REVERSE="_2"

# Define organism genome ID
# NOTE: This ID is a label specific for this script and is made for the user to identify which genome version was used. It can contain any text
export ORGANISM_GENOME_ID=""

# Define the dir to find the genome files like the .fa and .gtf files
export GENOME_FILES_DIR="/mnt/genome"

# Define the genome/transcriptome reference file names
export GENOME_FILE_NAME="somefile.fa"

# If you would like to use a quasi-mapper (Salmon) instead of a full mapper (STAR), set this to "Yes" and define the DECOYS and GENTROME_FILE_NAME variables
export QUASI_MAP="No" #"Yes" or "No". Would you like to use a quasi-mapping approach instead?
#export DECOYS="decoys.txt"
#export GENTROME_FILE_NAME="FHM_gentrome.fa"
# Specify if you would like to the genome or transcriptome reference file for quantification
export MAPTO="genome" #Specify "genome" or "transcriptome"... "genome" option aligns reads to the STAR genome index and projects reads on to the transcriptome, then maps the aligned reads to a genomic RSEM index. "transcriptome" option aligns reads to the STAR genome index and projects reads on to the transcriptome, then maps reads to the Salmon decoy-aware transcriptome index
# Specify the GTF file name
export GTF_FILE_NAME="somefile.gtf" # This file is necessary for genome mapping and quasi-mapping (to output gene counts in addition to transcript counts), but not for transcriptome mapping. Make sure there are no empty "transcript_ids" in the "transcript_id" field of the GTF file.

# Specify if your reads are below the length of 50bp on average (Separate indexes need to be made for very small reads - STAR). "Yes" or "No".
export BP50ORLESS="No"
# Define the expected length of your reads
export EXPECTEDLENGTHOFREADS=100

# Define the strandedness of the library prep. Specify "none", "forward", or "reverse"
export STRANDEDNESS="reverse"

# Whether the genome indexing has already been done. When "Yes" is specified, the indexing will be skipped. If "No" The index will be made
export GENOME_INDEX_ALREADY_AVAILABLE="No" #Specify "Yes" or "No"
export RSEM_INDEX_ALREADY_AVAILABLE="No"   #Specify "Yes" or "No"
export SKIP_TRIMMING="No" #Should the trimming step be skipped, if you have already trimmed the reads? "Yes" or "No"
export SKIP_SORTMERNA="No" #Should the sortmerna step be skipped, if you have already removed rRNA from reads in a previous run? "Yes" or "No"

#Specify whether you are working with a large (=human) genome. Specify "Yes" when working with human or "No"
export LARGE_GENOME="Yes" #Zebrafish genome index keeps failing due to insufficient memory

# Set CPU allocation based on SLURM settings
export CPU_FOR_ALIGNMENT=${SLURM_CPUS_PER_TASK}
export CPU_FOR_OTHER=${SLURM_CPUS_PER_TASK}
export MEMORY_IN_BYTES=$((${SLURM_MEM_PER_NODE} * 1048576)) #Converting memory in MB to bytes

# Debugging messages:
echo "============================================================"
echo "DEBUGGING VARIABLES"
echo "============================================================"
echo "Current Directory:          ${CURRENT_DIR}"
echo "Selected Study ID Directory:${STUDY_ID_DIR}"
echo "Target Script:              ${TARGET_SCRIPT}"
echo "------------------------------------------------------------"
echo "Genome Files Directory:     ${GENOME_FILES_DIR}"
echo "Genome File Name:           ${GENOME_FILE_NAME}"
echo "GTF File Name:              ${GTF_FILE_NAME}"
echo "------------------------------------------------------------"
echo "Suffix of Input Files:      ${SUFFIX_INPUTFILES}"
echo "Sequencing Type:            ${SEQTYPE}"
echo "Sequencing Mode:            ${SEQMODE}"
echo "Paired-end Suffix Forward:  ${PAIRED_END_SUFFIX_FORWARD}"
echo "Paired-end Suffix Reverse:  ${PAIRED_END_SUFFIX_REVERSE}"
echo "------------------------------------------------------------"
echo "Organism Genome ID:         ${ORGANISM_GENOME_ID}"
echo "Mapping Reference:          ${MAPTO}"
echo "Using Quasi-mapper?:        ${QUASI_MAP}"
echo "Genome Indexing Available:  ${GENOME_INDEX_ALREADY_AVAILABLE}"
echo "RSEM Indexing Available:    ${RSEM_INDEX_ALREADY_AVAILABLE}"
echo "SortMeRNA Skipped:          ${SKIP_SORTMERNA}"
echo "Trimming Skipped:           ${SKIP_TRIMMING}"
echo "------------------------------------------------------------"
echo "Library Type Strandedness:  ${STRANDEDNESS}"
echo "Reads Below 50bp:           ${BP50ORLESS}"
echo "Expected Length of Reads:   ${EXPECTEDLENGTHOFREADS}"
echo "Large Genome:               ${LARGE_GENOME}"
echo "------------------------------------------------------------"
echo "Memory in Bytes:            ${MEMORY_IN_BYTES}"
echo "CPU for Alignment:          ${CPU_FOR_ALIGNMENT}"
echo "CPU for Other:              ${CPU_FOR_OTHER}"
echo "============================================================"

# Create an empty directory to mask /root/.local
mkdir -p ${SLURM_TMPDIR}/empty_local

# Run the apptainer command with the current directory bound
# Note: ~/projects/def-someaccount/group_writable/r-odaf_default.sif should point to wherever your apptainer image is stored
# chmod 750 is owner can read, write and execute, group can read and execute... first apply to folders, then also to files
apptainer exec --fakeroot --bind ${CURRENT_DIR}:/mnt --bind ${SLURM_TMPDIR}/empty_local:/root/.local ~/projects/def-someaccount/group_writable/r-odaf_default.sif \
bash -c "find /mnt -type d -exec chmod 750 {} \; && \
find /mnt -type f -exec chmod 750 {} \; && \
/mnt/${STUDY_ID_DIR}/scripts/${TARGET_SCRIPT}"