#!/bin/bash
#SBATCH --job-name=r-odaf_DESeq2_Report
#SBATCH --account=def-someaccount
#SBATCH --cpus-per-task=6
#SBATCH --mem=12G
#SBATCH --time=1:00:00
#SBATCH --array=1-3

echo "Starting task $SLURM_ARRAY_TASK_ID"
# First create a file in the pwd that has a list of the dirs to be entered in the array
# Select the study ID directory based on the array index
# Set the STUDY_ID_DIR variable using the sed command
STUDY_ID_DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" array_dir_list.txt)
# Print the STUDY_ID_DIR to verify
echo "Processing directory: $STUDY_ID_DIR"
export STUDY_ID_DIR=${STUDY_ID_DIR}

# Load Apptainer module
module load apptainer/1.2.4

# Get the current working directory (note, this script must be executed from the appropriate root project directory, i.e., the root cloned directory if using cloned git repo)
CURRENT_DIR=$(pwd)

# Define the alignement/quantification option chosen in the nextflow pipeline
# This should match the option used in the nf-params.json file
# and should be one of "star_salmon", "star_rsem", "salmon", or "kallisto".
# The default is "star_salmon" but you can change it to match your choice
# in the nf-params.json file.
export ALIGN_QUANT_METHOD="star_salmon"

# Define the target script to run
TARGET_SCRIPT="R-ODAF_3_render_DESeq2_report.R"

# Create an empty directory to mask /root/.local and avoid conda environment contamination
# This is a workaround to avoid issues with conda environments in Apptainer
# Note: This directory will be created in the SLURM temporary directory
# and will be removed automatically after the job is done
mkdir -p ${SLURM_TMPDIR}/empty_local

# Run the apptainer command with the current directory bound
# Note: ~/projects/def-someaccount/group_writable/r-odaf_default.sif should point to wherever your apptainer image is stored
# chmod 750 is owner can read, write and execute, group can read and execute... first apply to folders, then also to files
# source /opt/miniconda3/etc/profile.d/conda.sh loads the conda profile
# conda activate DESeq2_report_env activates the appropriate environment containing all the necessary packages
apptainer exec --fakeroot --bind ${CURRENT_DIR}:/mnt --bind ${SLURM_TMPDIR}/empty_local:/root/.local ~/projects/def-someaccount/group_writable/r-odaf_default.sif \
bash -c "find /mnt/${STUDY_ID_DIR} -type d -maxdepth 10 -exec chmod 750 {} \; 2>/dev/null || true && \
find /mnt/${STUDY_ID_DIR} -type f -maxdepth 10 -exec chmod 750 {} \; 2>/dev/null || true && \
source /opt/miniconda3/etc/profile.d/conda.sh && \
conda activate DESeq2_report_env && \
Rscript /mnt/${STUDY_ID_DIR}/scripts/${TARGET_SCRIPT}"