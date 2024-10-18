#!/bin/bash
#SBATCH --job-name=r-odaf_QC_X
#SBATCH --account=def-someaccount
#SBATCH --cpus-per-task=6
#SBATCH --mem=4G
#SBATCH --time=1:00:00

# Load Apptainer module
module load apptainer/1.2.4

# Get the current working directory (note, this script must be executed from the appropriate root project directory, i.e., the root cloned directory if using cloned git repo)
CURRENT_DIR=$(pwd)

# Define the target script to run
TARGET_SCRIPT="R-ODAF_2_render_studywide_QC_report.R"

# Define the study ID directory
export STUDY_ID_DIR="Study_dir"

# Run the apptainer command with the current directory bound
# Note: ~/projects/def-someaccount/group_writable/r-odaf_default.sif should point to wherever your apptainer image is stored
# chmod 750 is owner can read, write and execute, group can read and execute... first apply to folders, then also to files
# source /opt/miniconda3/etc/profile.d/conda.sh loads the conda profile
# conda activate updated_my_r_pkgs activates the appropriate environment containing all the necessary packages
apptainer exec --fakeroot --bind ${CURRENT_DIR}:/mnt ~/projects/def-someaccount/group_writable/r-odaf_default.sif \
bash -c "find /mnt -type d -exec chmod 750 {} \; && \
find /mnt -type f -exec chmod 750 {} \; && \
source /opt/miniconda3/etc/profile.d/conda.sh && \
conda activate updated_my_r_pkgs && \
Rscript /mnt/${STUDY_ID_DIR}/scripts/${TARGET_SCRIPT}"

#Interactively using shell example:
#apptainer shell --fakeroot --bind ~/scratch/2024_tPOD_Curation_Project/GSE228670/:/mnt ~/projects/def-someaccount/group_writable/r-odaf_default.sif
#source /opt/miniconda3/etc/profile.d/conda.sh
#conda activate updated_my_r_pkgs