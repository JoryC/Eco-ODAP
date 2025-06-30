#!/bin/bash
#SBATCH --job-name=r-odaf_MultiQC_Report
#SBATCH --account=def-someaccount
#SBATCH --cpus-per-task=6
#SBATCH --mem=6G
#SBATCH --time=0:30:00
#SBATCH --array=1-3

module load apptainer/1.2.4

CURRENT_DIR=$(pwd)

echo "Starting task $SLURM_ARRAY_TASK_ID"
# First create a file in the pwd that has a list of the dirs to be entered in the array
# Select the study ID directory based on the array index
# Set the STUDY_ID_DIR variable using the sed command
STUDY_ID_DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" array_dir_list.txt)
# Print the STUDY_ID_DIR to verify
echo "Processing directory: $STUDY_ID_DIR"
export STUDY_ID_DIR=${STUDY_ID_DIR}
#STUDY_ID_DIR=Study_id_N_Compound_n_day_exposure

# Create reports directory if it doesn't exist
mkdir -p "${CURRENT_DIR}/${STUDY_ID_DIR}/reports/"
echo "Created or confirmed reports directory at: ${CURRENT_DIR}/${STUDY_ID_DIR}/reports/"

#CONFIGURATION_FILE=custom_multiqc_config.yaml
#If using a custom configuration file, set the correct argumnent for the multiqc command:
#multiqc -c /mnt/${STUDY_ID_DIR}/scripts/${CONFIGURATION_FILE} /mnt/${STUDY_ID_DIR}/output/ -f -o /mnt/${STUDY_ID_DIR}/reports/

mkdir -p ${SLURM_TMPDIR}/empty_local

apptainer exec --fakeroot --bind ${CURRENT_DIR}:/mnt --bind ${SLURM_TMPDIR}/empty_local:/root/.local ~/projects/def-someaccount/group_writable/r-odaf_default.sif \
bash -c "find /mnt/${STUDY_ID_DIR} -type d -maxdepth 10 -exec chmod 750 {} \; 2>/dev/null || true && \
find /mnt/${STUDY_ID_DIR} -type f -maxdepth 10 -exec chmod 750 {} \; 2>/dev/null || true && \
source /opt/miniconda3/etc/profile.d/conda.sh && \
conda activate multiqc && \
multiqc /mnt/${STUDY_ID_DIR}/output/ -f -o /mnt/${STUDY_ID_DIR}/reports/"