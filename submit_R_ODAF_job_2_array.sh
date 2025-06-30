#!/bin/bash
#SBATCH --job-name=r-odaf_QC
#SBATCH --account=def-someaccount
#SBATCH --cpus-per-task=6
#SBATCH --mem=6G
#SBATCH --time=2:00:00
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

# Define the target script to run
TARGET_SCRIPT="R-ODAF_2_render_studywide_QC_report.R"

# Define the alignement/quantification option chosen in the nextflow pipeline
# This should match the option used in the nf-params.json file
# and should be one of "star_salmon", "star_rsem", "salmon", or "kallisto".
# The default is "star_salmon" but you can change it to match your choice
# in the nf-params.json file.
export ALIGN_QUANT_METHOD="star_salmon"
# Format the count matrix for input to the R script
# Note: This assumes that the count matrix is in the format produced by the nf-core/rnaseq pipeline
# and that the output files are named with the pattern "*.merged.gene_counts.tsv"
# Check if the count matrix files exist
if [ ! -d "${STUDY_ID_DIR}/output/${ALIGN_QUANT_METHOD}/" ]; then
    echo "Error: Directory ${STUDY_ID_DIR}/output/${ALIGN_QUANT_METHOD}/ does not exist."
    exit 1
fi
if [ -z "$(ls -A ${STUDY_ID_DIR}/output/${ALIGN_QUANT_METHOD}/*.merged.gene_counts.tsv 2>/dev/null)" ]; then
    echo "Error: No count matrix files found in ${STUDY_ID_DIR}/output/${ALIGN_QUANT_METHOD}/."
    exit 1
fi
# Check if the format_count_matrix.sh script exists
if [ ! -f "${CURRENT_DIR}/${STUDY_ID_DIR}/scripts/format_count_matrix.sh" ]; then
    echo "Error: format_count_matrix.sh script not found in the expected location."
    exit 1
fi
# Format the count matrix
cd "${STUDY_ID_DIR}/output/${ALIGN_QUANT_METHOD}/"
bash "${CURRENT_DIR}/${STUDY_ID_DIR}/scripts/format_count_matrix.sh" *.merged.gene_counts.tsv
cd "${CURRENT_DIR}"

# Calculate Uniquely mapped and Multi-mapped reads
if [[ "${ALIGN_QUANT_METHOD}" == "kallisto" ]]; then
    cd "${STUDY_ID_DIR}/output/${ALIGN_QUANT_METHOD}/"
    # Check if the run_info.json files exist
    if [ -z "$(ls -A ERR*/run_info.json 2>/dev/null)" ]; then
        echo "Error: No run_info.json files found in ${STUDY_ID_DIR}/output/${ALIGN_QUANT_METHOD}/."
        exit 1
    fi
    # Check if the calculate_unique_multimapped_reads_kallisto.sh script exists
    if [ ! -f "${CURRENT_DIR}/${STUDY_ID_DIR}/scripts/calculate_unique_multimapped_reads_kallisto.sh" ]; then
        echo "Error: calculate_unique_multimapped_reads_kallisto.sh script not found in the expected location."
        exit 1
    fi
    # Calculate the unique and multi-mapped reads
    bash "${CURRENT_DIR}/${STUDY_ID_DIR}/scripts/calculate_unique_multimapped_reads_kallisto.sh" > "${STUDY_ID_DIR}/output/${ALIGN_QUANT_METHOD}/unique_multimapped_reads.txt"
    cd "${CURRENT_DIR}"
elif [[ "${ALIGN_QUANT_METHOD}" == "salmon" ]]; then
    cd "${STUDY_ID_DIR}/output/${ALIGN_QUANT_METHOD}/"
    # Check if the ambig_info.tsv files exist
    if [ -z "$(ls -A ERR*/aux_info/ambig_info.tsv 2>/dev/null)" ]; then
        echo "Error: No ambig_info.tsv files found in ${STUDY_ID_DIR}/output/${ALIGN_QUANT_METHOD}/."
        exit 1
    fi
    # Check if the calculate_unique_multimapped_reads_salmon.sh script exists
    if [ ! -f "${CURRENT_DIR}/${STUDY_ID_DIR}/scripts/calculate_unique_multimapped_reads_salmon.sh" ]; then
        echo "Error: calculate_unique_multimapped_reads_salmon.sh script not found in the expected location."
        exit 1
    fi
    # Calculate the unique and multi-mapped reads
    bash "${CURRENT_DIR}/${STUDY_ID_DIR}/scripts/calculate_unique_multimapped_reads_salmon.sh" > "${STUDY_ID_DIR}/output/${ALIGN_QUANT_METHOD}/unique_multimapped_reads.txt"
    cd "${CURRENT_DIR}"
else
    echo "Mapping statistics located in multiqc directory."
fi


# Create an empty directory to mask /root/.local and avoid conda environment contamination
# This is a workaround to avoid issues with conda environments in Apptainer
# Note: This directory will be created in the SLURM temporary directory
# and will be removed automatically after the job is done
mkdir -p ${SLURM_TMPDIR}/empty_local

# Run the apptainer command with the current directory bound
# Note: ~/projects/def-someaccount/group_writable/r-odaf_default.sif should point to wherever your apptainer image is stored
# chmod 750 is owner can read, write and execute, group can read and execute... first apply to folders, then also to files
# source /opt/miniconda3/etc/profile.d/conda.sh loads the conda profile
# conda activate updated_my_r_pkgs activates the appropriate environment containing all the necessary packages
apptainer exec --fakeroot --bind ${CURRENT_DIR}:/mnt --bind ${SLURM_TMPDIR}/empty_local:/root/.local ~/projects/def-someaccount/group_writable/r-odaf_default.sif \
bash -c "find /mnt/${STUDY_ID_DIR} -type d -maxdepth 10 -exec chmod 750 {} \; 2>/dev/null || true && \
find /mnt/${STUDY_ID_DIR} -type f -maxdepth 10 -exec chmod 750 {} \; 2>/dev/null || true && \
source /opt/miniconda3/etc/profile.d/conda.sh && \
conda activate updated_my_r_pkgs && \
Rscript /mnt/${STUDY_ID_DIR}/scripts/${TARGET_SCRIPT}"

#Interactively using shell example:
#apptainer shell --fakeroot --bind ~/scratch/2024_tPOD_Curation_Project/GSE228670/:/mnt ~/projects/def-someaccount/group_writable/r-odaf_default.sif
#source /opt/miniconda3/etc/profile.d/conda.sh
#conda activate updated_my_r_pkgs