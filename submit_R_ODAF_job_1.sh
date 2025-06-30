#!/bin/bash
# No need to populate SBATCH params because nextflow will take care of that

#Note: please install nextflow according to the documentation - https://docs.alliancecan.ca/wiki/Nextflow

# Load Apptainer modules
module load python/3.11
module load rust
module load postgresql
#Activate the nexflow virtual environment (Wherever you stored it during installation)
source ~/projects/def-someaccount/group_writable/nf-core-env/bin/activate
# Load the required modules
module load nextflow/24.04.4
module load apptainer/1

# Get the current working directory (note, this script must be executed from the appropriate root project directory, i.e., the root cloned directory if using cloned git repo)
CURRENT_DIR=$(pwd)
# Define the Study id
STUDY_ID="Study_id_N_Compound_n_day_exposure"
export STUDY_ID_DIR=${STUDY_ID}
# Check if the STUDY_ID_DIR is set
if [ -z "$STUDY_ID_DIR" ]; then
    echo "Error: STUDY_ID_DIR is not set. Please set it to the appropriate study ID directory."
    exit 1
fi
# Define the target script to run
TARGET_SCRIPT="R-ODAF_1_sequencing_DataPreProcess.sh"
# Check if the target script exists
if [ ! -f "${CURRENT_DIR}/${STUDY_ID_DIR}/scripts/${TARGET_SCRIPT}" ]; then
    echo "Error: Target script ${TARGET_SCRIPT} not found in the expected location."
    exit 1
fi
# Define the Nextflow parameters file
# Make sure to populate the nf-params.json file with appropriate paramaters.
# The nf-params.json file should be in the same directory as TARGET_SCRIPT
# and should contain the necessary parameters for the Nextflow pipeline.
# To create params use this interface: https://nf-co.re/launch/?pipeline=rnaseq&release=3.18.0
export PARAMS_JSON="nf-params.json"
# Check if the nf-params.json file exists
if [ ! -f "${CURRENT_DIR}/${STUDY_ID_DIR}/scripts/${PARAMS_JSON}" ]; then
    echo "Error: nf-params.json file not found in the expected location."
    exit 1
fi   

# Print the STUDY_ID_DIR to verify
echo "Current Directory: $CURRENT_DIR"
echo "Processing directory: $STUDY_ID_DIR"
 
# Run the Nextflow pipeline
bash "${CURRENT_DIR}/${STUDY_ID_DIR}/scripts/${TARGET_SCRIPT}"