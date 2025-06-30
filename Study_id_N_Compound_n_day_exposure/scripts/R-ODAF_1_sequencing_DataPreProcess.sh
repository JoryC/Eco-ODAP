#!/bin/bash

# write your nextflow run command here
nextflow run nf-core-rnaseq_3.18.0/3_18_0/ \
    -profile singularity,narval \
    -params-file "${CURRENT_DIR}/${STUDY_ID_DIR}/scripts/${PARAMS_JSON}" \
    --outdir "${CURRENT_DIR}/${STUDY_ID_DIR}/output" \

# E.g., resume a previous run
# nextflow run -resume SESSION_ID nf-core-rnaseq_3.18.0/3_18_0/ \
#    -profile singularity,narval \
#    -params-file "${CURRENT_DIR}/${STUDY_ID_DIR}/scripts/${PARAMS_JSON}" \
#    --outdir "${CURRENT_DIR}/${STUDY_ID_DIR}/output" \

# Provide a custom configuration file, e.g.:
# nextflow run -resume SESSION_ID nf-core-rnaseq_3.18.0/3_18_0/ \
#    -c "${CURRENT_DIR}/${STUDY_ID_DIR}/scripts/bracken_custom.config" \
#    -profile singularity,narval \
#    -params-file "${CURRENT_DIR}/${STUDY_ID_DIR}/scripts/${PARAMS_JSON}" \
#    --outdir "${CURRENT_DIR}/${STUDY_ID_DIR}/output" \