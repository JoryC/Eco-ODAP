#!/bin/bash
####################################################################################
### Sequencing R-ODAF, Omics Data Analysis Framework for Regulatory application  ###
####################################################################################

# Optional alignment to genome or transcriptome for full-aligners (STAR) or genome-transcriptome-hybrid for quasi-mapping (Salmon).
# RNASeq - Standard quantification: Fastp -> SortMeRNA -> STAR -> RSEM -> RSeQC -> MultiQC
# [Not yet fully supported] TempOSeq - Standard quantification: Fastp -> SortMeRNA -> STAR -> RSEM -> RSeQC -> MultiQC
# 3PrimeRNASeq - Standard quantification: Fastp -> SortMeRNA -> STAR -> RSEM -> RSeQC -> MultiQC
# 3PrimeRNASeq - Special quantification: Fastp -> SortMeRNA -> Salmon -> RSeQC -> MultiQC
# RNASeq - Special quantification: Fastp -> SortMeRNA -> Salmon -> RSeQC -> MultiQC

# To add later: Auto-infer strandedness before trimming, and verify strandedness using infer_strandedness.py
# To add later: genomic dna contamination removal
# To add later: Automatically create the decoy and gentrome files for salmon

# To improve: clean up STAR commands

# Major changes: Create snakemake for the pipeline? Use Nextflow?

####################################################
#### Settings which need to be adapted by user #####
####################################################
source /opt/miniconda3/etc/profile.d/conda.sh
project=${STUDY_ID_DIR} #Make this the directory name corresponding to the project.
PROJECT_DIR="/mnt/${STUDY_ID_DIR}" #Make this the full directory path corresponding to the project.
# Specify the directory for the output
OUTPUT_DIR="${PROJECT_DIR}/output/"
# Specify location of input fastq files. ALL FILES IN THE FOLDER WILL BE PROCESSED
RAW_SAMPLE_DIR="${PROJECT_DIR}/raw/"
# Specify extention of input files (".fastq" or ".fastq.gz")
SUFFIX_INPUTFILES=${SUFFIX_INPUTFILES}
# Specify the sequencing type (RNASeq, TempOSeq, or 3PrimeRNASeq)
SEQTYPE=${SEQTYPE}
if [ "${SEQTYPE}" != "RNASeq" ] && [ "${SEQTYPE}" != "TempOSeq" ] && [ "${SEQTYPE}" != "3PrimeRNASeq" ]; then
    echo "Error: SEQTYPE must be RNASeq, TempOSeq, or 3PrimeRNASeq. Exiting."
    exit 1
fi
# Specify the sequencing mode used to obtain the data
SEQMODE=${SEQMODE} #specify "paired" or "single" end mode

# Specify if you want to use Quasi-mapping instead of STAR-RSEM. Specify "Yes" or "No".
QUASI_MAP=${QUASI_MAP}
if [ "${QUASI_MAP}" == "Yes" ]; then
    echo "Quasi-mapping will be used."
    # Specify the decoy file for Quasi-mapping
    DECOYS=${DECOYS} #Specify the name of the decoy text file
    GENTROME=${GENTROME_FILE_NAME} #Specify the name of the gentrome file
else
    echo "STAR-RSEM full-mapping will be used."
    # Specify if you want to map to genome or transcriptome. Specify "genome" or "transcriptome"
    MAPTO=${MAPTO} 
    if [ "${MAPTO}" != "genome" ] &&
       [ "${MAPTO}" != "transcriptome" ]; then
        echo "Error: MAPTO must be genome or transcriptome. Exiting."
        exit 1
    fi
fi

# Specify the strandedness of the library type
STRANDEDNESS=${STRANDEDNESS} #Specify "none", "forward", "reverse"

# Specify the suffixes for paired-end reads
# Check if SEQMODE is "paired"
if [ "$SEQMODE" == "paired" ]; then
  # Define parameters for paired-end mode
  PAIRED_END_SUFFIX_FORWARD=${PAIRED_END_SUFFIX_FORWARD} #Customize suffix
  PAIRED_END_SUFFIX_REVERSE=${PAIRED_END_SUFFIX_REVERSE} #Customize suffix
  # Print a message for debugging (optional)
  echo "Paired-end mode detected. Using suffixes: $PAIRED_END_SUFFIX_FORWARD, $PAIRED_END_SUFFIX_REVERSE"
else
  # Handle other modes or provide a default action
  echo "SEQMODE is 'single'. No paired-end parameters set."
  echo "Please determine the mean and standard deviation of the fragment length of your single-end reads before continuing."
  FLDMEAN=${FLDMEAN} #Mean fragment length
  FLDSD=${FLDSD} #Standard deviation of fragment length
fi

# Choose the main organism for genome alignment (e.g "Rat_6.0.97"). {NOTE: This ID is a label specific for this script and is made for the user to identify which genome version was used. It can contain any text}.
# hg38 | Rat_6.0.84 | S1500 | HumanWT
ORGANISM_GENOME_ID=${ORGANISM_GENOME_ID}
# PATH/Directory in which the genome files are located
# ${HOME}/shared/dbs/human/hg38/ | ${HOME}/shared/dbs/rat/ensembl/rnor6_0/v84/genome
# ${HOME}/shared/dbs/biospyder/R-Scripts/Human_S1500_Surrogate/TSQR_Scripts_Human_Surrogate_1.2/reference/humansurrogate1_2
# ${HOME}/shared/dbs/biospyder/R-Scripts/Human_Whole_Transcriptome/TSQR_Scripts_HumanWT_1.1/reference/humanwt1_1
if [ "${QUASI_MAP}" == "Yes" ]; then
    echo "Quasi-mapping will be used."
    # Specify the decoy file for Quasi-mapping
    DECOYS=${DECOYS} #Specify the name of the decoy text file
    GENTROME_FILE_NAME=${GENTROME_FILE_NAME} #Specify the name of the gentrome file
else
    echo "STAR-RSEM will be used."
    # Filename of genome/transcriptome fasta file (without path)
    # Homo_sapiens_assembly38.fasta | Rnor_6.0.fa
    # S1500: humansurrogate1_2.fa
    # Human WT: humanwt1_1.fa
    GENOME_FILE_NAME=${GENOME_FILE_NAME}
fi
GENOME_FILES_DIR=${GENOME_FILES_DIR} #Copy the genome files (.fa and .gtf) into the mounted directory ad put into '/mnt/genome'. Mount the ACC_ID directory, not the Study_id

# Filename of GTF file (without path)
# hg38.ensGene.gtf | Rattus_norvegicus.Rnor_6.0.84.andERCC.gtf
# S1500: humansurrogate1_2.gtf
# Human WT: humanwt1_1.gtf
GTF_FILE_NAME=${GTF_FILE_NAME}
# Filename of the BED file (without path)
BED_FILE_NAME="$(echo "${GTF_FILE_NAME}" | sed 's/.gtf.gz$//; s/.gtf$//').bed"
# Filename of the genePred file (without path)
PRED_FILE_NAME="$(echo "${GTF_FILE_NAME}" | sed 's/.gtf.gz$//; s/.gtf$//').genePred"
#Specify if your reads are below the length of 50bp on average (Separate indexes need to be made for very small reads - STAR). "Yes" or "No".
BP50ORLESS=${BP50ORLESS}
#Specify the expected length of your reads
EXPECTEDLENGTHOFREADS=${EXPECTEDLENGTHOFREADS}
# Whether the genome indexing has already been done. When "Yes" is specified, the indexing will be skipped. If "No" The index will be made
GENOME_INDEX_ALREADY_AVAILABLE=${GENOME_INDEX_ALREADY_AVAILABLE} #Specify "Yes" or "No"
RSEM_INDEX_ALREADY_AVAILABLE=${RSEM_INDEX_ALREADY_AVAILABLE}   #Specify "Yes" or "No"
#Specify whether you are working with a large (=human) genome. Specify "Yes" when working with human or "No"
LARGE_GENOME=${LARGE_GENOME} #Zebrafish genome index keeps failing due to insufficient memory
SKIP_TRIMMING=${SKIP_TRIMMING} #Should the trimming step be skipped, if you have already trimmed the reads? "Yes" or "No"
SKIP_SORTMERNA=${SKIP_SORTMERNA} #Should the SortMeRNA step be skipped? "Yes" or "No"


################################
# System parameters from SLURM #
################################

# Specify amount of CPUs to use for alignment step (use 20 or 30)
#CPU_FOR_ALIGNMENT=${SLURM_CPUS_PER_TASK} #Requesting 40 CPUs with slurm total
# Specify amount of CPUs to use (use 6 or higher)
#CPU_FOR_OTHER=${SLURM_CPUS_PER_TASK}

# Commenting out the above because I am exporting the variables in the submit_R_ODAF_job_1.sh script

### No other input required ###

#################################
# Running the sequencing R-ODAF #
#################################

declare BASEDIR=${OUTPUT_DIR}            #specify working directory
declare SEQMODE=${SEQMODE}                 #specify "paired" or "single" end mode
declare QUASI_MAP=${QUASI_MAP}             #specify "Yes" or "No" to use Quasi-mapping instead of STAR-RSEM
declare RAW_SAMPLE_DIR=${RAW_SAMPLE_DIR} #specify location of fastq files

if [[ "${SUFFIX_INPUTFILES}" =~ ^\. ]]; then
    echo "Suffix has a leading period, continuing."
else
    echo "Suffix did not have a leading period, adding one now..."
    SUFFIX_INPUTFILES=${SUFFIX_INPUTFILES/#/.}
    echo "${SUFFIX_INPUTFILES}"
fi

declare SUFFIX_IN=${SUFFIX_INPUTFILES} # specify extension of input files (".fastq", ".fastq.gz", etc.)

### PAIRED ONLY ###
# Check if SEQMODE is "paired"
if [ "$SEQMODE" == "paired" ]; then
  # Define parameters for paired-end mode
  declare PAIR1=${PAIRED_END_SUFFIX_FORWARD}
  declare PAIR2=${PAIRED_END_SUFFIX_REVERSE}
  # Print a message for debugging (optional)
  echo "Paired-end mode detected. Declaring suffixes: $PAIRED_END_SUFFIX_FORWARD, $PAIRED_END_SUFFIX_REVERSE"
else
  # Handle other modes or provide a default action
  echo "SEQMODE is 'single'. No paired-end parameters set."
fi

#Specify reference genome files
declare GENOMEDIR=${GENOME_FILES_DIR}             #location of genome files
declare GTF="${GENOMEDIR}/${GTF_FILE_NAME}"       #genome GTF file
declare BED="${GENOMEDIR}/${BED_FILE_NAME}"       #genome BED file
declare PRED="${GENOMEDIR}/${PRED_FILE_NAME}"     #genome genePred file
declare GenomeID=${ORGANISM_GENOME_ID}            #Specify the genome name (e.g. Species+GenomeVersion or Species+DownloadDate) to prevent overwriting other indexed genomes

if [ "${QUASI_MAP}" == "Yes" ]; then
    declare DECOY_FILE="${GENOMEDIR}/${DECOYS}" #Specify the decoy file for Quasi-mapping
    declare GENTROME="${GENOMEDIR}/${GENTROME_FILE_NAME}" #Specify the gentrome file for Quasi-mapping
else
    declare GENOME="${GENOMEDIR}/${GENOME_FILE_NAME}" #genome fasta file
fi

declare LargeGenome=${LARGE_GENOME} #Specify "Yes" or "No" to indicate if Human=Large Genome is used
#Specify "Yes" or "No" to indicate if GenomeIndexing has already been done for STAR (if no is specified, the index will be made)
GenomeIndexDone=${GENOME_INDEX_ALREADY_AVAILABLE}
#Specify "Yes" or "No" to indicate if GenomeIndexing has already been done for RSEM (if no is specified, the index will be made)
GenomeIndexRSEMDone=${RSEM_INDEX_ALREADY_AVAILABLE}
declare TinyReads=${BP50ORLESS}             #Specify if the reads are smaller than or equal to 50BP. "Yes" or "No"
declare ReadLength=${EXPECTEDLENGTHOFREADS} #Specify the read length of the RNAseq data, integer.
#Should the trimming step be skipped, if you have already trimmed the reads?
declare SKIP_TRIMMING=${SKIP_TRIMMING} #"Yes" skips the trimming process
declare SKIP_SORTMERNA=${SKIP_SORTMERNA} #"Yes" skips the SortMeRNA process

################################
# System parameters from SLURM #
################################

declare CPUs=${CPU_FOR_OTHER}           #Specify amount of CPUs to use
declare CPUs_align=${CPU_FOR_ALIGNMENT} #Specify amount of CPUs to use for alignment step
declare RAM_LIMIT=${MEMORY_IN_BYTES}    #Specify the amount of RAM to use for indexing the genome

#######################################################################
### Defining parameters for script to run (no user input necessary) ###
#######################################################################

declare SOURCEDIR=${RAW_SAMPLE_DIR}
declare TRIMM_DIR="${BASEDIR}/Trimmed_reads/"
declare OUTPUTDIR=${TRIMM_DIR}
declare SORTMERNA_DIR="${OUTPUTDIR}/sortmerna_workdir/"
declare QC_DIR_fastp="${OUTPUTDIR}/fastpQCoutput/"
declare QC_DIR_multiQC="${OUTPUTDIR}/MultiQC/"
declare align_DIR="${OUTPUTDIR}/STAR/"
declare Quant_DIR="${OUTPUTDIR}/RSEM/"
declare RSEM_GENOMEDIR="${GENOMEDIR}/RSEM/"
declare QC_DIR_RSeQC="${OUTPUTDIR}/RSeQC/"

if [ "${QUASI_MAP}" == "Yes" ]; then
    declare QUASI_MAP_DIR="${OUTPUTDIR}/Salmon/"
fi
### TEMPOSEQ ONLY ###
#declare TEMPOSEQR="${HOME}/shared/projects/${project}/scripts/pete.star.script_v3.1.R" #This line absolutely needs user input. Put the TempO-SeqR_v3.1.R script in the main dir. Or you can make a 'scripts' dir and point to the R script in that dir
#declare TEMPOSEQFIGS="${HOME}/shared/projects/${project}/scripts/generate_figures.R"   #This line too... point to the generate_figures.R script... also consider renaming that script to something like 'TempO-Seq_generate_figures.R' because it is only ever used for TempO-Seq data...

declare SUFFIX1=${SUFFIX_IN}
declare SUFFIX_out="_trimmed${SUFFIX_IN}"

echo "Intializing directories..."
mkdir -p "${TRIMM_DIR}"
mkdir -p "${QC_DIR_fastp}"
mkdir -p "${QC_DIR_multiQC}"
mkdir -p "${align_DIR}"
mkdir -p "${Quant_DIR}"
mkdir -p "${RSEM_GENOMEDIR}"
mkdir -p "${QC_DIR_RSeQC}"
mkdir -p "${SORTMERNA_DIR}"
if [ "${QUASI_MAP}" == "Yes" ]; then
    mkdir -p "${QUASI_MAP_DIR}"
fi
echo "Done."

echo "Starting logging..."
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
mydate="$(date +'%d.%m.%Y.%H-%M')"
exec 1>"${OUTPUT_DIR}"/log_"${mydate}".out 2>&1

echo $SHELL
echo "Processing SEQTYPE: ${SEQTYPE}"

#--------------------------------------------------------------------------------------------------------------------------------------------------------#

##################################
### Trimming raw reads : Fastp ###
##################################

# Check if trimming should be skipped
if [ "${SKIP_TRIMMING}" == "Yes" ]; then
    echo "Reads have already been trimmed. Skipping Trimming."
else
    echo "Activating required software."
    conda activate fastp

    #############################
    # Trimming single end reads #
    #############################
    if [ "${SEQMODE}" == "single" ]; then
        declare FILES1=(${SOURCEDIR}*${SUFFIX1}) # Array initialization

        # Iterate over each input file
        for FILENAME in "${FILES1[@]}"; do
            BASENAME="${FILENAME:${#SOURCEDIR}:-${#SUFFIX1}}"  # Extract the base name without directory and suffix
            OUTPUT_FILE="${TRIMM_DIR}${BASENAME}${SUFFIX_out}"
            JSON_REPORT="${QC_DIR_fastp}${BASENAME}_fastp.json"
            HTML_REPORT="${QC_DIR_fastp}${BASENAME}_fastp.html"

            echo -e "[TRIMMING] fastp: [${BASENAME}]"

            # Skip processing if the output file already exists
            if [ -e "${OUTPUT_FILE}" ] && [ "${SKIP_TRIMMING}" == "Yes" ]; then
                echo "File ${OUTPUT_FILE} exists, skipping trimming."
                continue
            fi

            # Run fastp based on the sequencing type
            case "${SEQTYPE}" in
                RNASeq)
                    fastp --in1 "${FILENAME}" \
                        --out1 "${OUTPUT_FILE}" \
                        --json "${JSON_REPORT}" \
                        --html "${HTML_REPORT}" \
                        --cut_front --cut_front_window_size 1 \
                        --cut_front_mean_quality 3 \
                        --cut_tail \
                        --cut_tail_window_size 1 \
                        --cut_tail_mean_quality 3 \
                        --cut_right --cut_right_window_size 4 \
                        --cut_right_mean_quality 15 \
                        --length_required 36 \
                        --thread "${CPUs}"
                    ;;
                TempOSeq)
                    fastp --in1 "${FILENAME}" \
                        --out1 "${OUTPUT_FILE}" \
                        --json "${JSON_REPORT}" \
                        --html "${HTML_REPORT}" \
                        --trim_tail1 1 \
                        --cut_front --cut_front_window_size 1 \
                        --cut_front_mean_quality 3 \
                        --cut_tail \
                        --cut_tail_window_size 1 \
                        --cut_tail_mean_quality 3 \
                        --cut_right --cut_right_window_size 4 \
                        --cut_right_mean_quality 15 \
                        --length_required 50 \
                        --thread "${CPUs}"
                    ;;
                3PrimeRNASeq)
                    fastp --in1 "${FILENAME}" \
                        --out1 "${OUTPUT_FILE}" \
                        --json "${JSON_REPORT}" \
                        --html "${HTML_REPORT}" \
                        --trim_poly_x \
                        --poly_x_min_len 3 \
                        --cut_tail \
                        --cut_tail_window_size 1 \
                        --cut_tail_mean_quality 3 \
                        --cut_right --cut_right_window_size 4 \
                        --cut_right_mean_quality 15 \
                        --length_required 12 \
                        --correction \
                        --thread "${CPUs}"
                    ;;
                *)
                    echo "Sequencing type ${SEQTYPE} not recognized. Quitting."
                    exit 1
                    ;;
            esac
        done
    fi

    #############################
    # Trimming paired end reads #
    #############################
     if [ "${SEQMODE}" == "paired" ]; then
        declare FILES1=(${SOURCEDIR}*${PAIR1}${SUFFIX1}) # Array initialization

        # Iterate over each paired file
        for FILENAME in "${FILES1[@]}"; do
            READ1="${FILENAME}"
            READ2="${FILENAME:0:-${#SUFFIX1}-${#PAIR1}}${PAIR2}${SUFFIX1}"
            BASENAME="${READ1:${#SOURCEDIR}:-${#PAIR1}-${#SUFFIX1}}"  # Base name for logging and outputs
            OUTPUT_READ1="${TRIMM_DIR}${READ1:${#SOURCEDIR}:-${#SUFFIX1}}${SUFFIX_out}"
            OUTPUT_READ2="${TRIMM_DIR}${READ2:${#SOURCEDIR}:-${#SUFFIX1}}${SUFFIX_out}"
            JSON_REPORT="${QC_DIR_fastp}${BASENAME}_PE_fastp.json"
            HTML_REPORT="${QC_DIR_fastp}${BASENAME}_PE_fastp.html"

            echo -e "[TRIMMING] fastp: [${BASENAME}]"

            # Skip processing if both output files already exist
            if [ -e "${OUTPUT_READ1}" ] && [ -e "${OUTPUT_READ2}" ] && [ "${SKIP_TRIMMING}" == "Yes" ]; then
                echo "Files ${OUTPUT_READ1} and ${OUTPUT_READ2} exist, skipping trimming."
                continue
            fi

            # Run fastp based on the sequencing type
            case "${SEQTYPE}" in
                RNASeq)
                    fastp \
                        --in1 "${READ1}" \
                        --in2 "${READ2}" \
                        --out1 "${OUTPUT_READ1}" \
                        --out2 "${OUTPUT_READ2}" \
                        --json "${JSON_REPORT}" \
                        --html "${HTML_REPORT}" \
                        --cut_front --cut_front_window_size 1 \
                        --cut_front_mean_quality 3 \
                        --cut_tail \
                        --cut_tail_window_size 1 \
                        --cut_tail_mean_quality 3 \
                        --cut_right \
                        --cut_right_window_size 4 \
                        --cut_right_mean_quality 15 \
                        --length_required 36 \
                        --thread "${CPUs}"
                    ;;
                TempOSeq)
                    fastp \
                        --in1 "${READ1}" \
                        --in2 "${READ2}" \
                        --out1 "${OUTPUT_READ1}" \
                        --out2 "${OUTPUT_READ2}" \
                        --json "${JSON_REPORT}" \
                        --html "${HTML_REPORT}" \
                        --trim_tail1 1 \
                        --trim_tail2 1 \
                        --cut_front --cut_front_window_size 1 \
                        --cut_front_mean_quality 3 \
                        --cut_tail \
                        --cut_tail_window_size 1 \
                        --cut_tail_mean_quality 3 \
                        --cut_right --cut_right_window_size 4 \
                        --cut_right_mean_quality 15 \
                        --length_required 50 \
                        --thread "${CPUs}"
                    ;;
                3PrimeRNASeq)
                # There paramaters were tweaked specifically for Qiagen's 3' UPXome RNA-seq kit (2024_12)
                    fastp \
                        --in1 "${READ1}" \
                        --in2 "${READ2}" \
                        --out1 "${OUTPUT_READ1}" \
                        --out2 "${OUTPUT_READ2}" \
                        --json "${JSON_REPORT}" \
                        --html "${HTML_REPORT}" \
                        --detect_adapter_for_pe \
                        --trim_poly_x \
                        --poly_x_min_len 3 \
                        --trim_front2 6 \
                        --trim_front1 3 \
                        --length_required 12 \
                        --correction \
                        --thread "${CPUs}"
                    ;;
                *)
                    echo "Sequencing type ${SEQTYPE} not recognized. Quitting."
                    exit 1
                    ;;  
            esac
        done
    fi

    # Save the environment for reproducibility
    conda env export >"${OUTPUT_DIR}conda_environment_fastp.${mydate}.yml"

    # Deactivate the Conda environment
    conda deactivate
fi

################################
#INFORMATION ON TRIMMING PROCESS
#--cut_front --cut_front_window_size 1 --cut_front_mean_quality 3 == Trimmomatic  "LEADING:3"
#--cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3 == Trimmomatic  "TRAILING:3"
#--cut_right --cut_right_window_size 4 --cut_right_mean_quality 15 == Trimmomatic  "SLIDINGWINDOW:4:15"

#### If additional trimming is needed (see multiQCreport):
# add to R1: --trim_front1 {amount_bases} and/or --trim_tail1 {amount_bases}
# add to R2: --trim_front2 {amount_bases} and/or --trim_tail2 {amount_bases}

###########################################
### Removal of ribosomal RNA: SortMeRNA ###
###########################################

if [[ "${SKIP_SORTMERNA}" == "Yes" ]]; then
    echo "Skipping SortMeRNA."
else
    echo "Running SortMeRNA."
    echo "Activating required software."
    conda activate sortmerna

    echo "Building SortMeRNA index."
    if [ -z "$(ls -A "${SORTMERNA_DIR}" 2>/dev/null)" ]; then
        sortmerna --ref "${GENOMEDIR}/rRNA_database/smr_v4.3_default_db.fasta" \
                --index 1 \
                --workdir "${SORTMERNA_DIR}" \
                --threads ${CPUs} \
                -m ${RAM_LIMIT}
    else
        echo "SortMeRNA index already exists. Skipping indexing."
    fi

    #############################
    # Cleaning single end reads #
    #############################
    if [ "${SEQMODE}" == "single" ]; then
        echo "Running SortMeRNA for single-end reads."
        declare FILES1=(${OUTPUTDIR}*${SUFFIX1}) #Array intialization - WARNING this will break if filenames have spaces!
            for FILENAME in "${FILES1[@]}"; do
                READ1=${FILENAME}
                BASENAME=$(basename "${READ1}" "_trimmed${SUFFIX1}")
                WORKDIR="${SORTMERNA_DIR}/"

            echo -e "[rRNA_REMOVAL] sortmerna: [${BASENAME}]"

            sortmerna --ref "${GENOMEDIR}/rRNA_database/smr_v4.3_default_db.fasta" \
                        --reads "${READ1}" \
                        --workdir "${WORKDIR}" \
                        --index 0 \
                        --threads ${CPUs} \
                        --fastx \
                        --zip-out \
                        --other "${WORKDIR}/out/${BASENAME}_non_rRNA" \
                        -e 1e-3

            rm -rf "${WORKDIR}/kvdb/"
            rm -rf "${WORKDIR}/readdb/"
            done
    else
    #############################
    # Cleaning paired end reads #
    #############################
        echo "Running SortMeRNA for paired-end reads."
        for FILE1 in "${OUTPUTDIR}"*"${PAIR1}"*"${SUFFIX1}"; do
            FILE2=${FILE1/${PAIR1}/${PAIR2}}  # Get corresponding PAIR2
            BASENAME=$(basename "${FILE1}" "${PAIR1}_trimmed${SUFFIX1}")
            WORKDIR="${SORTMERNA_DIR}"

            echo -e "[rRNA_REMOVAL] sortmerna: [${BASENAME}]"

            sortmerna --ref "${GENOMEDIR}/rRNA_database/smr_v4.3_default_db.fasta" \
                        --reads "${FILE1}" \
                        --reads "${FILE2}" \
                        --workdir "${WORKDIR}" \
                        --index 0 \
                        --threads ${CPUs} \
                        --fastx \
                        --paired_in \
                        --zip-out \
                        --aligned "${WORKDIR}/out/${BASENAME}_rRNA" \
                        --other "${WORKDIR}/out/${BASENAME}_non_rRNA" \
                        --out2 \
                        -e 1e-3

            rm -rf "${WORKDIR}/kvdb/"
            rm -rf "${WORKDIR}/readdb/"
        done
    fi

    conda env export >"${OUTPUT_DIR}conda_environment_sortmerna.${mydate}.yml"
    conda deactivate
fi

cd ${BASEDIR}

#######################################################################
### Alignment and quantification of reads (paired end & single end) ###
#######################################################################

if [ "${QUASI_MAP}" == "Yes" ]; then
    
    echo "Running Salmon for quasi-mapping."
    
    conda activate quasi_map

    #########################
    # Creating Salmon index #
    #########################
    
    if [ ! -d "${GENOMEDIR}/salmon_index" ]; then
        mkdir -p "${GENOMEDIR}/salmon_index"
        #Assuming you have already built and uploaded the gentrome and decoy files...    
        salmon index -t "${GENTROME}" -d "${DECOY_FILE}" -i "${GENOMEDIR}/salmon_index" -p ${CPUs}
    else
        echo "Salmon index already exists. Skipping indexing."
    fi

    ######################################################
    # Quantifying reads using salmon quant quasi-mapping #
    ######################################################

    #############
    # 3' RNASeq #
    #############    
    if [ "${SEQMODE}" == "single" ] && [ "${SEQTYPE}" == "3PrimeRNASeq" ]; then
        echo "Running Salmon for single-end 3' reads."
        # Iterate over each file pattern for single-end reads
        for FILE in $(find "${SORTMERNA_DIR}/out/" -type f -regextype posix-extended -regex ".*_non_rRNA\.(fq|fastq|fq.gz|fastq.gz)" -print0 | xargs -0); do
            #Extract the basename from files ending in either .fq, .fastq, .fq.gz, or .fastq.gz
            BASENAME=$(basename "${FILE}" | sed -E 's/\.(fq|fastq)(\.gz)?$//')
            OUTPUT="${QUASI_MAP_DIR}${BASENAME}_quant"
            #--fldMean ${FLDMEAN} #Set the expected mean fragment length of the sequencing library which can be determine by looking at the fragment size distribution of the library from a bioanalyzer, tapestation or similar tool. This can slightly improve quantification results, especially for TPMs.
            #--fldSD ${FLDSD} #Set the standard deviation of the fragment length distribution of the library. This can slightly improve quantification results, especially for TPMs.
            #-w 200 #Reads "mapping" to more than this many locations are discarded
            #--noLengthCorrection #Important for 3' RNASeq data
            #--posBias #Since reads are 3' biased
            salmon quant -i "${GENOMEDIR}/salmon_index" -l A \
                         -r "${FILE}" --validateMappings --useVBOpt --noLengthCorrection --softclip --softclipOverhangs \
                         -p ${CPUs} -o "${OUTPUT}" -g "${GTF}" --writeMappings="${OUTPUT}_mappings.sam" -w 200
        done
        # More things to do go here
    elif [ "${SEQMODE}" == "paired" ] && [ "${SEQTYPE}" == "3PrimeRNASeq" ]; then
        echo "Running Salmon for paired-end 3' reads."
        # Iterate over each file pattern for paired-end reads
        for FILE1 in $(find "${SORTMERNA_DIR}/out/" -type f -regextype posix-extended -regex ".*_non_rRNA_fwd\.(fq|fastq|fq.gz|fastq.gz)" -print0 | xargs -0); do
            FILE2="${FILE1/_fwd/_rev}"  # Get corresponding PAIR2
            if [ -f "${FILE2}" ]; then  # Ensure FILE2 exists
                #Extract the basename from files ending in either .fq, .fastq, .fq.gz, or .fastq.gz AND remove the '_fwd' suffix
                BASENAME=$(basename "${FILE1}" | sed -E 's/\.(fq|fastq)(\.gz)?$//; s/_fwd$//')
                OUTPUT="${QUASI_MAP_DIR}${BASENAME}_quant"
                salmon quant -i "${GENOMEDIR}/salmon_index" -l A \
                             -1 "${FILE1}" -2 "${FILE2}" --validateMappings --useVBOpt --noLengthCorrection --softclip --softclipOverhangs \
                            -p ${CPUs} -o "${OUTPUT}" -g "${GTF}" --writeMappings="${OUTPUT}_mappings.sam" -w 200
            fi
        done
        # More things to do go here
    ##########
    # RNASeq #
    ##########
    elif [ "${SEQMODE}" == "single" ] && [ "${SEQTYPE}" == "RNASeq" ]; then
        echo "Running Salmon for single-end reads."
        # Iterate over each file pattern for single-end reads
        for FILE in $(find "${SORTMERNA_DIR}/out/" -type f -regextype posix-extended -regex ".*_non_rRNA\.(fq|fastq|fq.gz|fastq.gz)" -print0 | xargs -0); do
            #Extract the basename from files ending in either .fq, .fastq, .fq.gz, or .fastq.gz
            BASENAME=$(basename "${FILE}" | sed -E 's/\.(fq|fastq)(\.gz)?$//')
            OUTPUT="${QUASI_MAP_DIR}${BASENAME}_quant"
            salmon quant -i "${GENOMEDIR}/salmon_index" -l A \
                         -r "${FILE}" --validateMappings --seqBias --gcBias --useVBOpt --softclip --softclipOverhangs \
                         -p ${CPUs} -o "${OUTPUT}" -g "${GTF}" --writeMappings="${OUTPUT}_mappings.sam" -w 200
        done
    elif [ "${SEQMODE}" == "paired" ] && [ "${SEQTYPE}" == "RNASeq" ]; then
        echo "Running Salmon for paired-end reads."
        # Iterate over each file pattern for paired-end reads
        for FILE1 in $(find "${SORTMERNA_DIR}/out/" -type f -regextype posix-extended -regex ".*_non_rRNA_fwd\.(fq|fastq|fq.gz|fastq.gz)" -print0 | xargs -0); do
            FILE2="${FILE1/_fwd/_rev}"  # Get corresponding PAIR2
            if [ -f "${FILE2}" ]; then  # Ensure FILE2 exists
                #Extract the basename from files ending in either .fq, .fastq, .fq.gz, or .fastq.gz AND remove the '_fwd' suffix
                BASENAME=$(basename "${FILE1}" | sed -E 's/\.(fq|fastq)(\.gz)?$//; s/_fwd$//')
                OUTPUT="${QUASI_MAP_DIR}${BASENAME}_quant"
                salmon quant -i "${GENOMEDIR}/salmon_index" -l A \
                             -1 "${FILE1}" -2 "${FILE2}" --validateMappings --seqBias --gcBias --posBias --useVBOpt --noLengthCorrection --softclip --softclipOverhangs \
                            -p ${CPUs} -o "${OUTPUT}" -g "${GTF}" --writeMappings="${OUTPUT}_mappings.sam" -w 200
            fi
        done
    else
        echo "Sequencing type ${SEQTYPE} and/or sequencing mode ${SEQMODE} not recognized. Quitting."
        exit 1
    fi

    conda env export >"${OUTPUT_DIR}conda_environment_salmon.${mydate}".yml
    conda deactivate
    conda activate samtools

    #################################################################
    ### Generating BAM files from quasi-mapping for downstream QC ###
    #################################################################

    cd ${BASEDIR}

    # Test write permissions
    if ! touch "${align_DIR}/test" 2>/dev/null; then
        echo "ERROR: Cannot write to output directory: ${align_DIR}"
    else
        rm "${align_DIR}/test"
        echo "Write permissions confirmed on: ${align_DIR}"
    fi
    for samfile in $(find "${QUASI_MAP_DIR}" -type f -name "*_quant_mappings.sam"); do
        echo "Converting SAM to BAM for ${samfile}"
        BASENAME=$(basename "${samfile}" "_quant_mappings.sam")
        echo "[SAMPLE:] ${BASENAME}"
        OUTPUT_BAM="${align_DIR}/${BASENAME}Aligned.sortedByCoord.out.bam"
        echo "Will create: ${OUTPUT_BAM}"

        # Creating tmp directory for this file
        TMP_DIR="${QUASI_MAP_DIR}/samtools_tmp_${BASENAME}"
        mkdir -p "${TMP_DIR}"

        # Process sam file
        echo "Step 1: Filtering SAM file and removing decoy sequences..."
        samtools view -@ ${CPUS} -h -F 0x900 "${samfile}" | \
        awk '!/XT:A:D/ && !/DS:D/' | \
        samtools sort -@ ${CPUS} -T "${TMP_DIR}" -o "${OUTPUT_BAM}" -
    done
    conda deactivate

    ###########################################################
    ### Summarize the results of the Salmon quantification ####
    ###########################################################

    # Extract the gene names for each transcript from the GTF annotation file
    bash ${PROJECT_DIR}/scripts/extract_tx2gene.sh "${GTF}" "${QUASI_MAP_DIR}/tx2gene.tsv"

    conda activate salmon_R
    Rscript ${PROJECT_DIR}/scripts/summarize_salmon_quant_results.R ${QUASI_MAP_DIR} ${QUASI_MAP_DIR}/tx2gene.tsv ${Quant_DIR} ${SEQTYPE}
    conda deactivate

else


    if [ "${SEQTYPE}" == "RNASeq" ] || [ "${SEQTYPE}" == "3PrimeRNASeq" ]; then
    echo "Activating required alignment software - STAR."
    conda activate star

        #####################################
        # Indexing reference genome for STAR #
        #####################################
        cd ${GENOMEDIR}
        if [ ${GenomeIndexDone} == "No" ]; then
            if [ ${TinyReads} == "Yes" ]; then
                # If reads are smaller than or equal to 50bp, adjust sjdbOverhang based on expected read length
                overhangValue=$((${ReadLength} - 1))
                echo "Indexing genome with adjusted sjdbOverhang for reads smaller than 50bp"
            else
                # Default sjdbOverhang value
                overhangValue=99
            fi

            if [ ${LargeGenome} == "Yes" ]; then
                echo "Indexing Large(=Human) genome with adjusted parameters"
                STAR \
                    --runMode genomeGenerate \
                    --genomeDir ${GENOMEDIR} \
                    --genomeFastaFiles ${GENOME} \
                    --sjdbGTFfile ${GTF} \
                    --sjdbOverhang ${overhangValue} \
                    --runThreadN ${CPUs} \
                    --genomeSAsparseD 1 \
                    --genomeChrBinNbits 15 \
                    --limitGenomeGenerateRAM=${RAM_LIMIT}
            else
                echo "Indexing Small(=non-human) genome with adjusted parameters"
                # You may need to adjust param --genomeSAindexNbases based on the genome size
                STAR \
                    --runMode genomeGenerate \
                    --genomeDir ${GENOMEDIR} \
                    --genomeFastaFiles ${GENOME} \
                    --sjdbGTFfile ${GTF} \
                    --sjdbOverhang ${overhangValue} \
                    --runThreadN ${CPUs}
            fi
        fi
        cd ${BASEDIR}

        #################################################
        # STAR Alignment Script for RNA-Seq Experiments #
        #################################################

        # Define common STAR parameters for each scenario
        COMMON_PARAMS="--alignIntronMax 1000000 --alignIntronMin 20 --alignMatesGapMax 1000000 \
                    --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
                    --outFilterMismatchNmax 999 --outFilterMultimapNmax 20"

        # 3' RNA-Seq specific parameters
        THREEPRIME_PARAMS="--twopassMode Basic --outFilterMismatchNoverLmax 0.1 --outFilterType BySJout"

        # Function to run STAR alignment for transcriptome
        run_star_transcriptome() {
            local read1="$1"
            local read2="$2"  # Empty for single-end
            local prefix="$3"
            local is_compressed="$4"
            local seq_type="$5"
            local map_to="$6"
    
            echo "Creating Aligned.toTranscriptome.out.bam file..."
    
            local cmd="STAR --runThreadN ${CPUs_align} --genomeDir ${GENOMEDIR} --readFilesIn \"${read1}\""
    
            # Add read2 for paired-end
            if [[ -n "$read2" ]]; then
                cmd+=" \"${read2}\""
            fi
    
            # Common parameters
            cmd+=" ${COMMON_PARAMS}"
    
            # Add compression parameter
            if [[ "$is_compressed" == "yes" ]]; then
                cmd+=" --readFilesCommand zcat"
            fi
    
            # Add sequence type specific parameters
            if [[ "$seq_type" == "3PrimeRNASeq" ]]; then
                cmd+=" ${THREEPRIME_PARAMS}"
            fi
    
            # Parameters specific to transcriptome output
            if [[ "$map_to" == "transcriptome" ]]; then
                cmd+=" --quantMode TranscriptomeSAM --outSAMtype BAM Unsorted \
                    --outSAMattributes NH HI AS NM MD --quantTranscriptomeBan Singleend"
            else
                cmd+=" --quantMode TranscriptomeSAM"
            fi
    
            # Add output prefix
            cmd+=" --outFileNamePrefix \"${prefix}\""
    
            # Execute the command
            eval $cmd
        }

        # Function to run STAR alignment for genome
        run_star_genome() {
            local read1="$1"
            local read2="$2"  # Empty for single-end
            local prefix="$3"
            local is_compressed="$4"
            local seq_type="$5"
    
            echo "Creating Aligned.sortedByCoord.out.bam file..."
    
            local cmd="STAR --runThreadN ${CPUs_align} --genomeDir ${GENOMEDIR} --readFilesIn \"${read1}\""
    
            # Add read2 for paired-end
            if [[ -n "$read2" ]]; then
                cmd+=" \"${read2}\""
            fi
    
            # Common parameters
            cmd+=" ${COMMON_PARAMS}"
    
            # Add compression parameter
            if [[ "$is_compressed" == "yes" ]]; then
                cmd+=" --readFilesCommand zcat"
            fi
    
            # Add sequence type specific parameters
            if [[ "$seq_type" == "3PrimeRNASeq" ]]; then
                cmd+=" ${THREEPRIME_PARAMS}"
            fi
    
            # Parameters specific to genome output
            cmd+=" --outSAMtype BAM SortedByCoordinate"
    
            # Add output prefix
            cmd+=" --outFileNamePrefix \"${prefix}\""
    
            # Execute the command
            eval $cmd
        }

        # Process reads based on sequence mode
        if [ ${SEQMODE} == "single" ]; then
            echo "Processing single-end reads"
            FILES1=($(find "${SORTMERNA_DIR}/out/" -type f -regextype posix-extended -regex ".*_non_rRNA\.(fq|fastq|fq.gz|fastq.gz)" -print0 | xargs -0))
    
            for FILENAME in "${FILES1[@]}"; do
                READ1=${FILENAME}
                BASENAME=$(basename "${READ1}" | sed -E 's/\.(fq|fastq)(\.gz)?$//')
                FILE_PREFIX="${align_DIR}${BASENAME}"
                BAM="${FILE_PREFIX}Aligned.toTranscriptome.out.bam"
        
                # Skip if output exists
                if [ -e "${BAM}" ]; then
                    echo "BAM file exists for ${READ1}, skipping alignment."
                    continue
                fi
        
                echo -e "[ALIGNING] STAR : ${BASENAME}"
        
                # Determine compression status
                EXT="${FILENAME##*.}"
                IS_COMPRESSED="no"
                [[ "${EXT}" == "gz" ]] && IS_COMPRESSED="yes"
                [[ "${IS_COMPRESSED}" == "yes" ]] && echo "Gzipped file detected, using zcat to read"
        
                # Run appropriate STAR command(s) based on mapping target
                if [[ "${MAPTO}" == "genome" ]]; then
                    # For genome mapping, run both transcriptome and genome alignments
                    run_star_transcriptome "${READ1}" "" "${FILE_PREFIX}" "${IS_COMPRESSED}" "${SEQTYPE}" "${MAPTO}"
                    run_star_genome "${READ1}" "" "${FILE_PREFIX}" "${IS_COMPRESSED}" "${SEQTYPE}"
                else
                    # For transcriptome mapping, run only transcriptome alignment
                    run_star_transcriptome "${READ1}" "" "${FILE_PREFIX}" "${IS_COMPRESSED}" "${SEQTYPE}" "${MAPTO}"
                fi
            done
    
        elif [ ${SEQMODE} == "paired" ]; then
            echo "Processing paired-end reads"
            FILES1=($(find "${SORTMERNA_DIR}/out/" -type f -regextype posix-extended -regex ".*_non_rRNA_fwd.*\.(fq|fastq|fq.gz|fastq.gz)" -print0 | xargs -0))
    
            for FILENAME in "${FILES1[@]}"; do
                # Get corresponding read2
                READ1="${FILENAME}"
                READ2="${READ1/_fwd/_rev}"
        
                # Extract basename
                BASENAME=$(basename "${READ1}" | sed -E 's/\.(fq|fastq)(\.gz)?$//; s/_fwd$//')
                FILE_PREFIX="${align_DIR}${BASENAME}"
                BAM1="${FILE_PREFIX}Aligned.toTranscriptome.out.bam"
                BAM2="${FILE_PREFIX}Aligned.sortedByCoord.out.bam"
        
                # Skip if outputs exist
                if [ -e "${BAM1}" ] && [ -e "${BAM2}" ]; then
                    echo "Both BAM files exist, skipping alignment."
                    continue
                fi
        
                echo -e "[ALIGNING] STAR : [${BASENAME}]"

                # Determine compression status
                EXT="${FILENAME##*.}"
                IS_COMPRESSED="no"
                [[ "${EXT}" == "gz" ]] && IS_COMPRESSED="yes"
                [[ "${IS_COMPRESSED}" == "yes" ]] && echo "Gzipped file detected, using zcat to read"

                # Run appropriate STAR command(s) based on mapping target
                if [[ "${MAPTO}" == "genome" ]]; then
                    # For genome mapping, run both transcriptome and genome alignments
                    run_star_transcriptome "${READ1}" "${READ2}" "${FILE_PREFIX}" "${IS_COMPRESSED}" "${SEQTYPE}" "${MAPTO}"
                    run_star_genome "${READ1}" "${READ2}" "${FILE_PREFIX}" "${IS_COMPRESSED}" "${SEQTYPE}"
                else
                    # For transcriptome mapping, run only transcriptome alignment
                    run_star_transcriptome "${READ1}" "${READ2}" "${FILE_PREFIX}" "${IS_COMPRESSED}" "${SEQTYPE}" "${MAPTO}"
                fi
            done
        else
            echo "Error: SEQMODE must be 'single' or 'paired'"
            exit 1
        fi

        conda env export >"${OUTPUT_DIR}conda_environment_star.${mydate}.yml"
        conda deactivate
        
        if [[ "${MAPTO}" == "genome" ]]; then

            #######################
            # QUANTIFICATION RSEM #
            #######################
            echo "Activating required software."
            conda activate rsem

            #Indexing reference genome for RSEM
            echo "Indexing genome"
            cd ${GENOMEDIR}
            if [ ${GenomeIndexRSEMDone} == "No" ]; then
                rsem-prepare-reference --gtf "${GTF}" "${GENOME}" "${RSEM_GENOMEDIR}${GenomeID}"
                GenomeIndexRSEMDone="Yes"
            fi

            cd ${BASEDIR}

            # Convert STRANDEDNESS to RSEM-compatible flag
            if [[ "$STRANDEDNESS" == "none" ]]; then
                STRAND_OPTION="--strandedness none"
            elif [[ "$STRANDEDNESS" == "forward" ]]; then
                STRAND_OPTION="--strandedness forward"
            elif [[ "$STRANDEDNESS" == "reverse" ]]; then
                STRAND_OPTION="--strandedness reverse"
            fi

            if [ "${SEQTYPE}" == "3PrimeRNASeq" ]; then
                echo "Quantifying 3' RNASeq reads."
                # Quantify reads single end
                if [[ "${SEQMODE}" == "single" ]]; then
                    echo "Quantifying single end"
                    FILES1=($(find "${SORTMERNA_DIR}/out/" -type f -regextype posix-extended -regex ".*_non_rRNA\.(fq|fastq|fq.gz|fastq.gz)" -print0 | xargs -0))
                    for FILENAME in "${FILES1[@]}"; do
                        READ1=${FILENAME}
                        #Extract the basename from files ending in either .fq, .fastq, .fq.gz, or .fastq.gz
                        BASENAME=$(basename "${READ1}" | sed -E 's/\.(fq|fastq)(\.gz)?$//')
                        BAM="${align_DIR}${BASENAME}Aligned.toTranscriptome.out.bam"
                        OUTPUT_PREFIX="${Quant_DIR}${BASENAME}"
                        echo -e "[QUANTIFYING] RSEM : [${BASENAME}]"

                            #--fragment-length-mean 200 \ # Define the mean fragment length from library QC electropheerogram off of the Bioanalyzer/Tapestation 
                            #--fragment-length-sd 20 \ # Define the standard deviation of the fragment length from library QC electropheerogram off of the Bioanalyzer/Tapestation

                        rsem-calculate-expression \
                            -p ${CPUs} \
                            --bam "${BAM}" \
                            --no-bam-output \
                            --estimate-rspd \
                            $STRAND_OPTION \
                            ${RSEM_GENOMEDIR}${GenomeID} ${OUTPUT_PREFIX}

                    done
                # Quantify reads paired end
                elif [[ "${SEQMODE}" == "paired" ]]; then
                    echo "Quantifying paired end"
                    FILES1=($(find "${SORTMERNA_DIR}/out/" -type f -regextype posix-extended -regex ".*_non_rRNA_fwd.*\.(fq|fastq|fq.gz|fastq.gz)" -print0 | xargs -0))
                    for FILENAME in "${FILES1[@]}"; do
                        #Extract the basename from files ending in either .fq, .fastq, .fq.gz, or .fastq.gz AND remove the '_fwd' suffix
                        BASENAME=$(basename "${FILENAME}" | sed -E 's/\.(fq|fastq)(\.gz)?$//; s/_fwd$//')
                        READ1="${FILENAME}"
                        READ2="${READ1/_fwd/_rev}"  # Get corresponding PAIR2
                        BAM1="${align_DIR}${BASENAME}Aligned.toTranscriptome.out.bam"
                        #BAM2="${align_DIR}${BASENAME}Aligned.sortedByCoord.out.bam"
                        OUTPUT_PREFIX="${Quant_DIR}${BASENAME}"

                        echo -e "[QUANTIFYING] RSEM : [${BASENAME}]"

                        rsem-calculate-expression \
                            -p ${CPUs} \
                            --paired-end \
                            --bam "${BAM1}" \
                            --no-bam-output \
                            --estimate-rspd \
                            $STRAND_OPTION \
                            ${RSEM_GENOMEDIR}${GenomeID} ${OUTPUT_PREFIX}

                    done
                fi
            elif [[ "${SEQTYPE}" == "RNASeq" ]]; then
                # Quantify reads single end
                if [[ "${SEQMODE}" == "single" ]]; then
                    echo "Quantifying single end"
                    FILES1=($(find "${SORTMERNA_DIR}/out/" -type f -regextype posix-extended -regex ".*_non_rRNA\.(fq|fastq|fq.gz|fastq.gz)" -print0 | xargs -0))
                    for FILENAME in "${FILES1[@]}"; do
                        READ1=${FILENAME}
                        #Extract the basename from files ending in either .fq, .fastq, .fq.gz, or .fastq.gz
                        BASENAME=$(basename "${READ1}" | sed -E 's/\.(fq|fastq)(\.gz)?$//')
                        BAM="${align_DIR}${BASENAME}Aligned.toTranscriptome.out.bam"
                        OUTPUT_PREFIX="${Quant_DIR}${BASENAME}"

                        echo -e "[QUANTIFYING] RSEM : [${BASENAME}]"

                            #--fragment-length-mean 200 \ # Define the mean fragment length from library QC electropheerogram off of the Bioanalyzer/Tapestation 
                            #--fragment-length-sd 20 \ # Define the standard deviation of the fragment length from library QC electropheerogram off of the Bioanalyzer/Tapestation

                        rsem-calculate-expression \
                            -p ${CPUs} \
                            --bam "${BAM}" \
                            --no-bam-output \
                            $STRAND_OPTION \
                            ${RSEM_GENOMEDIR}${GenomeID} ${OUTPUT_PREFIX}
                    done
                # Quantify reads paired end
                elif [[ "${SEQMODE}" == "paired" ]]; then
                    echo "Quantifying paired end"
                    FILES1=($(find "${SORTMERNA_DIR}/out/" -type f -regextype posix-extended -regex ".*_non_rRNA_fwd.*\.(fq|fastq|fq.gz|fastq.gz)" -print0 | xargs -0))
                    for FILENAME in "${FILES1[@]}"; do
                        READ1="${FILENAME}"
                        READ2="${READ1/_fwd/_rev}"  # Get corresponding PAIR2
                        #Extract the basename from files ending in either .fq, .fastq, .fq.gz, or .fastq.gz AND remove the '_fwd' suffix
                        BASENAME=$(basename "${READ1}" | sed -E 's/\.(fq|fastq)(\.gz)?$//; s/_fwd$//')
                        BAM1="${align_DIR}${BASENAME}Aligned.toTranscriptome.out.bam"
                        #BAM2="${align_DIR}${BASENAME}Aligned.sortedByCoord.out.bam"
                        OUTPUT_PREFIX="${Quant_DIR}${BASENAME}"

                        echo -e "[QUANTIFYING] RSEM : [${BASENAME}]"

                        rsem-calculate-expression \
                            -p ${CPUs} \
                            --paired-end \
                            --bam "${BAM1}" \
                            --no-bam-output \
                            $STRAND_OPTION \
                            ${RSEM_GENOMEDIR}${GenomeID} ${OUTPUT_PREFIX}
                    done
                fi
            fi

            # Output summary file from RSEM results
            echo "Outputting RSEM results to .data.tsv files"
            cd ${Quant_DIR} #Better solution to cd before executing commands
            declare FILELIST=$(find ${Quant_DIR} -name "*genes.results" -printf "%f\t")
            rsem-generate-data-matrix $FILELIST >${Quant_DIR}/genes.data.tsv
            sed -i 's/\.genes.results//g' ${Quant_DIR}/genes.data.tsv #Removes '.genes.results' from the header

            declare FILELIST=$(find ${Quant_DIR} -name "*isoforms.results" -printf "%f\t") # Potential solution: changing the -printf argument from %f\t to %p\t so that it prints out the full filepath... the RSEM command was not working because of the lack of the full path
            rsem-generate-data-matrix $FILELIST >${Quant_DIR}/isoforms.data.tsv
            sed -i 's/\.genes.results//g' ${Quant_DIR}/isoforms.data.tsv

            conda env export >"${OUTPUT_DIR}conda_environment_rsem.${mydate}.yml"
            conda deactivate

            cd ${BASEDIR}

        else
            #########################
            # MAPTO = TRANSCRIPTOME #
            #########################
            ######################################################
            # Quantifying reads using salmon quant quasi-mapping #
            ######################################################

            #Note: One one .bam file is output from STAR for both single and paired end reads
            conda activate quasi_map

            #############
            # 3' RNASeq #
            #############    
            if [[ "${SEQTYPE}" == "3PrimeRNASeq" ]]; then
                echo "Running Salmon for single-end 3' reads."
                # Iterate over each file pattern for single-end reads
                for FILE in $(find "${align_DIR}" -type f -name "*_non_rRNAAligned.toTranscriptome.out.bam"); do
                    #Extract the basename from files
                    BASENAME=$(basename "${FILE}" | sed -E 's/_non_rRNAAligned\.toTranscriptome\.out\.bam$//')
                    OUTPUT="${QUASI_MAP_DIR}${BASENAME}_quant"
                    #--fldMean ${FLDMEAN} #Set the expected mean fragment length of the sequencing library which can be determine by looking at the fragment size distribution of the library from a bioanalyzer, tapestation or similar tool. This can slightly improve quantification results, especially for TPMs.
                    #--fldSD ${FLDSD} #Set the standard deviation of the fragment length distribution of the library. This can slightly improve quantification results, especially for TPMs.
                    #-w 200 #Reads "mapping" to more than this many locations are discarded
                    #--noLengthCorrection #Important for 3' RNASeq data
                    salmon quant -t "${GENOME}" -l A \
                                -a "${FILE}" -useVBOpt --noLengthCorrection \
                                -p ${CPUs} -o "${OUTPUT}" -g "${GTF}" --sampleOut -w 200
                done
            ##########
            # RNASeq #
            ##########
            elif [[ "${SEQTYPE}" == "RNASeq" ]]; then
                echo "Running Salmon for single-end reads."
                # Iterate over each file pattern for single-end reads
                for FILE in $(find "${align_DIR}" -type f -name "*_non_rRNAAligned.toTranscriptome.out.bam"); do
                    #Extract the basename from files 
                    BASENAME=$(basename "${FILE}" | sed -E 's/_non_rRNAAligned\.toTranscriptome\.out\.bam$//')
                    OUTPUT="${QUASI_MAP_DIR}${BASENAME}_quant"
                    salmon quant -t "${GENOME}" -l A \
                                -a "${FILE}" --seqBias --gcBias --useVBOpt \
                                -p ${CPUs} -o "${OUTPUT}" -g "${GTF}" --sampleOut --w 200
                done
            else
                echo "Sequencing type ${SEQTYPE} and/or sequencing mode ${SEQMODE} not recognized. Quitting."
                exit 1
            fi

            conda env export >"${OUTPUT_DIR}conda_environment_salmon.${mydate}".yml
            conda deactivate
            conda activate samtools

            #################################################################
            ### Generating BAM files from quasi-mapping for downstream QC ###
            #################################################################
            # Test write permissions
            if ! touch "${align_DIR}/test" 2>/dev/null; then
                echo "ERROR: Cannot write to output directory: ${align_DIR}"
            else
                rm "${align_DIR}/test"
                echo "Write permissions confirmed on: ${align_DIR}"
            fi
            # Find all postSample.bam files in subdirectories
            find "${QUASI_MAP_DIR}" -type f -name "postSample.bam" | while read -r bamfile; do
                echo "Processing ${bamfile}"
                # Extract directory name and remove _quant suffix
                dir_name=$(basename "$(dirname "$bamfile")")
                BASENAME=${dir_name%_quant}
                echo "Sample name: ${BASENAME}"
                # Create temp directory/prefix for this process
                TEMP_PREFIX="${QUASI_MAP_DIR}/samtools_tmp_${BASENAME}"
                # Process the BAM file
                samtools view -@ ${CPUS} -h -F 0x900 "${bamfile}" | \
                awk '!/XT:A:D/ && !/DS:D/' | \
                samtools sort -@ ${CPUS} -T "${TEMP_PREFIX}" -o "${align_DIR}/${BASENAME}Aligned.sortedByCoord.out.bam" -
            done
            conda deactivate

            ###########################################################
            ### Summarize the results of the Salmon quantification ####
            ###########################################################
            
            # Extract the gene names for each transcript from the GTF annotation file
            bash ${PROJECT_DIR}/scripts/extract_tx2gene.sh "${GTF}" "${QUASI_MAP_DIR}/tx2gene.tsv"

            conda activate salmon_R
            Rscript ${PROJECT_DIR}/scripts/summarize_salmon_quant_results.R ${QUASI_MAP_DIR} ${QUASI_MAP_DIR}/tx2gene.tsv ${Quant_DIR} ${SEQTYPE}
            conda deactivate

        fi

    elif [[ "${SEQTYPE}" == "TempOSeq" ]]; then
        #Okay so now we are going way back... If you have tempoSeq data this pipeline needs to be mega-editted, especially the TempO-SeqR_v3.1.R script
        #If using TempO-Seq data, the pipeline is using STAR to index, and align reads in the .R script, and then it is using the QuasR R package to quantify the aligned reads. This cutom TempO-Seq R script is going to need its own conda envirnment! The environment needs to have R installed and  the QuasR R package... If you need to know what packages to have installed in this environment you are going to have to look at the R-ODAF health Canada documentation
        # Command Line Arguments:
        # 1: FASTA Reference File
        # 2: Directory of FASTQ files to align
        # 3: Number of CPUs to use
        cd ${BASEDIR}
        echo "Parameters passed to Rscript for TempO-Seq Alignments..."
        echo ${TEMPOSEQR}
        echo ${GENOME}
        echo ${TRIMM_DIR}
        echo ${CPU_FOR_ALIGNMENT}
        Rscript ${TEMPOSEQR} ${GENOME} ${TRIMM_DIR} ${CPU_FOR_ALIGNMENT}
        declare READCOUNTS="${BASEDIR}/count_table.csv"
        declare MAPPEDUNMAPPED="${BASEDIR}/mapped_unmapped.csv"
        echo ${READCOUNTS}
        echo ${MAPPEDUNMAPPED}
        Rscript ${TEMPOSEQFIGS} ${READCOUNTS} ${MAPPEDUNMAPPED}
    fi
fi

#########################################################################################
### Quality control raw reads: Fastp + sortmerna + RSEM + STAR + RSeQC MultiQC report ###
#########################################################################################

echo "Activating required software."

cd ${align_DIR}

if [ "${SEQTYPE}" == "TempOSeq" ]; then
    #conda activate mymultiqc
    # conda activate odaf-star2.7.1 # Uncomment if genome index version differs.
    echo "TempOSeq not yet fully supported..."
elif [ "${SEQTYPE}" == "RNASeq" ]; then
    
    echo "Converting GTF file to BED file..."
    conda activate ucsc-gtf2bed
    gtfToGenePred -ignoreGroupsWithoutExons "${GTF}" "${PRED}"
    genePredToBed "${PRED}" "${BED}"
    conda deactivate

    echo "propcessing RNASeq data with RSeQC..."

    # Process BAM files with RSeQC tools for RNASeq
    for BAM in "${align_DIR}"*Aligned.sortedByCoord.out.bam; do
        BAM_BASENAME=$(basename "${BAM}" Aligned.sortedByCoord.out.bam)
        OUTPUT_FILE="${QC_DIR_RSeQC}/${BAM_BASENAME}_read_distribution.txt"
        
        echo "Calculating transcript coverage using RSeQC for ${BAM_BASENAME}..."
        conda activate RSeQC_env
        read_distribution.py -r "${BED}" -i "${BAM}" > "${OUTPUT_FILE}"
        conda deactivate

        echo "Indexing BAM file ${BAM_BASENAME}..."
        conda activate samtools
        samtools index "${BAM}" -@ ${CPUs}
        conda deactivate
    done
    echo "Calculating gene body coverage using RSeQC..."
    
    conda activate RSeQC_env
    geneBody_coverage.py -r "${BED}" -i "${align_DIR}" -o "${QC_DIR_RSeQC}/gene_body_coverage_output"
    conda deactivate

elif [ "${SEQTYPE}" == "3PrimeRNASeq" ]; then
    
    echo "Converting GTF file to BED file..."
    conda activate ucsc-gtf2bed
    gtfToGenePred -ignoreGroupsWithoutExons "${GTF}" "${PRED}"
    genePredToBed "${PRED}" "${BED}"
    conda deactivate

    echo "processing 3PrimeRNASeq data with RSeQC..."
    # Process BAM files with RSeQC tools for 3' RNASeq
    for BAM in "${align_DIR}"*Aligned.sortedByCoord.out.bam; do
        BAM_BASENAME=$(basename "${BAM}" Aligned.sortedByCoord.out.bam)
        OUTPUT_FILE="${QC_DIR_RSeQC}/${BAM_BASENAME}_read_distribution.txt"
        
        echo "[FILE ${BAM_BASENAME}]: Calculating biotype coverage for 3' end of transcripts using RSeQC..."
        conda activate RSeQC_env
        read_distribution.py -r "${BED}" -i "${BAM}" > "${OUTPUT_FILE}"
        conda deactivate
        
        echo "Indexing BAM file ${BAM_BASENAME}..."
        conda activate samtools
        samtools index "${BAM}" -@ ${CPUs}
        conda deactivate
        #A large proportion of reads should map to the 3' untranslated region (3' UTR) exons and CDS of genes. Given that 3' RNASeq targets the 3' end of transcripts, you should expect a high percentage of reads in this region and exons.
    done
    echo "Calculating gene body coverage using RSeQC..."
    
    conda activate RSeQC_env
    geneBody_coverage.py -r "${BED}" -i "${align_DIR}" -o "${QC_DIR_RSeQC}/gene_body_coverage_output"
    conda deactivate

else
    echo "Sequencing type not recognized. Quitting."
    exit 1
fi

cd ${BASEDIR}

# Running multiQC on fastp-output
# multiqc ${QC_DIR_fastp} --filename MultiQC_Report.html --outdir ${QC_DIR_multiQC}
conda activate multiqc
multiqc "${BASEDIR}" --cl-config "extra_fn_clean_exts: { '_fastp.json' }" --filename MultiQC_Report.html --outdir ${QC_DIR_multiQC}
###################################################################################################
conda env export >"${OUTPUT_DIR}conda_environment_multiqc.${mydate}.yml"
conda deactivate

echo "Pre-processing of data complete."