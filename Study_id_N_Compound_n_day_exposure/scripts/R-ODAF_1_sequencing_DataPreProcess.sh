#!/bin/bash
####################################################################################
### Sequencing R-ODAF, Omics Data Analysis Framework for Regulatory application  ###
####################################################################################

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
# Specify the sequencing type (RNASeq, TempOSeq)
SEQTYPE=${SEQTYPE}
# Specify the sequencing mode used to obtain the data
SEQMODE=${SEQMODE} #specify "paired" or "single" end mode

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
fi

# Choose the main organism for genome alignment (e.g "Rat_6.0.97"). {NOTE: This ID is a label specific for this script and is made for the user to identify which genome version was used. It can contain any text}.
# hg38 | Rat_6.0.84 | S1500 | HumanWT
ORGANISM_GENOME_ID=${ORGANISM_GENOME_ID}
# PATH/Directory in which the genome files are located
# ${HOME}/shared/dbs/human/hg38/ | ${HOME}/shared/dbs/rat/ensembl/rnor6_0/v84/genome
# ${HOME}/shared/dbs/biospyder/R-Scripts/Human_S1500_Surrogate/TSQR_Scripts_Human_Surrogate_1.2/reference/humansurrogate1_2
# ${HOME}/shared/dbs/biospyder/R-Scripts/Human_Whole_Transcriptome/TSQR_Scripts_HumanWT_1.1/reference/humanwt1_1
GENOME_FILES_DIR=${GENOME_FILES_DIR} #Copy the genome files (.fa and .gtf) into the mounted directory ad put into '/mnt/genome'. Mount the ACC_ID directory, not the Study_id
# Filename of genome fasta file (without path)
# Homo_sapiens_assembly38.fasta | Rnor_6.0.fa
# S1500: humansurrogate1_2.fa
# Human WT: humanwt1_1.fa
GENOME_FILE_NAME=${GENOME_FILE_NAME}
# Filename of GTF file (without path)
# hg38.ensGene.gtf | Rattus_norvegicus.Rnor_6.0.84.andERCC.gtf
# S1500: humansurrogate1_2.gtf
# Human WT: humanwt1_1.gtf
GTF_FILE_NAME=${GTF_FILE_NAME}
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
declare RAW_SAMPLE_DIR=${RAW_SAMPLE_DIR} #specify location of fastq files
if [[ "${SUFFIX_INPUTFILES}" =~ ^\. ]]; then
    echo "Suffix has a leading period, continuing."
else
    echo "Suffix did not have a leading period, adding one now..."
    SUFFIX_INPUTFILES=${SUFFIX_INPUTFILES/#/.}
    echo ${SUFFIX_INPUTFILES}
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
declare GENOME="${GENOMEDIR}/${GENOME_FILE_NAME}" #genome fasta file
declare GTF="${GENOMEDIR}/${GTF_FILE_NAME}"       #genome GTF file
declare GenomeID=${ORGANISM_GENOME_ID}            #Specify the genome name (e.g. Species+GenomeVersion or Species+DownloadDate) to prevent overwriting other indexed genomes

declare LargeGenome=${LARGE_GENOME} #Specify "Yes" or "No" to indicate if Human=Large Genome is used
#Specify "Yes" or "No" to indicate if GenomeIndexing has already been done for STAR (if no is specified, the index will be made)
GenomeIndexDone=${GENOME_INDEX_ALREADY_AVAILABLE}
#Specify "Yes" or "No" to indicate if GenomeIndexing has already been done for RSEM (if no is specified, the index will be made)
GenomeIndexRSEMDone=${RSEM_INDEX_ALREADY_AVAILABLE}
declare TinyReads=${BP50ORLESS}             #Specify if the reads are smaller than or equal to 50BP. "Yes" or "No"
declare ReadLength=${EXPECTEDLENGTHOFREADS} #Specify the read length of the RNAseq data, integer.
#Should the trimming step be skipped, if you have already trimmed the reads?
declare SKIP_TRIMMING=${SKIP_TRIMMING} #"Yes" skips the trimming process

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
declare QC_DIR_fastp="${OUTPUTDIR}/fastpQCoutput/"
declare QC_DIR_multiQC="${OUTPUTDIR}/MultiQC/"
declare align_DIR="${OUTPUTDIR}/STAR/"
declare Quant_DIR="${OUTPUTDIR}/RSEM/"
declare RSEM_GENOMEDIR="${GENOMEDIR}/RSEM/"
### TEMPOSEQ ONLY ###
#declare TEMPOSEQR="${HOME}/shared/projects/${project}/scripts/pete.star.script_v3.1.R" #This line absolutely needs user input. Put the TempO-SeqR_v3.1.R script in the main dir. Or you can make a 'scripts' dir and point to the R script in that dir
#declare TEMPOSEQFIGS="${HOME}/shared/projects/${project}/scripts/generate_figures.R"   #This line too... point to the generate_figures.R script... also consider renaming that script to something like 'TempO-Seq_generate_figures.R' because it is only ever used for TempO-Seq data...

declare SUFFIX1=${SUFFIX_IN}
declare SUFFIX_out="_trimmed${SUFFIX_IN}"

echo "Intializing directories..."
mkdir -p ${TRIMM_DIR}
mkdir -p ${QC_DIR_fastp}
mkdir -p ${QC_DIR_multiQC}
mkdir -p ${align_DIR}
mkdir -p ${Quant_DIR}
mkdir -p ${RSEM_GENOMEDIR}
echo "Done."

echo "Starting logging..."
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
mydate="$(date +'%d.%m.%Y.%H-%M')"
exec 1>${OUTPUT_DIR}/log_${mydate}.out 2>&1

echo $SHELL

##################################
### Trimming raw reads : Fastp ###
##################################

# Check if trimming should be skipped
if [ "${SKIP_TRIMMING}" == "Yes" ]; then
    echo "Reads have already been trimmed. Skipping Trimming."
else
    echo "Activating required software."
    if [ "${SEQTYPE}" == "TempOSeq" ]; then
        conda activate fastp
        #conda activate odaf-star2.7.1 # If genome index is a different version...
    elif [ "${SEQTYPE}" == "RNASeq" ]; then
        conda activate fastp
    else
        echo "Sequencing type not recognized. There could be a typo in the SEQTYPE, please double-check. Quitting."
        exit 1
    fi

    #############################
    # Trimming single end reads #
    #############################
    if [ "${SEQMODE}" == "single" ]; then
        declare FILES1="${SOURCEDIR}*${SUFFIX1}"
        for FILENAME in ${FILES1[@]}; do
            echo -e "[TRIMMING] fastp: [${FILENAME:${#SOURCEDIR}:-${#SUFFIX1}}]"
            #Single end
            #Prevent overwrite:
            if [ -e ${TRIMM_DIR}${FILENAME:${#SOURCEDIR}:-${#SUFFIX1}}${SUFFIX_out} ]; then
                echo "File exists, continuing"
            else
                if [ "${SEQTYPE}" == "RNASeq" ]; then
                    fastp --in1 ${FILENAME} \
                        --out1 ${TRIMM_DIR}${FILENAME:${#SOURCEDIR}:-${#SUFFIX1}}${SUFFIX_out} \
                        --json "${QC_DIR_fastp}${FILENAME:${#SOURCEDIR}:-${#SUFFIX1}}_fastp.json" \
                        --html ${QC_DIR_fastp}"${FILENAME:${#SOURCEDIR}:-${#SUFFIX1}}fastp.html" \
                        --cut_front --cut_front_window_size 1 \
                        --cut_front_mean_quality 3 \
                        --cut_tail \
                        --cut_tail_window_size 1 \
                        --cut_tail_mean_quality 3 \
                        --cut_right --cut_right_window_size 4 \
                        --cut_right_mean_quality 15 \
                        --length_required 36 \
                        --thread ${CPUs}
                elif [ "${SEQTYPE}" == "TempOSeq" ]; then
                    fastp --in1 ${FILENAME} \
                        --out1 ${TRIMM_DIR}${FILENAME:${#SOURCEDIR}:-${#SUFFIX1}}${SUFFIX_out} \
                        --json "${QC_DIR_fastp}${FILENAME:${#SOURCEDIR}:-${#SUFFIX1}}_fastp.json" \
                        --html ${QC_DIR_fastp}"${FILENAME:${#SOURCEDIR}:-${#SUFFIX1}}fastp.html" \
                        --trim_tail1 1 \
                        --cut_front --cut_front_window_size 1 \
                        --cut_front_mean_quality 3 \
                        --cut_tail \
                        --cut_tail_window_size 1 \
                        --cut_tail_mean_quality 3 \
                        --cut_right --cut_right_window_size 4 \
                        --cut_right_mean_quality 15 \
                        --length_required 50 \
                        --thread ${CPUs}
                else
                    echo "Sequencing type not recognized. Quitting"
                    exit 1
                fi
            fi
        done
    fi

    #############################
    # Trimming paired end reads #
    #############################
    if [ "${SEQMODE}" == "paired" ]; then
        declare FILES1="${SOURCEDIR}*${PAIR1}${SUFFIX1}"
        for FILENAME in ${FILES1[@]}; do
            READ1=${FILENAME}
            READ2=${FILENAME:0:-${#SUFFIX1}-${#PAIR1}}${PAIR2}${SUFFIX1}
            echo -e "[TRIMMING] fastp: [${READ1:${#SOURCEDIR}:-${#PAIR1}-${#SUFFIX1}}]"
            #Prevent overwrite:
            if [ -e ${TRIMM_DIR}${READ1:${#SOURCEDIR}:-${#SUFFIX1}}${SUFFIX_out} ] &&
                [ -e ${TRIMM_DIR}${READ2:${#SOURCEDIR}:-${#SUFFIX1}}${SUFFIX_out} ]; then
                echo "Files exist, continuing"
            else
                if [ "${SEQTYPE}" == "RNASeq" ]; then
                    fastp \
                        --in1 ${READ1} \
                        --in2 ${READ2} \
                        --out1 ${TRIMM_DIR}${READ1:${#SOURCEDIR}:-${#SUFFIX1}}${SUFFIX_out} \
                        --out2 ${TRIMM_DIR}${READ2:${#SOURCEDIR}:-${#SUFFIX1}}${SUFFIX_out} \
                        --json "${QC_DIR_fastp}${READ1:${#SOURCEDIR}:-${#SUFFIX1}-${#PAIR1}}PE_fastp.json" \
                        --html "${QC_DIR_fastp}${READ1:${#SOURCEDIR}:-${#SUFFIX1}-${#PAIR1}}PE_fastp.html" \
                        --cut_front --cut_front_window_size 1 \
                        --cut_front_mean_quality 3 \
                        --cut_tail \
                        --cut_tail_window_size 1 \
                        --cut_tail_mean_quality 3 \
                        --cut_right \
                        --cut_right_window_size 4 \
                        --cut_right_mean_quality 15 \
                        --length_required 36 \
                        --thread ${CPUs}
                else
                    echo "Sequencing type not recognized. Quitting"
                    exit 1
                fi
            fi
        done
    fi

    conda env export >${OUTPUT_DIR}conda_environment_fastp.${mydate}.yml

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

####################################################
### Alignment of reads (paired end & single end) ###
####################################################

echo "Activating required software."
if [ ${SEQTYPE} == "TempOSeq" ]; then
    conda activate star
    #conda activate odaf-star2.7.1 # If genome index is a different version...
elif [ ${SEQTYPE} == "RNASeq" ]; then
    conda activate star
else
    echo "Sequencing type not recognized. Quitting."
    break
fi

if [ ${SEQTYPE} == "RNASeq" ]; then

    #####################################
    #Indexing reference genome for STAR #
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
            STAR \
                --runMode genomeGenerate \
                --genomeDir ${GENOMEDIR} \
                --genomeFastaFiles ${GENOME} \
                --sjdbGTFfile ${GTF} \
                --sjdbOverhang ${overhangValue} \
                --runThreadN ${CPUs}
        fi
    fi
    cd ${SOURCEDIR}
    #############################
    # Aligning reads single end #
    #############################
    if [ ${SEQMODE} == "single" ]; then
        declare FILES1=${OUTPUTDIR}*${PAIR1}*${SUFFIX1}
        for FILENAME in ${FILES1[@]}; do
            READ1=${FILENAME:0:-${#SUFFIX1}}
            #Prevent overwrite:
            if [ -e ${align_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1} + 8)}Log.final.out ]; then
                echo "File exists, continuing"
            else
                echo -e "[ALIGNING] STAR : [${READ1:${#OUTPUTDIR}:-${#PAIR1}}]"
                if [[ "${SUFFIX1}" == *".gz" ]]; then
                    echo "gzipped file detected, using zcat to read"
                    STAR \
                        --runThreadN ${CPUs_align} \
                        --genomeDir ${GENOMEDIR} \
                        --readFilesIn "${READ1}${SUFFIX1}" \
                        --quantMode TranscriptomeSAM \
                        --readFilesCommand zcat \
                        --outFileNamePrefix ${align_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1} + 8)}
                else
                    echo "uncompressed FASTQ format detected"
                    STAR \
                        --runThreadN ${CPUs_align} \
                        --genomeDir ${GENOMEDIR} \
                        --readFilesIn "${READ1}${SUFFIX1}" \
                        --quantMode TranscriptomeSAM \
                        --outFileNamePrefix ${align_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1} + 8)}
                fi
            fi
        done
    fi
    #############################
    # Aligning reads paired end #
    #############################
    if [ ${SEQMODE} == "paired" ]; then
        declare FILES1=${OUTPUTDIR}*${PAIR1}*${SUFFIX1}
        for FILENAME in ${FILES1[@]}; do
            READ1=${FILENAME:0:-${#SUFFIX1}}
            READ2=${FILENAME:0:-${#SUFFIX1}-(${#PAIR1} + 8)}${PAIR2}"_trimmed"
            #Prevent overwrite:
            if [ -e ${align_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1} + 8)}Log.final.out ]; then
                echo "File exists, continuing"
            else
                echo -e "[ALIGNING] STAR : [${READ1:${#OUTPUTDIR}:-(${#PAIR1} + 8)}]"
                if [[ "${SUFFIX1}" == *".gz" ]]; then
                    echo "gzipped file detected, using zcat to read"
                    STAR \
                        --runThreadN ${CPUs_align} \
                        --genomeDir ${GENOMEDIR} \
                        --readFilesIn "${READ1}${SUFFIX1}" "${READ2}${SUFFIX1}" \
                        --quantMode TranscriptomeSAM \
                        --readFilesCommand zcat \
                        --outFileNamePrefix ${align_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1} + 8)}
                else
                    echo "uncompressed FASTQ format detected"
                    STAR \
                        --runThreadN ${CPUs_align} \
                        --genomeDir ${GENOMEDIR} \
                        --readFilesIn "${READ1}${SUFFIX1}" "${READ2}${SUFFIX1}" \
                        --quantMode TranscriptomeSAM \
                        --outFileNamePrefix ${align_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1} + 8)}
                fi
            fi
        done
    fi

    conda env export >${OUTPUT_DIR}conda_environment_star.${mydate}.yml

    conda deactivate

    #######################
    # QUANTIFICATION RSEM #
    #######################
    #This if statement is not necessary
    echo "Activating required software."
    if [ ${SEQTYPE} == "TempOSeq" ]; then
        conda activate rsem
        #conda activate odaf-star2.7.1 # If genome index is a different version...
    elif [ ${SEQTYPE} == "RNASeq" ]; then
        conda activate rsem
    else
        echo "Sequencing type not recognized. Quitting."
        break
    fi

    #Indexing reference genome for RSEM
    echo "Indexing genome"
    cd ${GENOMEDIR}
    if [ ${GenomeIndexRSEMDone} == "No" ]; then
        rsem-prepare-reference --gtf ${GTF} ${GENOME} ${RSEM_GENOMEDIR}${GenomeID}
        GenomeIndexRSEMDone="Yes"
    fi

    echo "Quantifying single end"
    cd ${SOURCEDIR}
    # Quantify reads single end
    if [ ${SEQMODE} == "single" ]; then
        declare FILES1=${OUTPUTDIR}*${PAIR1}*${SUFFIX1}
        for FILENAME in ${FILES1[@]}; do
            READ1=${FILENAME:0:-${#SUFFIX1}}
            echo -e "[QUANTIFYING] RSEM : [${READ1:${#OUTPUTDIR}:-(${#PAIR1} + 8)}]"

            rsem-calculate-expression \
                -p ${CPUs} \
                --bam "${align_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1} + 8)}Aligned.toTranscriptome.out.bam" \
                --no-bam-output ${RSEM_GENOMEDIR}${GenomeID} ${Quant_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1} + 8)}
        done
    fi

    echo "Quantifying paired end"
    # Quantify reads paired end
    if [ ${SEQMODE} == "paired" ]; then
        declare FILES1=${OUTPUTDIR}*${PAIR1}*${SUFFIX1}
        for FILENAME in ${FILES1[@]}; do
            READ1=${FILENAME:0:-${#SUFFIX1}}
            READ2=${FILENAME:0:-${#SUFFIX1}-(${#PAIR1} + 8)}${PAIR2}"_trimmed"
            echo -e "[QUANTIFYING] RSEM : [${READ1:${#OUTPUTDIR}:-(${#PAIR1} + 8)}]"

            rsem-calculate-expression \
                -p ${CPUs} \
                --paired-end \
                --bam "${align_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1} + 8)}Aligned.toTranscriptome.out.bam" \
                --no-bam-output ${RSEM_GENOMEDIR}${GenomeID} ${Quant_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1} + 8)}
        done
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

else
    if [ ${SEQTYPE} == "TempOSeq" ]; then
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

conda env export >${OUTPUT_DIR}conda_environment_rsem.${mydate}.yml

conda deactivate

#####################################################################
### Quality control raw reads: Fastp + RSEM + STAR MultiQC report ###
#####################################################################

echo "Activating required software."
if [ ${SEQTYPE} == "TempOSeq" ]; then
    conda activate mymultiqc
    #conda activate odaf-star2.7.1 # If genome index is a different version...
elif [ ${SEQTYPE} == "RNASeq" ]; then
    conda activate mymultiqc
else
    echo "Sequencing type not recognized. Quitting."
    break
fi

# Running multiQC on fastp-output
# multiqc ${QC_DIR_fastp} --filename MultiQC_Report.html --outdir ${QC_DIR_multiQC}
multiqc --cl_config "extra_fn_clean_exts: { '_fastp.json' }" ${BASEDIR} --filename MultiQC_Report.html --outdir ${QC_DIR_multiQC}

###################################################################################################

conda env export >${OUTPUT_DIR}conda_environment_multiqc.${mydate}.yml

conda deactivate
echo "Pre-processing of data complete."