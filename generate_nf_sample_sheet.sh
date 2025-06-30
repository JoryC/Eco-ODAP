#!/bin/bash

# Default values
FASTQ_DIR="."
OUTPUT_FILE="samples.csv"
PATTERN="*fastq.gz"
HELP=false

# Function to display usage information
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Generate a sample sheet CSV file for Nextflow RNA-seq pipeline."
    echo
    echo "Options:"
    echo "  -d, --dir DIR       Directory containing FASTQ files (default: current directory)"
    echo "  -o, --output FILE   Output CSV filename (default: samples.csv)"
    echo "  -p, --pattern PAT   FASTQ file pattern to match (default: *fastq.gz)"
    echo "  -h, --help          Display this help message and exit"
    echo
    echo "Output format: sample,fastq_1,fastq_2,strandedness"
    echo "All strandedness values will be set to 'auto'"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--dir)
            FASTQ_DIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        -p|--pattern)
            PATTERN="$2"
            shift 2
            ;;
        -h|--help)
            HELP=true
            shift
            ;;
        *)
            echo "Error: Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Display help message if requested
if [ "$HELP" = true ]; then
    usage
    exit 0
fi

# Convert FASTQ directory to absolute path
if [[ "$FASTQ_DIR" != /* ]]; then
    FASTQ_DIR="$(pwd)/$FASTQ_DIR"
fi
FASTQ_DIR=$(readlink -f "$FASTQ_DIR")

# Convert output file to absolute path if not already
if [[ "$OUTPUT_FILE" != /* ]]; then
    OUTPUT_FILE="$(pwd)/$OUTPUT_FILE"
fi

# Check if FASTQ directory exists
if [ ! -d "$FASTQ_DIR" ]; then
    echo "Error: Directory '$FASTQ_DIR' does not exist."
    exit 1
fi

# Write header to CSV file
echo "sample,fastq_1,fastq_2,strandedness" > "$OUTPUT_FILE"

# Track samples to avoid duplicates
declare -A SAMPLE_R1_FILES
declare -A SAMPLE_R2_FILES

echo "Searching for FASTQ files in '$FASTQ_DIR' matching pattern '$PATTERN'..."

# Find all FASTQ files matching the pattern
FASTQ_FILES=($(find "$FASTQ_DIR" -type f -name "*$PATTERN" | sort))

# First pass: Sort files into R1 and R2 groups by sample name
for FASTQ_FILE in "${FASTQ_FILES[@]}"; do
    # Extract base filename without path
    BASE_NAME=$(basename "$FASTQ_FILE")
    
    # Extract sample name based on different naming patterns
    if [[ "$BASE_NAME" =~ (.*)\.R1\. ]]; then
        # For files named like: SAMPLE.R1.fastq.gz
        SAMPLE_NAME="${BASH_REMATCH[1]}"
        SAMPLE_R1_FILES["$SAMPLE_NAME"]="$FASTQ_FILE"
    elif [[ "$BASE_NAME" =~ (.*)\.R2\. ]]; then
        # For files named like: SAMPLE.R2.fastq.gz
        SAMPLE_NAME="${BASH_REMATCH[1]}"
        SAMPLE_R2_FILES["$SAMPLE_NAME"]="$FASTQ_FILE"
    elif [[ "$BASE_NAME" =~ (.*)_R1 ]]; then
        # For files named like: SAMPLE_R1.fastq.gz
        SAMPLE_NAME="${BASH_REMATCH[1]}"
        SAMPLE_R1_FILES["$SAMPLE_NAME"]="$FASTQ_FILE"
    elif [[ "$BASE_NAME" =~ (.*)_R2 ]]; then
        # For files named like: SAMPLE_R2.fastq.gz
        SAMPLE_NAME="${BASH_REMATCH[1]}"
        SAMPLE_R2_FILES["$SAMPLE_NAME"]="$FASTQ_FILE"
    elif [[ "$BASE_NAME" =~ (.*)_1\. ]]; then
        # For files named like: SAMPLE_1.fastq.gz
        SAMPLE_NAME="${BASH_REMATCH[1]}"
        SAMPLE_R1_FILES["$SAMPLE_NAME"]="$FASTQ_FILE"
    elif [[ "$BASE_NAME" =~ (.*)_2\. ]]; then
        # For files named like: SAMPLE_2.fastq.gz
        SAMPLE_NAME="${BASH_REMATCH[1]}"
        SAMPLE_R2_FILES["$SAMPLE_NAME"]="$FASTQ_FILE"
    else
        # For single-end or other naming formats
        SAMPLE_NAME="${BASE_NAME%.*}"
        SAMPLE_NAME="${SAMPLE_NAME%.fastq}"
        SAMPLE_R1_FILES["$SAMPLE_NAME"]="$FASTQ_FILE"
    fi
done

# Second pass: Create the sample sheet
for SAMPLE_NAME in "${!SAMPLE_R1_FILES[@]}"; do
    R1_FILE="${SAMPLE_R1_FILES[$SAMPLE_NAME]}"
    R2_FILE="${SAMPLE_R2_FILES[$SAMPLE_NAME]:-}"
    
    if [[ -n "$R2_FILE" ]]; then
        # Paired-end sample
        echo "$SAMPLE_NAME,$R1_FILE,$R2_FILE,auto" >> "$OUTPUT_FILE"
    else
        # Single-end sample
        echo "$SAMPLE_NAME,$R1_FILE,,auto" >> "$OUTPUT_FILE"
    fi
done

# Count samples
SAMPLE_COUNT=${#SAMPLE_R1_FILES[@]}

echo "----------------------------------------"
echo "Summary:"
echo "  Created '$OUTPUT_FILE' with $SAMPLE_COUNT samples"
echo "  Found ${#SAMPLE_R1_FILES[@]} R1 files and ${#SAMPLE_R2_FILES[@]} R2 files"
echo "----------------------------------------"
echo "Example usage with Nextflow:"
echo "  nextflow run nf-core/rnaseq --input $OUTPUT_FILE [other options]"
echo "----------------------------------------"

# Show sample of the generated file
if [ $SAMPLE_COUNT -gt 0 ]; then
    echo "Preview of generated file:"
    head -n 3 "$OUTPUT_FILE"
    if [ $SAMPLE_COUNT -gt 2 ]; then
        echo "... plus $(($SAMPLE_COUNT - 2)) more samples"
    fi
fi

# Exit with error if no samples were found
if [ $SAMPLE_COUNT -eq 0 ]; then
    echo "Error: No FASTQ files found. Check your directory and pattern."
    exit 1
fi

exit 0
