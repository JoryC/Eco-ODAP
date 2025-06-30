#!/bin/bash

# format_count_matrix.sh - Reformat RNA-Seq count matrices

# Check if an input file is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 <input_count_matrix>"
    echo "Output will be saved as 'genes.data.tsv' in the current directory"
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="genes.data.tsv"

# Make sure the input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found"
    exit 1
fi

# Process the file using awk
awk 'BEGIN {FS="\t"; OFS="\t"}
     NR==1 {
         # Start from column 3 to skip gene_id and transcript_id/gene_name columns
         start_col = 3
         
         # Print header with empty first cell and sample IDs in quotes
         printf ""
         for (i=start_col; i<=NF; i++) {
             printf "%s\"%s\"", (i==start_col ? OFS : OFS), $i
         }
         printf "\n"
     }
     NR>1 {
         # Print gene_id in quotes followed by expression values (skip second column)
         printf "\"%s\"", $1
         for (i=start_col; i<=NF; i++) {
             printf "%s%s", OFS, $i
         }
         printf "\n"
     }' "$INPUT_FILE" > "$OUTPUT_FILE"

echo "Processed count matrix from '$INPUT_FILE' to '$OUTPUT_FILE'"

# Validate the output
echo "First few lines of the formatted count matrix:"
head -n 3 "$OUTPUT_FILE"