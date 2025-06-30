#!/bin/bash

# Check if argument is supplied
if [ -z "$1" ]; then
    echo "Please provide a GTF file as an argument."
    exit 1
fi

if [ -z "$2" ]; then
    echo "Please provide an output filename as an argument."
    exit 1
fi

# Extract fields from GTF file
gtf_file=$1
output_file=$2

# Use awk to extract relevant information
awk '
    BEGIN { FS = "\t"; OFS = "\t" }
    $3 == "transcript" {
        gene_id = ""; transcript_id = ""; gene_name = ""; gene_biotype = "";
        for (i = 1; i <= NF; i++) {
            if ($i ~ /gene_id/) {
                match($i, /gene_id "([^"]+)";/, arr);
                if (arr[1] != "") gene_id = arr[1];
            }
            if ($i ~ /transcript_id/) {
                match($i, /transcript_id "([^"]+)";/, arr);
                if (arr[1] != "") transcript_id = arr[1];
            }
            if ($i ~ /gene_name/ || $i ~ /gene /) {
                match($i, /(gene_name|gene) "([^"]+)";/, arr);
                if (arr[2] != "") gene_name = arr[2];
            }
            if ($i ~ /transcript_biotype/ || $i ~ /gene_biotype/) {
                match($i, /(transcript_biotype|gene_biotype) "([^"]+)";/, arr);
                if (arr[2] != "") gene_biotype = arr[2];
            }
        }
        if(gene_id != "" && transcript_id != "") {
            if(gene_name == "") gene_name = "NA";
            if(gene_biotype == "") gene_biotype = "NA";
            print transcript_id, gene_id, gene_name, gene_biotype;
        }
    }
' $gtf_file > $output_file

echo "Extraction completed. Output saved to $output_file."
