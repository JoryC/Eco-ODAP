#!/bin/bash

# RSeQC-Compatible BED File Generator
# Specially designed to avoid the "Count (-1) must be non-negative" error

if [ -z "$1" ]; then
    echo "Please provide a GTF file as an argument."
    exit 1
fi

GTF_FILE=$1
OUTPUT_BED="${GTF_FILE%.*}_rseqc_safe.bed"

echo "Generating RSeQC-compatible BED file from $GTF_FILE..."

# Extract exons and convert to BED format with strict validation
awk -F'\t' 'BEGIN{OFS="\t"}
  # Only process exon features and non-comment lines
  $3 == "exon" && !/^#/ {
    # Ensure end > start (after adjusting for 0-based BED coordinates)
    start = $4 - 1
    end = $5
    
    # Skip if end <= start (would create invalid range)
    if (end <= start) {
      print "Skipping invalid range at line: " NR > "/dev/stderr"
      next
    }
    
    # Extract gene_id
    if (match($9, /gene_id "([^"]+)";/, gid)) {
      gene_id = gid[1]
      
      # Extract transcript_id
      if (match($9, /transcript_id "([^"]+)";/, tid)) {
        transcript_id = tid[1]
        
        # Validate chromosome name (remove any whitespace)
        chrom = $1
        gsub(/[[:space:]]/, "", chrom)
        
        # Create a valid and safe BED entry
        # BED format: chrom, start, end, name, score, strand
        print chrom, start, end, gene_id"_"transcript_id, "0", $7
      }
    }
  }' "${GTF_FILE}" | \
  # Sort by chromosome and start position
  sort -k1,1 -k2,2n | \
  # Final validation - ensure start < end
  awk 'BEGIN{OFS="\t"} $3 > $2 {print}' > "${OUTPUT_BED}"

# Verify the BED file
echo "Validating BED file..."
invalid_count=$(awk 'BEGIN{c=0} $3 <= $2 {c++} END{print c}' "${OUTPUT_BED}")

if [ "$invalid_count" -gt 0 ]; then
    echo "ERROR: Found $invalid_count lines with invalid ranges in the output BED file."
    exit 1
else
    echo "SUCCESS: Generated valid BED file at ${OUTPUT_BED}"
    echo "Total entries: $(wc -l < ${OUTPUT_BED})"
fi
