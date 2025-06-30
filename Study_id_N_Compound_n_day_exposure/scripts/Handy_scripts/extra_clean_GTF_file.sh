#!/bin/bash

# Enhanced script to clean GTF files and generate valid RSeQC-compatible BED files
# - Filters entries with missing transcript IDs
# - Optionally appends transcript versions to transcript IDs
# - Removes features with negative or zero-length after coordinate conversion
# - Standardizes chromosome names
# - Validates feature types
# - Generates a properly formatted BED file for RSeQC with strict validation

# Default options
APPEND_VERSIONS=false

# Parse command-line options
usage() {
  echo "Usage: $0 [-v] gtf_file"
  echo "  -v    Append transcript versions to transcript IDs"
  echo "  -h    Show this help message"
  exit 1
}

while getopts "vh" opt; do
  case $opt in
    v) APPEND_VERSIONS=true ;;
    h) usage ;;
    \?) usage ;;
  esac
done

# Shift the options so we can get the positional parameter
shift $((OPTIND-1))

# Check if argument is supplied
if [ -z "$1" ]; then
    echo "Please provide a GTF file as an argument."
    usage
fi

GTF=$1

# Create output filename
if [[ "$GTF" == *.* ]]; then
    ext="${GTF##*.}"
    base="${GTF%.*}"
    output_gtf="${base}_cleaned.${ext}"
else
    output_gtf="${GTF}_cleaned"
fi

# Temp files
sorted_gtf="${output_gtf}.sorted.tmp"
stats_file="${output_gtf}.stats.txt"

echo "=== Cleaning GTF file: $GTF ===" > $stats_file
echo "Started at $(date)" >> $stats_file
echo "Append versions to transcript IDs: $APPEND_VERSIONS" >> $stats_file

# Step 1: Process the GTF file with comprehensive cleaning
awk -v append_versions="$APPEND_VERSIONS" -F'\t' 'BEGIN{OFS="\t"; neg_count=0; no_transcript=0; chr_fixed=0;}
  # Print header lines and skip processing
  /^#/ {print; next}

  # Main processing block
  {
    line=$0;

    # Check for negative or zero length features (for RSeQC safety, need end > start + 1)
    if ($5 <= $4) {
      neg_count++;
      next; # Skip this line
    }

    # Standardize chromosome names (remove "chr" prefix if present)
    orig_chr = $1;
    if ($1 ~ /^chr/) {
      $1 = substr($1, 4);
      chr_fixed++;
    }
    
    # Remove any whitespace from chromosome names
    gsub(/[[:space:]]/, "", $1);

    # Check if the line has a transcript_id attribute
    if (match(line, /transcript_id "([^"]*)";/, tid)) {
      transcript_id = tid[1];

      # Only process lines where transcript_id is not empty
      if (transcript_id != "") {
        # Handle transcript versioning if enabled
        if (append_versions == "true") {
          # Check if transcript_version exists
          if (match(line, /transcript_version "([^"]*)";/, tver)) {
            transcript_version = tver[1];

            # Replace the transcript_id with versioned format
            sub(/transcript_id "[^"]*";/, "transcript_id \""transcript_id"."transcript_version"\";");
          }
        }

        # Print the cleaned line
        print;
      } else {
        no_transcript++;
      }
    } else {
      no_transcript++;
    }
  }

  END {
    print "Removed " neg_count " features with negative or zero lengths" > "/dev/stderr";
    print "Skipped " no_transcript " entries without valid transcript_id" > "/dev/stderr";
    print "Fixed " chr_fixed " chromosome names" > "/dev/stderr";
  }' "${GTF}" > "${output_gtf}"

# Step 2: Sort the GTF file
sort -k1,1 -k4,4n -k5,5n "${output_gtf}" > "${sorted_gtf}"
mv "${sorted_gtf}" "${output_gtf}"

# Step 3: Generate statistics for the cleaned file
echo "=== Statistics for Cleaned GTF ===" >> $stats_file
echo "Total entries: $(grep -v "^#" ${output_gtf} | wc -l)" >> $stats_file
echo "Unique chromosomes/scaffolds: $(cut -f1 ${output_gtf} | grep -v "^#" | sort | uniq | wc -l)" >> $stats_file
echo "Feature types:" >> $stats_file
cut -f3 ${output_gtf} | grep -v "^#" | sort | uniq -c | sort -nr >> $stats_file
echo "Transcript count: $(grep -c 'transcript_id' ${output_gtf})" >> $stats_file

if [ "$APPEND_VERSIONS" = true ]; then
  echo "Versioned transcripts: $(grep -E 'transcript_id "[^"]+\.[0-9]+"' ${output_gtf} | wc -l)" >> $stats_file
fi

echo "Finished at $(date)" >> $stats_file

echo "Created cleaned GTF file at ${output_gtf}"
echo "Cleaning statistics saved to ${stats_file}"

# Step 4: Generate a BED file for RSeQC with strict validation
output_bed="${base}_cleaned.bed"
echo "Generating BED file for RSeQC..."

# Extract exon features and convert to BED format with careful validation
awk -F'\t' 'BEGIN{OFS="\t"; valid=0; skipped=0;}
  $3 == "exon" {
    # Calculate BED coordinates (0-based start, 1-based end)
    start = $4 - 1;
    end = $5;
    
    # Strict validation to ensure positive length after conversion
    if (end <= start) {
      skipped++;
      next;
    }
    
    # Extract gene_id
    if (!match($9, /gene_id "([^"]+)";/, gid)) {
      skipped++;
      next;
    }
    gene_id = gid[1];
    
    # Extract transcript_id
    if (!match($9, /transcript_id "([^"]+)";/, tid)) {
      skipped++;
      next;
    }
    transcript_id = tid[1];
    
    # Create a valid name field without special characters that might confuse RSeQC
    name = gene_id"_"transcript_id;
    gsub(/[^a-zA-Z0-9_.-]/, "_", name);
    
    # BED format: chrom, start, end, name, score, strand
    # Use a fixed score of 0 for better compatibility
    valid++;
    print $1, start, end, name, "0", $7;
  }
  END {
    print "Valid BED entries: " valid > "/dev/stderr";
    print "Skipped entries: " skipped > "/dev/stderr";
  }' "${output_gtf}" | \
  sort -k1,1 -k2,2n | \
  # Final validation pass - ensure all entries have positive length
  awk 'BEGIN{OFS="\t"} $3 > $2 {print}' > "${output_bed}"

# Final verification of the BED file
echo "=== BED File Validation ===" >> $stats_file
echo "Total BED entries: $(wc -l < ${output_bed})" >> $stats_file
echo "Unique chromosomes: $(cut -f1 ${output_bed} | sort | uniq | wc -l)" >> $stats_file
invalid_count=$(awk 'BEGIN{c=0} $3 <= $2 {c++} END{print c}' "${output_bed}")
echo "Invalid entries (should be 0): ${invalid_count}" >> $stats_file

if [ "$invalid_count" -gt 0 ]; then
    echo "WARNING: BED file contains ${invalid_count} invalid entries. RSeQC may fail." | tee -a $stats_file
else 
    echo "BED file validation passed. All entries have positive length." | tee -a $stats_file
fi

echo "Created BED file at ${output_bed}"
echo "Additional validation statistics in ${stats_file}"
