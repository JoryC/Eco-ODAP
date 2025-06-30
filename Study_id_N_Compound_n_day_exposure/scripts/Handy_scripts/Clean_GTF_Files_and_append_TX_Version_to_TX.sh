#!/bin/bash

# Create a filtered GTF file that skips any entries with missing transcript IDs
# and appends transcript versions to transcript IDs

# Check if argument is supplied
if [ -z "$1" ]; then
    echo "Please provide a GTF file as an argument."
    exit 1
fi

GTF=$1

# Create output filename by appending "_filtered_versioned" to original filename
if [[ "$GTF" == *.* ]]; then
    # Get the base name and extension
    ext="${GTF##*.}"
    base="${GTF%.*}"
    # Create new filename with _filtered_versioned before the extension
    output_gtf="${base}_filtered_versioned.${ext}"
else
    # No extension case
    output_gtf="${GTF}_filtered_versioned"
fi

# Process the GTF file: filter entries without transcript_id and append versions
awk -F'\t' 'BEGIN{OFS="\t"}
  /^#/ {print; next}  # Print header lines and skip to next line
  {
    line=$0;
    
    # Check if the line has a transcript_id attribute
    if (match(line, /transcript_id "([^"]*)";/, tid)) {
      transcript_id = tid[1];

      # Only process lines where transcript_id is not empty
      if (transcript_id != "") {
        # Check if transcript_version exists
        if (match(line, /transcript_version "([^"]*)";/, tver)) {
          transcript_version = tver[1];
          
          # Replace the transcript_id with versioned format
          versioned_line = gensub(/transcript_id "([^"]*)";/, "transcript_id \"\\1."transcript_version"\";", "g", line);
          print versioned_line;
        } else {
          # No version found, use the original line
          print line;
        }
      } else {
        print "WARNING: Skipping line with empty transcript_id: " line > "/dev/stderr";
      }
    } else {
      # No transcript_id attribute found at all
      print "WARNING: Skipping line without transcript_id attribute: " line > "/dev/stderr";
    }
  }' "${GTF}" > "${output_gtf}"

echo "Created filtered and versioned GTF file at ${output_gtf}"

# Count processed entries
TRANSCRIPT_COUNT=$(grep -c 'transcript_id' "${output_gtf}")
VERSIONED_COUNT=$(grep -E 'transcript_id "[^"]+\.[0-9]+"' "${output_gtf}" | wc -l)

echo "Processed ${TRANSCRIPT_COUNT} transcript entries"
echo "Added version numbers to approximately ${VERSIONED_COUNT} transcript IDs"