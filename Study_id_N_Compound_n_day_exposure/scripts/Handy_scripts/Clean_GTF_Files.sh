#!/bin/bash

# Create a cleaned GTF file with simplified attributes
fixed_gtf="${GENOMEDIR}/fixed_annotation.gtf"

# Keep header lines and ensure proper attribute formatting
awk -F'\t' 'BEGIN{OFS="\t"} 
  /^#/ {print; next} 
  {
    line=$0;
    if ($3 == "transcript" || $3 == "exon") {
      # Extract the gene_id and transcript_id
      match(line, /gene_id "([^"]*)";/, gid);
      gene_id = gid[1];
      match(line, /transcript_id "([^"]*)";/, tid);
      transcript_id = tid[1];
      
      if (transcript_id != "") {
        # Reformat the attributes to ensure compatible formatting
        attrs = "gene_id \"" gene_id "\"; transcript_id \"" transcript_id "\";";
        
        # Print the line with simplified attributes
        print $1,$2,$3,$4,$5,$6,$7,$8,attrs;
      } else if ($3 == "transcript") {
        print "WARNING: Found transcript without ID at line: " line > "/dev/stderr";
      }
    } else {
      # For other feature types, just print the line
      print line;
    }
  }' "${GTF}" > "${fixed_gtf}"

echo "Created fixed GTF file at ${fixed_gtf}"