#!/bin/bash

# Create a temporary file to store results
temp_file=$(mktemp)

echo -e "Sample\t%Unique\t%Multi-mapped\t%Mapped\t%True-Unique\t%True-Multi" > "$temp_file"

# Process each ERR directory
for dir in ERR*; do
  if [ -d "$dir" ] && [ -f "$dir/aux_info/ambig_info.tsv" ] && [ -f "$dir/aux_info/meta_info.json" ]; then
    # Extract counts
    unique_sum=$(awk 'NR>1 {sum+=$1} END {print sum}' "$dir/aux_info/ambig_info.tsv")
    ambig_sum=$(awk 'NR>1 {sum+=$2} END {print sum}' "$dir/aux_info/ambig_info.tsv")
    total=$((unique_sum + ambig_sum))
    
    # Get percent_mapped from meta_info.json
    percent_mapped=$(grep -o '"percent_mapped":[^,]*' "$dir/aux_info/meta_info.json" | cut -d':' -f2)
    # Remove any quotes or whitespace
    percent_mapped=$(echo "$percent_mapped" | tr -d ' "')
    
    # Calculate percentages
    if [ $total -gt 0 ]; then
      unique_pct=$(echo "scale=2; 100*$unique_sum/$total" | bc)
      ambig_pct=$(echo "scale=2; 100*$ambig_sum/$total" | bc)
      
      # Calculate true percentages (relative to all reads, not just mapped ones)
      true_unique_pct=$(echo "scale=2; $unique_pct*$percent_mapped/100" | bc)
      true_ambig_pct=$(echo "scale=2; $ambig_pct*$percent_mapped/100" | bc)
    else
      unique_pct="0.00"
      ambig_pct="0.00"
      true_unique_pct="0.00"
      true_ambig_pct="0.00"
    fi
    
    echo -e "$dir\t$unique_pct\t$ambig_pct\t$percent_mapped\t$true_unique_pct\t$true_ambig_pct" >> "$temp_file"
  else
    echo "Warning: Missing required files for $dir" >&2
  fi
done

# Display the results in a nice table
column -t "$temp_file"

# Calculate summary statistics
echo -e "\nSummary Statistics:"
tail -n +2 "$temp_file" | awk '
  BEGIN {
    min_unique=101; max_unique=0; 
    min_ambig=101; max_ambig=0;
    min_true_unique=101; max_true_unique=0;
    min_true_ambig=101; max_true_ambig=0;
    min_mapped=101; max_mapped=0;
    count=0
  } 
  {
    count++;
    if ($2 < min_unique) min_unique=$2;
    if ($2 > max_unique) max_unique=$2;
    if ($3 < min_ambig) min_ambig=$3;
    if ($3 > max_ambig) max_ambig=$3;
    if ($4 < min_mapped) min_mapped=$4;
    if ($4 > max_mapped) max_mapped=$4;
    if ($5 < min_true_unique) min_true_unique=$5;
    if ($5 > max_true_unique) max_true_unique=$5;
    if ($6 < min_true_ambig) min_true_ambig=$6;
    if ($6 > max_true_ambig) max_true_ambig=$6;
  } 
  END {
    if (count == 0) {
      print "No data found. Check if paths to required files are correct.";
    } else {
      printf "Original calculations (relative to mapped reads only):\n";
      printf "  Uniquely mapped reads: min=%.2f%%, max=%.2f%%\n", min_unique, max_unique;
      printf "  Multi-mapped reads:    min=%.2f%%, max=%.2f%%\n", min_ambig, max_ambig;
      printf "Total mapped reads:      min=%.2f%%, max=%.2f%%\n", min_mapped, max_mapped;
      printf "True calculations (relative to all reads):\n";
      printf "  Uniquely mapped reads: min=%.2f%%, max=%.2f%%\n", min_true_unique, max_true_unique;
      printf "  Multi-mapped reads:    min=%.2f%%, max=%.2f%%\n", min_true_ambig, max_true_ambig;
    }
  }'

# Clean up
rm "$temp_file"