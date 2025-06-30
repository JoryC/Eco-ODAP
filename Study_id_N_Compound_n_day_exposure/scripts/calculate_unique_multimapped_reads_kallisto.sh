#!/bin/bash

# Create a temporary file to store results
temp_file=$(mktemp)

echo -e "Sample\t%Unique\t%Multi-mapped" > "$temp_file"

# Process each ERR directory
for dir in ERR*; do
  if [ -d "$dir" ] && [ -f "$dir/run_info.json" ]; then
    # Extract counts using jq if available, or grep+sed if not
    if command -v jq &> /dev/null; then
      # Using jq for proper JSON parsing (recommended)
      p_unique=$(jq -r '.p_unique // 0' "$dir/run_info.json")
      p_pseudoaligned=$(jq -r '.p_pseudoaligned // 0' "$dir/run_info.json")
      
      # Check if these values exist in the JSON
      if [[ "$p_unique" == "null" || "$p_pseudoaligned" == "null" ]]; then
        # Fallback for older Kallisto versions
        n_unique=$(jq -r '.n_unique // 0' "$dir/run_info.json")
        n_processed=$(jq -r '.n_processed // 0' "$dir/run_info.json")
        n_pseudoaligned=$(jq -r '.n_pseudoaligned // 0' "$dir/run_info.json")
        
        if [[ "$n_processed" != "0" && "$n_processed" != "null" ]]; then
          p_unique=$(echo "scale=2; 100*$n_unique/$n_processed" | bc)
          p_pseudoaligned=$(echo "scale=2; 100*$n_pseudoaligned/$n_processed" | bc)
        else
          p_unique="0.00"
          p_pseudoaligned="0.00"
        fi
      fi
    else
      # Fallback to grep+sed if jq is not available
      p_unique=$(grep -o '"p_unique":[^,]*' "$dir/run_info.json" | sed 's/"p_unique"://')
      p_pseudoaligned=$(grep -o '"p_pseudoaligned":[^,]*' "$dir/run_info.json" | sed 's/"p_pseudoaligned"://')
      
      # If p_unique doesn't exist, calculate from n_unique and n_processed
      if [ -z "$p_unique" ]; then
        n_unique=$(grep -o '"n_unique":[^,]*' "$dir/run_info.json" | sed 's/"n_unique"://')
        n_processed=$(grep -o '"n_processed":[^,]*' "$dir/run_info.json" | sed 's/"n_processed"://')
        n_pseudoaligned=$(grep -o '"n_pseudoaligned":[^,]*' "$dir/run_info.json" | sed 's/"n_pseudoaligned"://')
        
        if [ -n "$n_processed" ] && [ "$n_processed" != "0" ]; then
          p_unique=$(echo "scale=2; 100*$n_unique/$n_processed" | bc)
          p_pseudoaligned=$(echo "scale=2; 100*$n_pseudoaligned/$n_processed" | bc)
        else
          p_unique="0.00"
          p_pseudoaligned="0.00"
        fi
      fi
    fi
    
    # Calculate multi-mapped percentage
    p_multi=$(echo "scale=2; $p_pseudoaligned - $p_unique" | bc)
    
    # Handle negative values (could happen due to rounding errors)
    if (( $(echo "$p_multi < 0" | bc -l) )); then
      p_multi="0.00"
    fi
    
    echo -e "$dir\t$p_unique\t$p_multi" >> "$temp_file"
  fi
done

# Display the results in a nice table
column -t "$temp_file"

# Calculate summary statistics
echo -e "\nSummary Statistics:"
tail -n +2 "$temp_file" | awk '
  BEGIN {min_unique=101; max_unique=0; min_multi=101; max_multi=0; count=0} 
  {
    count++;
    if ($2 < min_unique) min_unique=$2;
    if ($2 > max_unique) max_unique=$2;
    if ($3 < min_multi) min_multi=$3;
    if ($3 > max_multi) max_multi=$3;
  } 
  END {
    if (count == 0) {
      print "No data found. Check if paths to run_info.json are correct.";
    } else {
      printf "Uniquely mapped reads: min=%.2f%%, max=%.2f%%\n", min_unique, max_unique;
      printf "Multi-mapped reads:    min=%.2f%%, max=%.2f%%\n", min_multi, max_multi;
    }
  }'

# Clean up
rm "$temp_file"