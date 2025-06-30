#!/bin/bash
# bed_validator.sh

BED_FILE=$1
echo "=== BED File Validation: $BED_FILE ==="

# Check basic format (6+ columns)
COLUMN_CHECK=$(awk '{if (NF < 6) print NR, NF}' $BED_FILE | head -5)
if [ -n "$COLUMN_CHECK" ]; then
    echo "ERROR: Lines with fewer than 6 columns found:"
    echo "$COLUMN_CHECK"
fi

# Check coordinates (end > start)
COORD_CHECK=$(awk '$3 <= $2 {print NR, $1, $2, $3}' $BED_FILE | head -5)
if [ -n "$COORD_CHECK" ]; then
    echo "ERROR: Lines with invalid coordinates (end <= start):"
    echo "$COORD_CHECK"
fi

# Check strands (must be +, -, or .)
STRAND_CHECK=$(awk '$6 !~ /^[+\-.]$/ {print NR, $6}' $BED_FILE | head -5)
if [ -n "$STRAND_CHECK" ]; then
    echo "ERROR: Lines with invalid strand (not +, -, or .):"
    echo "$STRAND_CHECK"
fi

# Check for chromosome problems
CHR_ISSUES=$(awk '{gsub(/[[:space:]]/, "", $1); if($1 == "") print NR}' $BED_FILE | head -5)
if [ -n "$CHR_ISSUES" ]; then
    echo "ERROR: Lines with empty chromosome names:"
    echo "$CHR_ISSUES"
fi

# RSeQC-specific checks
echo "Running RSeQC-specific validations..."

# Check for 0-length features after RSeQC processing
echo "Checking for potential RSeQC merge issues..."
MERGE_CHECK=$(bedtools merge -i $BED_FILE | awk '$3 <= $2 {print NR, $0}' | head -5)
if [ -n "$MERGE_CHECK" ]; then
    echo "WARNING: Potential RSeQC merge issues detected:"
    echo "$MERGE_CHECK"
fi

# Summary statistics
echo -e "\n=== BED File Statistics ==="
echo "Total entries: $(wc -l < $BED_FILE)"
echo "Unique chromosomes: $(cut -f1 $BED_FILE | sort | uniq | wc -l)"
echo "Unique feature names: $(cut -f4 $BED_FILE | sort | uniq | wc -l)"
echo "Bases covered: $(bedtools genomecov -i $BED_FILE -g <(cut -f1,3 $BED_FILE | sort -k1,1 -k2,2n | bedtools groupby -g 1 -c 2 -o max) | awk '$2 > 0 {sum += $3} END {print sum}')"
