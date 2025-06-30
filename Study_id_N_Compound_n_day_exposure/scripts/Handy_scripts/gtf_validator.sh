#!/bin/bash
# gtf_validator.sh

GTF_FILE=$1
echo "=== GTF File Validation: $GTF_FILE ==="

# Check basic format (9 columns)
COLUMN_CHECK=$(awk -F'\t' '{if (NF != 9) print NR, NF}' $GTF_FILE | head -5)
if [ -n "$COLUMN_CHECK" ]; then
    echo "ERROR: Lines with incorrect number of columns:"
    echo "$COLUMN_CHECK"
fi

# Check coordinates (end >= start)
COORD_CHECK=$(awk -F'\t' '$5 < $4 {print NR, $1, $4, $5}' $GTF_FILE | head -5)
if [ -n "$COORD_CHECK" ]; then
    echo "ERROR: Lines with invalid coordinates (end < start):"
    echo "$COORD_CHECK"
fi

# Check essential attributes
ATTR_CHECK=$(awk -F'\t' '$3 == "exon" && !($9 ~ /gene_id/) {print NR}' $GTF_FILE | head -5)
if [ -n "$ATTR_CHECK" ]; then
    echo "ERROR: Exon entries missing gene_id attribute:"
    echo "$ATTR_CHECK"
fi

TRANSCRIPT_CHECK=$(awk -F'\t' '$3 == "exon" && !($9 ~ /transcript_id/) {print NR}' $GTF_FILE | head -5)
if [ -n "$TRANSCRIPT_CHECK" ]; then
    echo "ERROR: Exon entries missing transcript_id attribute:"
    echo "$TRANSCRIPT_CHECK"
fi

# Check for feature consistency
echo -e "\n=== Feature Consistency Check ==="
echo "Genes: $(grep -c $'\tgene\t' $GTF_FILE)"
echo "Transcripts: $(grep -c $'\ttranscript\t' $GTF_FILE)"
echo "Exons: $(grep -c $'\texon\t' $GTF_FILE)"
echo "CDS: $(grep -c $'\tCDS\t' $GTF_FILE)"

# Find orphan exons (exons without parent transcripts)
echo -e "\nChecking for exons without parent transcripts..."
awk -F'\t' '$3=="exon" {match($9, /transcript_id "([^"]+)"/, tid); print tid[1]}' $GTF_FILE | \
  sort | uniq > exons.tmp
awk -F'\t' '$3=="transcript" {match($9, /transcript_id "([^"]+)"/, tid); print tid[1]}' $GTF_FILE | \
  sort | uniq > transcripts.tmp
comm -23 exons.tmp transcripts.tmp > orphans.tmp

if [ -s orphans.tmp ]; then
    echo "WARNING: Found $(wc -l < orphans.tmp) exons without parent transcripts"
    echo "First 5 orphan transcript IDs:"
    head -5 orphans.tmp
fi
