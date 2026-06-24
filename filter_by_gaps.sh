#!/bin/bash
# The script scans a folder with aligned genes in FASTA format,
# counts gaps (dashes) and reports how many genes remain when using a
# "minimum __ percent of nucleotides" filter across the whole file.

# Input directory
INPUT_DIR="./astral_pipeline/perfect_aligned_genes"


echo "=== Alignment quality analysis ==="

# Initialize counters
total=0
c65=0
c50=0
c40=0
c30=0
c20=0
c10=0
c5=0

for file in "$INPUT_DIR"/*.fasta; do
    [ -e "$file" ] || continue

    # 1. Count total nucleotides (ignoring header lines starting with '>')
    total_chars=$(grep -v "^>" "$file" | tr -d '\n' | wc -c)

    # 2. Count dashes ('-')
    gaps=$(grep -v "^>" "$file" | tr -d -c '-' | wc -c)

    # Guard against division by zero (in case the file is empty)
    if [ "$total_chars" -eq 0 ]; then continue; fi

    # 3. Calculate gap percentage
    pct=$(( gaps * 100 / total_chars ))

    # 4. Assign to categories
    ((total++))
    if [ "$pct" -le 65 ]; then ((c65++)); fi
    if [ "$pct" -le 50 ]; then ((c50++)); fi
    if [ "$pct" -le 40 ]; then ((c40++)); fi
    if [ "$pct" -le 30 ]; then ((c30++)); fi
    if [ "$pct" -le 20 ]; then ((c20++)); fi
    if [ "$pct" -le 10 ]; then ((c10++)); fi
    if [ "$pct" -le 5 ];  then ((c5++));  fi
done

echo ""
echo "------------------------------------------------"
echo "Total genes analyzed: $total"
echo "------------------------------------------------"
echo "Remaining at threshold <= 65% gaps:  $c65"
echo "Remaining at threshold <= 50% gaps:  $c50"
echo "Remaining at threshold <= 40% gaps:  $c40"
echo "------------------------------------------------"
echo "Remaining at threshold <= 30% gaps:  $c30"
echo "Remaining at threshold <= 20% gaps:  $c20"
echo "Remaining at threshold <= 10% gaps:  $c10"
echo "Remaining at threshold <= 5%  gaps:  $c5"
echo "------------------------------------------------"
