#!/bin/bash
# Select loci using modern samples only. The same locus list is then used for
# the modern alignments and for the alignments containing ancient samples.

set -euo pipefail

# Configuration
WORK_DIR="./subset_parallel/pipeline_bulletproof_final"
THRESHOLD_FOR_LIST=30
MAX_SEQ_LENGTH=5000
MIN_ATGC_PERCENT=30
MIN_NUCLEOTIDES_PER_SAMPLE=30
DATASET_MODERN="modern"
DATASET_WITH_ALL="with_all_samples"
MODERN_LIST="$WORK_DIR/list_${THRESHOLD_FOR_LIST}_percent_gap_${DATASET_MODERN}.txt"
ALL_LIST="$WORK_DIR/list_${THRESHOLD_FOR_LIST}_percent_gap_${DATASET_WITH_ALL}.txt"

MODERN_DIR="$WORK_DIR/final_alignments_${DATASET_MODERN}"
ALL_DIR="$WORK_DIR/final_alignments_${DATASET_WITH_ALL}"

if [ ! -d "$MODERN_DIR" ]; then
    echo "Error: modern alignment directory not found: $MODERN_DIR" >&2
    exit 1
fi
if [ ! -d "$ALL_DIR" ]; then
    echo "Error: all-sample alignment directory not found: $ALL_DIR" >&2
    exit 1
fi

: > "$MODERN_LIST"
for file in "$MODERN_DIR"/*.fasta; do
    [ -e "$file" ] || continue
    gene=$(basename "$file" .fasta)
    result=$(awk -v max_len="$MAX_SEQ_LENGTH" \
        -v min_atgc="$MIN_ATGC_PERCENT" \
        -v min_nuc="$MIN_NUCLEOTIDES_PER_SAMPLE" \
        -v max_gaps="$THRESHOLD_FOR_LIST" '
        /^>/ {
            if (id != "") sequences[id] = sequence
            id = substr($1, 2)
            sequence = ""
            next
        }
        { sequence = sequence $0 }
        END {
            if (id != "") sequences[id] = sequence
            all_length = 0
            all_atgc = 0
            failed = 0

            for (sample in sequences) {
                sequence = sequences[sample]
                seq_len = length(sequence)
                if (seq_len == 0 || seq_len > max_len) failed = 1
                all_length += seq_len
                gsub(/[^ATGCatgc]/, "", sequence)
                atgc = length(sequence)
                all_atgc += atgc
                if (atgc < min_nuc) failed = 1
                if (seq_len > 0 && atgc * 100 / seq_len < min_atgc) failed = 1
            }

            gap_percent = all_length > 0 ? 100 - all_atgc * 100 / all_length : 100
            if (gap_percent > max_gaps) failed = 1
            if (failed == 0) print "PASS"
        }
    ' "$file")

    if [ "$result" = "PASS" ]; then
        printf '%s\n' "$gene" >> "$MODERN_LIST"
    fi
done
sort -u "$MODERN_LIST" -o "$MODERN_LIST"

if [ ! -s "$MODERN_LIST" ]; then
    echo "Error: no modern loci passed filtering." >&2
    exit 1
fi

# Ancient samples do not determine locus inclusion. They are retained wherever
# the corresponding modern-selected alignment exists.
while read -r gene; do
    [ -z "$gene" ] && continue
    if [ ! -s "$ALL_DIR/${gene}.fasta" ]; then
        echo "Error: all-sample alignment missing for modern locus: $gene" >&2
        exit 1
    fi
done < "$MODERN_LIST"
cp "$MODERN_LIST" "$ALL_LIST"

echo "Modern loci passing filters: $(wc -l < "$MODERN_LIST")"
echo "Modern locus list: $MODERN_LIST"
echo "Shared locus list for ancient inclusion: $ALL_LIST"
