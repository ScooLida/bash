#!/bin/bash
# Filter aligned BUSCO genes by gap percentage and required sample content.
# The no-excluded dataset requires 4k; the all-sample dataset requires
# 1k, 3k, 4k, and 5kS8.

set -euo pipefail

# Configuration
WORK_DIR="./subset_parallel/pipeline_bulletproof_final"
THRESHOLD_FOR_LIST=30
MAX_SEQ_LENGTH=5000
MIN_ATGC_PERCENT=30
MIN_NUCLEOTIDES_PER_SAMPLE=30
DATASET_WITHOUT="without_1k_3k_5k_5kS8"
DATASET_WITH_ALL="with_all_samples"
DATASETS=("$DATASET_WITHOUT" "$DATASET_WITH_ALL")
REQUIRED_SAMPLES_WITHOUT=""
REQUIRED_SAMPLES_WITH_ALL="1k,3k,5k,5kS8"

filter_dataset() {
    local dataset=$1
    local required_samples=$2
    local input_dir="$WORK_DIR/final_alignments_${dataset}"
    local output_list="$WORK_DIR/list_${THRESHOLD_FOR_LIST}_percent_gap_${dataset}.txt"

    if [ ! -d "$input_dir" ]; then
        echo "Error: alignment directory not found: $input_dir"
        return 1
    fi
    : > "$output_list"

    for file in "$input_dir"/*.fasta; do
        [ -e "$file" ] || continue
        gene=$(basename "$file" .fasta)
        result=$(awk -v max_len="$MAX_SEQ_LENGTH" \
            -v min_atgc="$MIN_ATGC_PERCENT" \
            -v min_nuc="$MIN_NUCLEOTIDES_PER_SAMPLE" \
            -v max_gaps="$THRESHOLD_FOR_LIST" \
            -v req_samples="$required_samples" '
            BEGIN {
                if (req_samples != "") {
                    split(req_samples, required_list, ",")
                    for (i in required_list) required[required_list[i]] = 1
                }
            }
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
                    if (seq_len > max_len) failed = 1
                    all_length += seq_len
                    gsub(/[^ATGCatgc]/, "", sequence)
                    atgc = length(sequence)
                    all_atgc += atgc
                    if (atgc < min_nuc) failed = 1
                    if (sample in required && (seq_len == 0 || atgc * 100 / seq_len < min_atgc)) failed = 1
                }

                for (sample in required) if (!(sample in sequences)) failed = 1
                gap_percent = all_length > 0 ? 100 - all_atgc * 100 / all_length : 100
                if (gap_percent > max_gaps) failed = 1
                if (failed == 0) print "PASS"
            }
        ' "$file")

        if [ "$result" = "PASS" ]; then printf '%s\n' "$gene" >> "$output_list"; fi
    done
    sort -u "$output_list" -o "$output_list"
    echo "Genes passing $dataset: $(wc -l < "$output_list")"
    echo "List saved to: $output_list"
}

filter_dataset "$DATASET_WITHOUT" "$REQUIRED_SAMPLES_WITHOUT"
filter_dataset "$DATASET_WITH_ALL" "$REQUIRED_SAMPLES_WITH_ALL"
