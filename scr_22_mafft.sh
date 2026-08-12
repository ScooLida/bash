#!/bin/bash
# Align the two BUSCO FASTA datasets with MAFFT.

set -euo pipefail

# Configuration
THREADS=9
WORK_DIR="./subset_parallel/pipeline_bulletproof_final"
GENE_LIST="$WORK_DIR/extracted_genes_list.txt"
DATASET_WITHOUT="without_1k_3k_5k_5kS8"
DATASET_WITH_ALL="with_all_samples"
DATASETS=("$DATASET_WITHOUT" "$DATASET_WITH_ALL")

if [ ! -f "$GENE_LIST" ]; then
    echo "Error: gene list not found: $GENE_LIST"
    exit 1
fi

run_dataset() {
    local dataset=$1
    local source_dir="$WORK_DIR/combined_unaligned_${dataset}"
    local output_dir="$WORK_DIR/final_alignments_${dataset}"

    if [ ! -d "$source_dir" ]; then
        echo "Error: source directory not found: $source_dir"
        return 1
    fi
    mkdir -p "$output_dir"

    run_mafft() {
        local gene=$1
        local input="$source_dir/${gene}.fasta"
        local output="$output_dir/${gene}.fasta"
        [ -f "$input" ] || return 0
        mafft --auto "$input" > "$output" 2>/dev/null
    }
    export -f run_mafft
    export source_dir output_dir

    grep -v '^$' "$GENE_LIST" |
        xargs -P "$THREADS" -I {} bash -c 'run_mafft "{}"'
    echo "Alignment complete: $dataset"
}

for dataset in "${DATASETS[@]}"; do run_dataset "$dataset"; done
