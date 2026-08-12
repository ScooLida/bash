#!/bin/bash
# Build IQ-TREE gene trees and ASTRAL species trees for both BUSCO datasets.

set -euo pipefail

# Configuration
THREADS=9
WORK_DIR="./subset_parallel/pipeline_bulletproof_final"
MODEL="GTR+G"
ASTRAL_JAR="$HOME/Astral/astral.5.16.3.jar"
GAP_THRESHOLD=30
DATASET_WITHOUT="without_1k_3k_5k_5kS8"
DATASET_WITH_ALL="with_all_samples"
DATASETS=("$DATASET_WITHOUT" "$DATASET_WITH_ALL")

run_dataset() {
    local dataset=$1
    local align_dir="$WORK_DIR/final_alignments_${dataset}"
    local gene_list="$WORK_DIR/list_${GAP_THRESHOLD}_percent_gap_${dataset}.txt"
    local trees_dir="$WORK_DIR/gene_trees_${dataset}"
    local all_trees="$WORK_DIR/all_gene_trees_${dataset}.treefile"
    local astral_output="$WORK_DIR/astral_species_tree_${dataset}.tre"

    if [ ! -d "$align_dir" ] || [ ! -f "$gene_list" ]; then
        echo "Error: missing alignment directory or gene list for $dataset"
        return 1
    fi
    mkdir -p "$trees_dir"

    run_iqtree() {
        local gene=$1
        local input="$align_dir/${gene}.fasta"
        [ -f "$input" ] || return 0
        iqtree -s "$input" -st DNA -m "$MODEL" -nt 1 \
            --prefix "$trees_dir/$gene" -quiet
    }
    export -f run_iqtree
    export align_dir trees_dir MODEL
    grep -v '^$' "$gene_list" |
        xargs -P "$THREADS" -I {} bash -c 'run_iqtree "{}"'

    : > "$all_trees"
    for treefile in "$trees_dir"/*.treefile; do
        [ -e "$treefile" ] || continue
        cat "$treefile" >> "$all_trees"
    done
    if [ ! -s "$all_trees" ]; then
        echo "Error: no gene trees generated for $dataset"
        return 1
    fi

    if [ ! -f "$ASTRAL_JAR" ]; then
        echo "Error: ASTRAL jar not found: $ASTRAL_JAR"
        return 1
    fi
    java -jar "$ASTRAL_JAR" -i "$all_trees" -o "$astral_output" \
        2> "$WORK_DIR/astral_${dataset}.log"
    echo "Species tree saved to: $astral_output"
}

for dataset in "${DATASETS[@]}"; do run_dataset "$dataset"; done
