#!/bin/bash
# Align modern BUSCO loci first, then add ancient sequences to the fixed
# modern alignments without changing the modern alignment columns. Missing
# ancient sequences are represented by Ns across the full alignment.

set -euo pipefail

# Configuration
THREADS=9
WORK_DIR="./subset_parallel/pipeline_bulletproof_final"
GENE_LIST="$WORK_DIR/extracted_genes_list.txt"
ANCIENT_SAMPLES="1k,3k,4k,5kS8"
DATASET_MODERN="modern"
DATASET_WITH_ALL="with_all_samples"
SCRIPT_PATH="$(readlink -f "$0")"

align_gene() {
    local gene=$1
    local modern_source="$WORK_DIR/combined_unaligned_${DATASET_MODERN}/${gene}.fasta"
    local all_source="$WORK_DIR/combined_unaligned_${DATASET_WITH_ALL}/${gene}.fasta"
    local modern_output_dir="$WORK_DIR/final_alignments_${DATASET_MODERN}"
    local all_output_dir="$WORK_DIR/final_alignments_${DATASET_WITH_ALL}"
    local modern_output="$modern_output_dir/${gene}.fasta"
    local all_output="$all_output_dir/${gene}.fasta"

    if [ ! -s "$modern_source" ]; then
        # A locus without modern sequence cannot define the modern backbone.
        echo "Skipping locus without modern sequence: $gene" >&2
        return 0
    fi
    if [ ! -s "$all_source" ]; then
        echo "Error: all-sample FASTA not found or empty: $all_source" >&2
        return 1
    fi
    IFS=',' read -r -a ancient_samples <<< "$ANCIENT_SAMPLES"
    all_ancient_present=true
    if [ -s "$all_output" ]; then
        for sample in "${ancient_samples[@]}"; do
            if ! grep -qE "^>${sample}([[:space:]]|$)" "$all_output"; then
                all_ancient_present=false
                break
            fi
        done
    else
        all_ancient_present=false
    fi
    if [ -s "$modern_output" ] && [ -s "$all_output" ] && "$all_ancient_present"; then
        return 0
    fi

    mkdir -p "$modern_output_dir" "$all_output_dir"
    mafft --auto "$modern_source" > "$modern_output" 2>/dev/null

    ancient_file=$(mktemp "$WORK_DIR/.ancient_fragments.XXXXXX.fasta")
    trap 'rm -f "$ancient_file"' EXIT
    awk -v ancient="$ANCIENT_SAMPLES" '
        BEGIN {
            split(ancient, names, ",")
            for (i in names) ancient_sample[names[i]] = 1
        }
        /^>/ {
            sample = $1
            sub(/^>/, "", sample)
            keep = (sample in ancient_sample)
        }
        keep { print }
    ' "$all_source" > "$ancient_file"

    if [ -s "$ancient_file" ]; then
        mafft --addfragments "$ancient_file" "$modern_output" > "$all_output" 2>/dev/null
    else
        cp "$modern_output" "$all_output"
    fi

    alignment_length=$(awk '!/^>/{print length($0); exit}' "$modern_output")
    if [ -z "$alignment_length" ] || [ "$alignment_length" -le 0 ]; then
        echo "Error: could not determine alignment length: $modern_output" >&2
        return 1
    fi

    for sample in "${ancient_samples[@]}"; do
        if ! grep -qE "^>${sample}([[:space:]]|$)" "$all_output"; then
            printf ">%s\n" "$sample" >> "$all_output"
            printf '%*s\n' "$alignment_length" '' | tr ' ' 'N' >> "$all_output"
        fi
    done
}

if [ "${1:-}" = "__align_gene" ]; then
    align_gene "${2:-}"
    exit $?
fi

if [ ! -f "$GENE_LIST" ]; then
    echo "Error: gene list not found: $GENE_LIST" >&2
    exit 1
fi

if [ ! -d "$WORK_DIR/combined_unaligned_${DATASET_MODERN}" ]; then
    echo "Error: modern FASTA directory not found." >&2
    exit 1
fi
if [ ! -d "$WORK_DIR/combined_unaligned_${DATASET_WITH_ALL}" ]; then
    echo "Error: all-sample FASTA directory not found." >&2
    exit 1
fi

mkdir -p "$WORK_DIR/final_alignments_${DATASET_MODERN}" \
    "$WORK_DIR/final_alignments_${DATASET_WITH_ALL}"

grep -v '^$' "$GENE_LIST" |
    xargs -P "$THREADS" -I {} bash "$SCRIPT_PATH" __align_gene "{}"

echo "Modern-backbone and ancient-added alignments complete."
