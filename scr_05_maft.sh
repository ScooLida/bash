#!/bin/bash
# scr_05_maft.sh
# Aligns selected gene FASTA files using MAFFT in parallel.
# Takes a list of gene names, reads corresponding FASTA files from SOURCE_DIR,
# and writes aligned FASTA files to OUTPUT_DIR.

# ==============================================================================
# --- CONFIG (all variables in one place) --------------------------------------
# ==============================================================================
THREADS=9

WORK_DIR="./subset_parallel/pipeline_bulletproof_final"
GENE_LIST="$WORK_DIR/gold_genes_list.txt"
SOURCE_DIR="$WORK_DIR/combined_unaligned"
OUTPUT_DIR="$WORK_DIR/final_alignments"
# ==============================================================================

echo "========================================================================"
echo " MAFFT parallel alignment"
echo "========================================================================"

if [ ! -f "$GENE_LIST" ]; then
    echo "Error: gene list $GENE_LIST not found!"
    exit 1
fi

if [ ! -d "$SOURCE_DIR" ]; then
    echo "Error: source directory $SOURCE_DIR not found!"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

total_genes=$(wc -l < "$GENE_LIST" | awk '{print $1}')
echo "Genes to align: $total_genes"
echo "Running in $THREADS thread(s)..."
echo "------------------------------------------------------------------------"

export SOURCE_DIR
export OUTPUT_DIR

run_mafft() {
    gene=$1
    gene=$(echo "$gene" | tr -d '\r' | tr -d ' ')

    infile="$SOURCE_DIR/${gene}.fasta"
    outfile="$OUTPUT_DIR/${gene}.fasta"

    if [ -f "$infile" ]; then
        mafft --auto "$infile" > "$outfile" 2>/dev/null
        echo "Aligned: $gene"
    else
        echo "File not found: $infile"
    fi
}
export -f run_mafft

cat "$GENE_LIST" | grep -v '^$' | xargs -P "$THREADS" -I {} bash -c 'run_mafft "{}"'

echo "------------------------------------------------------------------------"
echo "Alignment complete."
echo "Output directory: $OUTPUT_DIR"
echo "========================================================================"
