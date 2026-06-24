#!/bin/bash
# Builds individual gene trees with IQ-TREE and combines them into a species tree with ASTRAL.
# Takes a list of gene names, reads aligned FASTA files from ALIGN_DIR,
# runs IQ-TREE for each gene in parallel, concatenates the resulting tree files,
# and runs ASTRAL to produce a species tree.

# ==============================================================================
THREADS=9

WORK_DIR="./subset_parallel/pipeline_bulletproof_final"
ALIGN_DIR="$WORK_DIR/final_alignments"
GENE_LIST="$WORK_DIR/list_30_percent_gap.txt"

GENE_TREES_DIR="$WORK_DIR/gene_trees"
ALL_GENE_TREES_FILE="$WORK_DIR/all_gene_trees.treefile"

# Path to ASTRAL jar file (update version/path if needed)
ASTRAL_JAR="~/Astral/astral.5.16.3.jar"
ASTRAL_OUTPUT="$WORK_DIR/astral_species_tree.tre"

# Standard substitution model for all genes
MODEL="GTR+G"
# ==============================================================================

if [ ! -d "$ALIGN_DIR" ]; then
    echo "Error: alignment directory $ALIGN_DIR not found!"
    exit 1
fi

if [ ! -f "$GENE_LIST" ]; then
    echo "Error: gene list $GENE_LIST not found!"
    exit 1
fi

mkdir -p "$GENE_TREES_DIR"

total_genes=$(wc -l < "$GENE_LIST" | awk '{print $1}')
echo "Found $total_genes genes in the list"
echo "Using model: $MODEL"
echo "------------------------------------------------------------------------"

# --- STEP 1: RUN IQ-TREE FOR EACH GENE IN PARALLEL ---
echo "Step 1: Building individual gene trees using $THREADS thread(s)..."

run_iqtree_single() {
    file=$1
    gene=$(basename "$file" .fasta)

    # -st DNA avoids sequence type detection errors
    iqtree -s "$file" -st DNA -m "$MODEL" -nt 1 --prefix "$GENE_TREES_DIR/$gene" -quiet

    if [ $? -eq 0 ]; then
        echo "Tree ready for gene: $gene"
    else
        echo "Error building tree for gene: $gene"
    fi
}
export -f run_iqtree_single
export ALIGN_DIR
export GENE_TREES_DIR
export MODEL

# Process only genes present in the list and existing as FASTA files
while read -r gene; do
    [ -z "$gene" ] && continue
    file="$ALIGN_DIR/${gene}.fasta"
    [ -f "$file" ] || continue
    echo "$file"
done < "$GENE_LIST" | xargs -P "$THREADS" -I {} bash -c 'run_iqtree_single "{}"'

# --- STEP 2: CONCATENATE ALL GENE TREES INTO ONE FILE ---
echo "------------------------------------------------------------------------"
echo "Step 2: Concatenating individual gene trees..."

> "$ALL_GENE_TREES_FILE"
for treefile in "$GENE_TREES_DIR"/*.treefile; do
    [ -e "$treefile" ] || continue
    cat "$treefile" >> "$ALL_GENE_TREES_FILE"
done

total_trees=$(wc -l < "$ALL_GENE_TREES_FILE" | awk '{print $1}')
echo "Collected gene trees: $total_trees (saved to $ALL_GENE_TREES_FILE)"

# --- STEP 3: RUN ASTRAL ---
echo "------------------------------------------------------------------------"
echo "Step 3: Running ASTRAL to build species tree..."

if [ ! -f "$ASTRAL_JAR" ]; then
    echo "Error: ASTRAL jar file not found at '$ASTRAL_JAR'"
    echo "Please copy ASTRAL to the expected location or update ASTRAL_JAR."
    echo "Gene trees are available in: $GENE_TREES_DIR"
    exit 1
fi

java -jar "$ASTRAL_JAR" -i "$ALL_GENE_TREES_FILE" -o "$ASTRAL_OUTPUT" 2>"$WORK_DIR/astral.log"

if [ $? -eq 0 ] && [ -f "$ASTRAL_OUTPUT" ]; then
    echo "ASTRAL species tree built successfully."
    echo "Output saved to: $ASTRAL_OUTPUT"
else
    echo "Error running ASTRAL. Check log: $WORK_DIR/astral.log"
fi
echo "========================================================================"
