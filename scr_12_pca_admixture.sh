#!/bin/bash
# Run PCA and ADMIXTURE twice: without 1k, 3k, 5k, and 5kS8, then with all samples.
# ADMIXTURE is tested across a K range and the K with the lowest CV error is reported.

set -euo pipefail

# Configuration
DATA_PREFIX="MyHare"
WITHOUT_VCF="${DATA_PREFIX}_without_1k_3k_5k_5kS8.vcf.gz"
ALL_VCF="${DATA_PREFIX}_with_all_samples.vcf.gz"
ANALYSIS_DIR="./population_analysis"
PLINK="$HOME/plink"
ADMIXTURE="$HOME/admixture/dist/admixture_linux-1.3.0/admixture"
PLINK_THREADS=8
ADMIXTURE_THREADS=4
MISSINGNESS=0.2
PCA_COMPONENTS=4
K_MIN=2
K_MAX=15

run_analysis() {
    local label=$1
    local input_vcf=$2
    local prefix="$ANALYSIS_DIR/$label"
    local cv_table="${prefix}_admixture_cv.tsv"

    if [ ! -f "$input_vcf" ]; then
        echo "Error: VCF not found for $label: $input_vcf"
        return 1
    fi

    "$PLINK" --vcf "$input_vcf" --geno "$MISSINGNESS" --make-bed \
        --threads "$PLINK_THREADS" --out "$prefix" --allow-extra-chr
    "$PLINK" --bfile "$prefix" --pca "$PCA_COMPONENTS" \
        --threads "$PLINK_THREADS" --allow-extra-chr --out "$prefix"

    printf "K\tCV_error\n" > "$cv_table"
    for K in $(seq "$K_MIN" "$K_MAX"); do
        log_file="${prefix}_admixture_K${K}.log"
        "$ADMIXTURE" -j"$ADMIXTURE_THREADS" --cv "$prefix.bed" "$K" > "$log_file" 2>&1
        cv_error=$(awk '/CV error/ { value=$NF } END { print value }' "$log_file")
        if [ -n "$cv_error" ]; then
            printf "%s\t%s\n" "$K" "$cv_error" >> "$cv_table"
        fi
    done

    best_k=$(awk 'NR > 1 && $2 != "" { print }' "$cv_table" | sort -k2,2n | head -n 1 | cut -f1)
    if [ -z "$best_k" ]; then
        echo "Error: could not determine optimal K for $label."
        return 1
    fi
    echo "Optimal K for $label: $best_k"
    echo "$best_k" > "${prefix}_optimal_K.txt"
}

mkdir -p "$ANALYSIS_DIR"
run_analysis "without_1k_3k_5k_5kS8" "$WITHOUT_VCF"
run_analysis "with_all_samples" "$ALL_VCF"

echo "PCA and ADMIXTURE analyses complete."
