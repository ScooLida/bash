#!/bin/bash
# Call and filter variants from Paleomix BAM files.
# Per-sample VCFs are generated once, then two filtered datasets are created:
# one without 1k, 3k, 5k, and 5kS8, and one with all samples.

set -euo pipefail

# Configuration
THREADS=7
REFERENCE="$HOME/hare_work/krol_g.fasta"
BAM_DIR="$HOME/hare_work/MyKrol2/my_genome"
SAMPLE_LIST="sample_list.txt"
CHROM_LIST="chr.txt"
VCF_DIR="./per_sample_vcfs"
MERGE_LIST="./vcf_to_merge.txt"
DATA_PREFIX="MyHare"
QUAL=20
DEPTH=3
EXCLUDED_SAMPLES=("1k" "3k" "5k" "5kS8")

if [ ! -f "$SAMPLE_LIST" ]; then
    echo "Error: sample list not found: $SAMPLE_LIST"
    exit 1
fi
if [ ! -f "$CHROM_LIST" ]; then
    echo "Error: chromosome list not found: $CHROM_LIST"
    exit 1
fi
if [ ! -f "$REFERENCE" ]; then
    echo "Error: reference genome not found: $REFERENCE"
    exit 1
fi

mapfile -t SAMPLES < <(awk 'NF { print $1 }' "$SAMPLE_LIST")
if [ "${#SAMPLES[@]}" -eq 0 ]; then
    echo "Error: sample list is empty."
    exit 1
fi

mkdir -p "$VCF_DIR"

for sample in "${SAMPLES[@]}"; do
    bam="$BAM_DIR/$sample/$sample.rescaled.bam"
    if [ ! -f "$bam" ]; then
        echo "Error: BAM file not found for $sample: $bam"
        exit 1
    fi
done

call_sample() {
    local sample=$1
    local bam="$BAM_DIR/$sample/$sample.rescaled.bam"
    local output="$VCF_DIR/$sample.vcf.gz"

    if [ -f "$output" ] && [ -f "$output.tbi" ]; then
        echo "VCF already exists: $sample"
        return
    fi

    if [ ! -f "$bam.bai" ] && [ ! -f "${bam%.bam}.bai" ] && [ ! -f "$bam.csi" ]; then
        samtools index "$bam"
    fi

    bcftools mpileup -a DP -d 250 -R "$CHROM_LIST" -f "$REFERENCE" "$bam" |
        bcftools call -mv -Oz -o "$output"
    tabix -p vcf "$output"
}
export -f call_sample
export BAM_DIR CHROM_LIST REFERENCE VCF_DIR

printf '%s\n' "${SAMPLES[@]}" |
    parallel --tmpdir . --bar -j "$THREADS" bash -c 'call_sample "$@"' _ {}

printf '%s\n' "${SAMPLES[@]/#/$VCF_DIR/}" | sed 's/$/.vcf.gz/' > "$MERGE_LIST"

MERGED_VCF="${DATA_PREFIX}_merged_all.vcf.gz"
bcftools merge -l "$MERGE_LIST" -Oz -o "$MERGED_VCF"
tabix -p vcf "$MERGED_VCF"

is_excluded() {
    local candidate=$1
    local excluded
    for excluded in "${EXCLUDED_SAMPLES[@]}"; do
        [ "$candidate" = "$excluded" ] && return 0
    done
    return 1
}

ALL_SAMPLE_FILE="${DATA_PREFIX}_samples_all.txt"
WITHOUT_SAMPLE_FILE="${DATA_PREFIX}_samples_without_1k_3k_5k_5kS8.txt"
printf '%s\n' "${SAMPLES[@]}" > "$ALL_SAMPLE_FILE"
{
    for sample in "${SAMPLES[@]}"; do
        if ! is_excluded "$sample"; then
            printf '%s\n' "$sample"
        fi
    done
} > "$WITHOUT_SAMPLE_FILE"

filter_dataset() {
    local label=$1
    local sample_file=$2
    local output="${DATA_PREFIX}_${label}.vcf.gz"

    if [ ! -s "$sample_file" ]; then
        echo "Error: no samples remain for dataset $label."
        return 1
    fi

    bcftools view -S "$sample_file" \
        -e "QUAL < ${QUAL} || MIN(FMT/DP) < ${DEPTH}" \
        -Oz -o "$output" "$MERGED_VCF"
    tabix -p vcf "$output"
    echo "Created: $output"
}

filter_dataset "without_1k_3k_5k_5kS8" "$WITHOUT_SAMPLE_FILE"
filter_dataset "with_all_samples" "$ALL_SAMPLE_FILE"

echo "Variant calling and filtering complete."
