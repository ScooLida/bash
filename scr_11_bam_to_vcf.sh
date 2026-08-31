#!/bin/bash
# Call and filter variants from Paleomix BAM files.
# Per-sample VCFs are generated once. Modern samples define the filtered variant
# sites; the all-sample VCF keeps exactly those sites and retains missing ancient
# genotypes instead of using ancient coverage to remove a site.

set -euo pipefail

# Configuration
THREADS=3
REFERENCE="$HOME/hare_work/krol_g.fasta"
BAM_DIR="$HOME/hare_work/My_new_krol/my_genome"
SAMPLE_LIST="sample_list.txt"
CHROM_LIST="chr.txt"
VCF_DIR="./per_sample_vcfs"
MERGE_LIST="./vcf_to_merge.txt"
DATA_PREFIX="MyHare"
SELECTED_SITES="./${DATA_PREFIX}_modern_selected_sites.tsv"
MODERN_MERGE_LIST="./vcf_to_merge_modern.txt"
QUAL=20
DEPTH=3
ANCIENT_SAMPLES=("1k" "3k" "4k" "5kS8")
SCRIPT_PATH="$(readlink -f "$0")"

call_sample() {
    local sample="${1:-}"
    if [ -z "$sample" ]; then
        echo "Error: empty sample identifier." >&2
        return 1
    fi

    local bam="$BAM_DIR/$sample/$sample.rescaled.bam"
    local output="$VCF_DIR/$sample.vcf.gz"

    if [ -f "$output" ] && [ -f "$output.tbi" ]; then
        echo "VCF already exists: $sample"
        return
    fi

    if [ ! -f "$bam.bai" ] && [ ! -f "${bam%.bam}.bai" ] && [ ! -f "$bam.csi" ]; then
        samtools index "$bam"
    fi

    echo "Starting variant calling: $sample"
    if ! bcftools mpileup -a DP -d 250 -R "$CHROM_LIST" -f "$REFERENCE" "$bam" |
        bcftools call -mv -Oz -o "$output"; then
        echo "Error: variant calling failed for sample: $sample" >&2
        return 1
    fi
    tabix -p vcf "$output"
    echo "Completed variant calling: $sample"
}

if [ "${1:-}" = "__call_sample" ]; then
    call_sample "${2:-}"
    exit $?
fi

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

mapfile -t SAMPLES < <(
    tr -d '\r' < "$SAMPLE_LIST" |
        awk 'NF { print $1 }'
)
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

printf '%s\n' "${SAMPLES[@]}" |
    parallel --tmpdir . --bar -j "$THREADS" \
        bash "$SCRIPT_PATH" __call_sample {}

is_ancient() {
    local candidate=$1
    local ancient
    for ancient in "${ANCIENT_SAMPLES[@]}"; do
        [ "$candidate" = "$ancient" ] && return 0
    done
    return 1
}

: > "$MERGE_LIST"
: > "$MODERN_MERGE_LIST"
for sample in "${SAMPLES[@]}"; do
    printf '%s\n' "$VCF_DIR/$sample.vcf.gz" >> "$MERGE_LIST"
    if ! is_ancient "$sample"; then
        printf '%s\n' "$VCF_DIR/$sample.vcf.gz" >> "$MODERN_MERGE_LIST"
    fi
done

MERGED_VCF="${DATA_PREFIX}_merged_all.vcf.gz"
bcftools merge -l "$MERGE_LIST" -Oz -o "$MERGED_VCF"
tabix -p vcf "$MERGED_VCF"

ALL_SAMPLE_FILE="${DATA_PREFIX}_samples_all.txt"
MODERN_SAMPLE_FILE="${DATA_PREFIX}_samples_modern.txt"
printf '%s\n' "${SAMPLES[@]}" > "$ALL_SAMPLE_FILE"
{
    for sample in "${SAMPLES[@]}"; do
        if ! is_ancient "$sample"; then
            printf '%s\n' "$sample"
        fi
    done
} > "$MODERN_SAMPLE_FILE"

filter_modern_dataset() {
    if [ ! -s "$MODERN_SAMPLE_FILE" ]; then
        echo "Error: no modern samples remain." >&2
        return 1
    fi

    MERGED_MODERN_VCF="${DATA_PREFIX}_merged_modern.vcf.gz"
    bcftools merge -l "$MODERN_MERGE_LIST" -Oz -o "$MERGED_MODERN_VCF"
    tabix -p vcf "$MERGED_MODERN_VCF"

    bcftools view \
        -e "QUAL < ${QUAL} || MIN(FMT/DP) < ${DEPTH}" \
        -Oz -o "${DATA_PREFIX}_modern.vcf.gz" "$MERGED_MODERN_VCF"
    tabix -p vcf "${DATA_PREFIX}_modern.vcf.gz"
    echo "Created: ${DATA_PREFIX}_modern.vcf.gz"
}

filter_modern_dataset

# Use modern-selected variant positions as the master list for the all-sample
# dataset. Merge the filtered modern VCF with ancient VCFs at those positions;
# missing ancient genotypes remain ./., rather than removing the site.
bcftools query -f '%CHROM\t%POS\t%POS\n' \
    "${DATA_PREFIX}_modern.vcf.gz" > "$SELECTED_SITES"
if [ ! -s "$SELECTED_SITES" ]; then
    echo "Error: no modern variants passed filtering." >&2
    exit 1
fi

ALL_MERGE_INPUTS=("${DATA_PREFIX}_modern.vcf.gz")
for sample in "${SAMPLES[@]}"; do
    if is_ancient "$sample"; then
        ALL_MERGE_INPUTS+=("$VCF_DIR/$sample.vcf.gz")
    fi
done

bcftools merge -R "$SELECTED_SITES" -Oz \
    -o "${DATA_PREFIX}_with_all_samples.vcf.gz" \
    "${ALL_MERGE_INPUTS[@]}"
tabix -p vcf "${DATA_PREFIX}_with_all_samples.vcf.gz"
echo "Created: ${DATA_PREFIX}_with_all_samples.vcf.gz"

echo "Variant calling and filtering complete."
