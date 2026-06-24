#!/bin/bash
# scr_04_bam_to_busco.sh
# Pipeline: extract gene loci from BAM files using a BUSCO BED file,
# merge per-sample FASTAs into one file per gene, and filter the merged FASTAs.
# Genes shorter than MIN_SEQ_LENGTH are discarded after extraction.

# ==============================================================================
# --- CONFIG (all variables in one place) --------------------------------------
# ==============================================================================
SAMPLES="SRR14535670,SRR32541919"
REF_GENOME="../krol_g.fasta"
BED_FILE="../krol_genes_fixed.bed"  # (из busco)
BAM_DIR="../MyKrol2/my_genome"
OLD_GENES_DIR="./subset_parallel/genes_by_locus"

MIN_COVERAGE=2
NUM_THREADS=4          # Optimal: 4 threads + SLEEP_INTERVAL protect the disk from overload
SLEEP_INTERVAL=0.1
MIN_SEQ_LENGTH=10      # Genes shorter than this are discarded after extraction
MAX_SEQ_LENGTH=5000    # Genes longer than this are discarded during filtering
MIN_ATGC_PERCENT=30    # Minimum ATGC percentage for a gene to pass filtering
REQUIRED_SAMPLES="1k,4k,5kS8,3k"
PROGRESS_INTERVAL=500  # Print progress every N files during filtering

WORK_DIR="./subset_parallel/pipeline_bulletproof_final"
EXTRACTED_DIR="$WORK_DIR/1_extracted_fasta"
COMBINED_DIR="$WORK_DIR/2_combined_unaligned"
EXTRACTED_GENES_LIST="$WORK_DIR/extracted_genes_list.txt"
GOLD_LIST="$WORK_DIR/gold_genes_list.txt"
# ==============================================================================

mkdir -p "$WORK_DIR" "$EXTRACTED_DIR" "$COMBINED_DIR"

# Helper: format seconds as human-readable time
display_time() {
    local T=$1
    local H=$((T/3600))
    local M=$((T%3600/60))
    local S=$((T%60))
    if [ $H -gt 0 ]; then printf "%dч %dм %dс" $H $M $S
    elif [ $M -gt 0 ]; then printf "%dм %dс" $M $S
    else printf "%dс" $S; fi
}

echo "========================================================================"
echo " BAM -> BUSCO FASTA extraction, merge and filter"
echo "========================================================================"

> "$EXTRACTED_GENES_LIST"
> "$GOLD_LIST"

if [ ! -f "${REF_GENOME}.fai" ]; then samtools faidx "$REF_GENOME"; fi

SAMPLE_ARRAY=($(echo "$SAMPLES" | tr ',' ' '))
TOTAL_GENES=$(wc -l < "$BED_FILE" | awk '{print $1}')

# ==============================================================================
# EXTRACTION (parallel)
# ==============================================================================
echo -e "\n=== Extracting loci from BAM ==="

for SAMPLE in "${SAMPLE_ARRAY[@]}"; do
    BAM_FILE="$BAM_DIR/${SAMPLE}/${SAMPLE}.rescaled.bam"
    SUB_EXTRACT_DIR="$EXTRACTED_DIR/$SAMPLE"
    mkdir -p "$SUB_EXTRACT_DIR"

    start_extract=$(date +%s)
    job_count=0

    # Read the whole BED file
    while read -r chrom start end gene; do
        [ -z "$gene" ] && continue

        chrom=$(echo "$chrom" | tr -d '\r' | tr -d ' ')
        start=$(echo "$start" | tr -d '\r' | tr -d ' ')
        end=$(echo "$end" | tr -d '\r' | tr -d ' ')
        region="${chrom}:${start}-${end}"

        (
            TEMP_VCF="$SUB_EXTRACT_DIR/${gene}_temp.vcf.gz"
            MASK_BED="$SUB_EXTRACT_DIR/${gene}_mask.bed"
            NEW_FASTA="$SUB_EXTRACT_DIR/${gene}.fasta"
            TEMP_REF="$SUB_EXTRACT_DIR/${gene}_ref.fasta"

            samtools faidx "$REF_GENOME" "$region" > "$TEMP_REF" 2>/dev/null
            bcftools mpileup -Ou -f "$REF_GENOME" -r "$region" "$BAM_FILE" 2>/dev/null | \
                bcftools call -c -Oz -o "$TEMP_VCF" 2>/dev/null
            bcftools index -t "$TEMP_VCF" 2>/dev/null
            samtools depth -a -r "$region" "$BAM_FILE" 2>/dev/null | \
                awk -v min_cov="$MIN_COVERAGE" '$3 < min_cov {print $1"\t"$2-1"\t"$2}' > "$MASK_BED"

            echo ">${SAMPLE}" > "$NEW_FASTA"
            bcftools consensus -m "$MASK_BED" -f "$TEMP_REF" "$TEMP_VCF" 2>/dev/null | grep -v "^>" >> "$NEW_FASTA"

            rm -f "$TEMP_VCF" "${TEMP_VCF}.tbi" "$MASK_BED" "$TEMP_REF"
        ) </dev/null &

        sleep "$SLEEP_INTERVAL"  # Protect the disk from choking
        ((job_count++))

        if (( job_count % NUM_THREADS == 0 )); then
            wait
            curr_time=$(date +%s)
            elapsed=$((curr_time - start_extract))
            if [ $elapsed -eq 0 ]; then elapsed=1; fi
            remaining_genes=$((TOTAL_GENES - job_count))
            eta_sec=$(( (elapsed * remaining_genes) / job_count ))
            echo "Progress: $job_count / $TOTAL_GENES | Elapsed: $(display_time $elapsed) | ETA: $(display_time $eta_sec)"
        fi
    done < "$BED_FILE"
    wait
done

for s in "${SAMPLE_ARRAY[@]}"; do
    for f in "$EXTRACTED_DIR/$s"/*.fasta; do
        [ -e "$f" ] || continue
        len=$(grep -v "^>" "$f" | tr -d '\n' | wc -c)
        if [ "$len" -gt "$MIN_SEQ_LENGTH" ]; then basename "$f" .fasta >> "$EXTRACTED_GENES_LIST"
        else rm -f "$f"; fi
    done
done
sort -u "$EXTRACTED_GENES_LIST" -o "$EXTRACTED_GENES_LIST"
echo "Extraction complete. Total genes extracted: $(wc -l < "$EXTRACTED_GENES_LIST")"

# ==============================================================================
# MERGE
# ==============================================================================
echo -e "\n=== Merging samples ==="
combined_count=0

while read -r gene; do
    [ -z "$gene" ] && continue
    old_fasta="$OLD_GENES_DIR/${gene}.fasta"
    combined_file="$COMBINED_DIR/${gene}.fasta"

    > "$combined_file"
    if [ -f "$old_fasta" ]; then cat "$old_fasta" >> "$combined_file"; fi

    for s in "${SAMPLE_ARRAY[@]}"; do
        if [ -f "$EXTRACTED_DIR/$s/${gene}.fasta" ]; then
            cat "$EXTRACTED_DIR/$s/${gene}.fasta" >> "$combined_file"
        fi
    done
    ((combined_count++))
done < "$EXTRACTED_GENES_LIST"
echo "Merge complete. Merged or created files: $combined_count"

# ==============================================================================
# FILTER (length first, then quality / required samples)
# ==============================================================================
echo -e "\n=== Filtering quality ==="

if [ ! -d "$COMBINED_DIR" ]; then
    echo "Error: directory $COMBINED_DIR not found!"
    exit 1
fi

total_files=$(ls -1 "$COMBINED_DIR"/*.fasta 2>/dev/null | wc -l)
if [ "$total_files" -eq 0 ]; then
    echo "No files to check in $COMBINED_DIR!"
    exit 1
fi

echo "Files to check: $total_files"
echo -e "Starting check...\n"

passed_length=0
passed_vip=0
current=0

for file in "$COMBINED_DIR"/*.fasta; do
    [ -e "$file" ] || continue
    gene=$(basename "$file" .fasta)

    res=$(awk -v max_len="$MAX_SEQ_LENGTH" -v min_atgc="$MIN_ATGC_PERCENT" -v req_samples="$REQUIRED_SAMPLES" '
    BEGIN {
        split(req_samples, req, ",");
        for(i in req) required[req[i]]=1;
    }
    /^>/ {
        if (id != "") seqs[id] = seq;
        id = substr($1, 2); seq = ""; next;
    }
    { seq = seq $0; }
    END {
        if (id != "") seqs[id] = seq;

        # --- STEP 1: FAST LENGTH CHECK ---
        for (sid in seqs) {
            if (length(seqs[sid]) > max_len) {
                print "FAIL_LENGTH";
                exit;
            }
        }

        # --- STEP 2: REQUIRED SAMPLES AND ATGC QUALITY CHECK ---
        all_len = 0; all_let = 0; missing_or_bad_req = 0;

        for (sid in seqs) {
            s = seqs[sid]; len = length(s); all_len += len;
            gsub(/[^ATGCatgc]/, "", s); let_cnt = length(s); all_let += let_cnt;

            if (sid in required) {
                required_found[sid] = 1;
                if (len == 0 || (let_cnt * 100 / len) < min_atgc) missing_or_bad_req = 1;
            }
        }

        for (sid in required) { if (!required_found[sid]) missing_or_bad_req = 1; }

        if (all_len > 0 && (all_let * 100 / all_len) >= min_atgc && missing_or_bad_req == 0) {
            print "PASS";
        } else {
            print "FAIL_VIP";
        }
    }' "$file")

    if [ "$res" == "PASS" ]; then
        ((passed_length++))
        ((passed_vip++))
        echo "$gene" >> "$GOLD_LIST"
    elif [ "$res" == "FAIL_VIP" ]; then
        ((passed_length++))
    fi

    ((current++))
    if (( current % PROGRESS_INTERVAL == 0 )); then
        echo "   ... checked $current / $total_files files"
    fi
done

echo -e "\n========================================================================"
echo "FILTERING STATISTICS:"
echo "========================================================================"
echo "Input genes:                $total_files"
echo "Passed length filter:       $passed_length (removed $((total_files - passed_length)) giants)"
echo "Passed VIP/quality filter:  $passed_vip"
echo "========================================================================"
echo "Result saved to: $GOLD_LIST"
