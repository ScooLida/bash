samples="1k,4k,3k"

# X-хромосома
X="NC_091453.1"

# Митохондриальная хромосома,
MT=""

# Общий файл с результатами
result_file="sex_results.txt"

# Заголовок
echo -e "Sample\tX_depth\tAutosomal_mean_depth\tX_autosome_ratio\tSex" > "$result_file"

# Разделяем список по запятым
IFS=',' read -ra sample_array <<< "$samples"

for sample in "${sample_array[@]}"
do
    echo "Processing: $sample"

    bam="${sample}.rescaled.bam"
    filtered_bam="${sample}.filtered.bam"
    coverage_file="${sample}.coverage.txt"

    # 1. Фильтрация BAM
    # MAPQ >= 20
    # убрать unmapped, secondary, QC-fail, duplicates, supplementary
#    samtools view -b \
 #       -q 20 \
  ##      -F 3844 \
    #    "$bam" \
    #    > "$filtered_bam"

    # 2. Индексация
   # samtools index "$filtered_bam"

    # 3. Coverage по хромосомам
    samtools coverage "$filtered_bam" \
        > "$coverage_file"
    # 4. X/autosome ratio и запись результата
    awk -v SAMPLE="$sample" \
        -v X="$X" \
        -v MT="$MT" '
    $1 == X {
        x=$7
    }

    $1 ~ /^NC_/ && $1 != X && (MT=="" || $1 != MT) {
        sum += $7
        n++
    }

    END {
        auto=sum/n
        ratio=x/auto

        if (ratio < 0.65)
            sex="MALE"
        else if (ratio > 0.80)
            sex="FEMALE"
        else
            sex="AMBIGUOUS"

        print SAMPLE "\t" x "\t" auto "\t" ratio "\t" sex
    }' "$coverage_file" >> "$result_file"

done

echo "Done. Results:"
column -t "$result_file"
C


