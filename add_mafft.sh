#!/bin/bash

# ==============================================================================
# --- БЛОК НАСТРОЕК (CONFIG) ---------------------------------------------------
# ==============================================================================
SAMPLE_NAME="SRR14535670"
BAM_FILE="../MyKrol2/my_genome/SRR14535670/SRR14535670.rescaled.bam"
REF_GENOME="../krol_g.fasta"
BED_FILE="../krol_genes_fixed.bed"

ALIGNED_DIR="./subset_parallel/aligned_genes"

MIN_TAXA=4              
MIN_LETTERS_PCT=30      

WORK_DIR="./subset_parallel/pipeline_${SAMPLE_NAME}"
EXTRACTED_DIR="$WORK_DIR/extracted_fasta"
FINAL_ALIGNED_DIR="$WORK_DIR/final_aligned_14_species"
INTERNAL_GENE_LIST="$WORK_DIR/internal_target_genes.txt" 
# ==============================================================================

echo "========================================================================"
echo " ЗАПУСК ПАЙПЛАЙНА (Фильтр: >= $MIN_LETTERS_PCT% букв и >= $MIN_TAXA видов)"
echo "========================================================================"

mkdir -p "$EXTRACTED_DIR"
mkdir -p "$FINAL_ALIGNED_DIR"
> "$INTERNAL_GENE_LIST"

if [ ! -f "${REF_GENOME}.fai" ]; then
    echo "⏳ Индекс для референса не найден. Создаю (samtools faidx)..."
    samtools faidx "$REF_GENOME"
fi

echo ""
echo "=== ЭТАП 1: Строгая проверка качества выравниваний в $ALIGNED_DIR ==="
count_target=0

for aln_file in "$ALIGNED_DIR"/*.fasta; do
    [ -e "$aln_file" ] || continue
    
    taxa_count=$(grep -c "^>" "$aln_file" || echo "0")
    total_chars=$(grep -v "^>" "$aln_file" 2>/dev/null | tr -d '\n' | wc -c)
    letters=$(grep -v "^>" "$aln_file" 2>/dev/null | tr -d '\n-' | wc -c)
    
    if [ "$total_chars" -eq 0 ]; then continue; fi
    pct_letters=$(( letters * 100 / total_chars ))
    
    if [ "$taxa_count" -ge "$MIN_TAXA" ] && [ "$pct_letters" -ge "$MIN_LETTERS_PCT" ]; then
        gene_name=$(basename "$aln_file" .fasta)
        echo "$gene_name" >> "$INTERNAL_GENE_LIST"
        count_target=$((count_target + 1))
    fi
done

echo "Найдено идеальных генов-кандидатов (>= 30% букв): $count_target"
echo ""

if [ "$count_target" -eq 0 ]; then
    echo "❌ Ошибка: Не найдено подходящих генов!"
    exit 1
fi

echo "=== ЭТАП 2: ИЗВЛЕЧЕНИЕ ИЗ BAM И ПРОФИЛЬНОЕ ВЫРАВНИВАНИЕ ==="
extracted=0
aligned_success=0

while read -r gene; do
    [ -z "$gene" ] && continue
    
    coord_line=$(grep -w "$gene" "$BED_FILE" | head -n 1)
    if [ -z "$coord_line" ]; then
        continue
    fi
    
    chrom=$(echo "$coord_line" | awk '{print $1}')
    start=$(echo "$coord_line" | awk '{print $2}')
    end=$(echo "$coord_line" | awk '{print $3}')
    region="${chrom}:${start}-${end}"
    
    TEMP_VCF="$EXTRACTED_DIR/${gene}_temp.vcf.gz"
    NEW_FASTA="$EXTRACTED_DIR/${gene}.fasta"
    OLD_FASTA="$ALIGNED_DIR/${gene}.fasta"
    FINAL_FASTA="$FINAL_ALIGNED_DIR/${gene}.fasta"
    
    # 1. Извлечение
    bcftools mpileup -Ou -f "$REF_GENOME" -r "$region" "$BAM_FILE" 2>/dev/null | \
    bcftools call -c -Oz -o "$TEMP_VCF" 2>/dev/null
    
    bcftools index "$TEMP_VCF"
    echo ">${SAMPLE_NAME}" > "$NEW_FASTA"
    bcftools consensus -f "$REF_GENOME" "$TEMP_VCF" 2>/dev/null | grep -v "^>" >> "$NEW_FASTA"
    rm -f "$TEMP_VCF" "${TEMP_VCF}.csi"
    
    # 2. ПРОВЕРКА НА ПУСТОТУ (наш фикс!)
    # Считаем количество реальных букв в полученном файле
    seq_len=$(grep -v "^>" "$NEW_FASTA" 2>/dev/null | tr -d '\n' | wc -c)
    
    if [ "$seq_len" -lt 10 ]; then
        echo "⚠️ Пропуск $gene: в BAM нет данных для региона $region"
        continue
    fi
    
    extracted=$((extracted + 1))
    
    # 3. Выравнивание
    if [ -f "$OLD_FASTA" ]; then
        mafft --quiet --add "$NEW_FASTA" "$OLD_FASTA" > "$FINAL_FASTA" 2>/dev/null
        aligned_success=$((aligned_success + 1))
        
        if (( aligned_success % 10 == 0 )); then
            echo "Обработано $aligned_success из $count_target генов..."
        fi
    fi

done < "$INTERNAL_GENE_LIST"

echo "========================================================================"
echo "✅ ГОТОВО! Успешно создано выравниваний: $aligned_success"
echo "========================================================================"
