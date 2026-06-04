#!/bin/bash

# 1. Основные настройки
SAMPLES="1k,5kS8,SRR11020265,SRR14535648,SRR32541919,3k,SRR17908655,SRR5949633,4k,SRR10011656,SRR14535637,SRR14535689,SRR28956761"
echo "$SAMPLES" | tr ',' '\n' | grep "." > sample_list.txt

THREADS=4
REFERENCE="$HOME/hare_work/krol_g.fasta"

# Укажите путь к вашему файлу со списком хромосом (если он в этой же папке, оставьте так)
CHROM_LIST="chr.txt"

echo "------------------------------------------------"
echo "Этап 0: Проверка исходных файлов..."

# 0.1 Проверка BAM-файлов
MISSING_FILES=0
for SAMPLE in $(cat sample_list.txt); do
  FILE_PATH="$HOME/hare_work/MyKrol2/my_genome/${SAMPLE}/${SAMPLE}.rescaled.bam"
  if [ ! -f "$FILE_PATH" ]; then
    echo -e "\e[31m[ОШИБКА]\e[0m BAM-файл не найден: $FILE_PATH"
    MISSING_FILES=$((MISSING_FILES + 1))
  fi
done

if [ "$MISSING_FILES" -gt 0 ]; then
  echo "Обнаружено отсутствующих файлов: $MISSING_FILES. Скрипт остановлен."
  exit 1
fi

# 0.2 Проверка файла с хромосомами
if [ ! -f "$CHROM_LIST" ]; then
    echo -e "\e[31m[ОШИБКА]\e[0m Файл со списком хромосом ($CHROM_LIST) не найден!"
    exit 1
else
    echo -e "\e[32m[OK]\e[0m Все файлы на месте. В списке хромосом строк: $(wc -l < $CHROM_LIST)"
fi
echo "------------------------------------------------"

echo "Этап 1: Запуск параллельного вызова вариантов (только по chr.txt)..."
# В mpileup добавлен флаг -R $CHROM_LIST
cat sample_list.txt | parallel --tmpdir . --bar -j $THREADS "
  if [ -f '{}.vcf.gz.tbi' ]; then
    echo '[ПРОПУСК] Образец {} уже готов.'
  else
    bcftools mpileup -a DP -d 50 -R $CHROM_LIST -f $REFERENCE $HOME/hare_work/MyKrol2/my_genome/{}/{}.rescaled.bam | \
    bcftools call -mv -Oz -o {}.vcf.gz && \
    tabix -p vcf {}.vcf.gz
  fi
"

echo "Этап 2: Объединение VCF файлов в один..."
ls *.vcf.gz | grep -v "merged" | grep -v "filtered" > vcf_to_merge.txt
bcftools merge -l vcf_to_merge.txt -Oz -o MyKrol2_selected_merged.vcf.gz
tabix -p vcf MyKrol2_selected_merged.vcf.gz

echo "Этап 3: Фильтрация генотипов (Качество < 20 ИЛИ Глубина < 3)..."
bcftools view MyKrol2_selected_merged.vcf.gz -e 'QUAL < 20 || FMT/DP < 3' -O z -o All_filtered_DP3_Q20.vcf.gz
tabix -p vcf All_filtered_DP3_Q20.vcf.gz

echo "------------------------------------------------"
echo "Успешно завершено! Итоговый файл: All_filtered_DP3_Q20.vcf.gz"
