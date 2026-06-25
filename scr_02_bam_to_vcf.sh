#!/bin/bash
# The script takes sample names from sample_list.txt and a list of chromosomes
# from chr.txt. It extracts BAM data produced by paleomix, keeps only reads
# mapped to the specified chromosomes, calls variants for each sample in parallel,
# merges all per-sample VCFs into one file, and finally filters the merged VCF
# by quality and read depth.

# 1. Your list of samples must be provided in sample_list.txt
#SAMPLES="1k,3k,4k,5kS8,SRR11020265,SRR14535615,SRR14535670,SRR14535681,SRR14535692,SRR17908653,SRR17908654,SRR17908655,SRR17908657,SRR5949621,SRR5949622,SRR5949623,SRR5949627"
#echo "$SAMPLES" | tr ',' '\n' | grep "." > sample_list.txt

# 2. Number of threads
THREADS=7

# Final output file prefix and filtering thresholds
DATA="MyHare"
QUAL=20
DEPTH=3
REFERENCE="/hare_work/krol_g.fasta"

# Path to your chromosome list file (leave as is if it is in the same folder)
CHROM_LIST="chr.txt"

echo "------------------------------------------------"
echo "Step 0: Checking input files..."

# 0.1 Check BAM files
MISSING_FILES=0
for SAMPLE in $(cat sample_list.txt); do
  FILE_PATH="$HOME/hare_work/MyKrol2/my_genome/${SAMPLE}/${SAMPLE}.rescaled.bam"
  if [ ! -f "$FILE_PATH" ]; then
    echo -e "\e[31m[ERROR]\e[0m BAM file not found: $FILE_PATH"
    MISSING_FILES=$((MISSING_FILES + 1))
  fi
done

if [ "$MISSING_FILES" -gt 0 ]; then
  echo "Missing files detected: $MISSING_FILES. Script stopped."
  exit 1
fi

# 0.2 Check chromosome list file
if [ ! -f "$CHROM_LIST" ]; then
    echo -e "\e[31m[ERROR]\e[0m Chromosome list file ($CHROM_LIST) not found!"
    exit 1
else
    echo -e "\e[32m[OK]\e[0m All files are in place. Chromosome list lines: $(wc -l < $CHROM_LIST)"
fi
echo "------------------------------------------------"

echo "Step 1: Running parallel variant calling (only for chromosomes in chr.txt)..."
# The -R $CHROM_LIST flag is added to mpileup
cat sample_list.txt | parallel --tmpdir . --bar -j $THREADS "
  if [ -f '{}.vcf.gz.tbi' ]; then
    echo '[SKIP] Sample {} is already ready.'
  else
    # Проверяем наличие индекса .bam.bai, если его нет — создаем
    if [ ! -f $HOME/hare_work/MyKrol2/my_genome/{}/{}.rescaled.bam.bai ]; then
      echo '[INDEX] Creating index for {}...'
      samtools index $HOME/hare_work/MyKrol2/my_genome/{}/{}.rescaled.bam
    fi

    # Запуск основного пайплайна
    bcftools mpileup -a DP -d 250 -R $CHROM_LIST -f $REFERENCE $HOME/hare_work/MyKrol2/my_genome/{}/{}.rescaled.bam | \
    bcftools call -mv -Oz -o {}.vcf.gz && \
    tabix -p vcf {}.vcf.gz
  fi
"
echo "Done, looking for .vcf.gz files"
# 2. Build a strict merge list (ONLY our required samples)
for SAMPLE in $(cat sample_list.txt); do
  echo "${SAMPLE}.vcf.gz" >> vcf_to_merge.txt
done

# 3. Merge files
echo "Merging..."
bcftools merge -l vcf_to_merge.txt -Oz -o ${DATA}_merged.vcf.gz
tabix -p vcf ${DATA}_merged.vcf.gz

# 4. Filter (Quality >= 20 AND Depth >= 3)
echo "Filtering..."
bcftools view ${DATA}_merged.vcf.gz -e "QUAL < ${QUAL} || FMT/DP < ${DEPTH}" -O z -o ${DATA}.vcf.gz
tabix -p vcf ${DATA}.vcf.gz

echo "Successfully completed! Final file: ${DATA}.vcf.gz"
