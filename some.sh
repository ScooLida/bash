#mileup - загонять все vcf

export PATH="$HOME/miniconda/bin:$PATH"

#проверка на наличие 
if test -f "$var.vcf"; then
    echo "$var exists."
else echo "$var not"
fi  

#показать строку
sed '2221,2222!d' SRR17908658.vcf

#архив
bgzip -c SRR10011655.vcf > SRR10011655.vcf.gz
gzip -dk SRR10011655.vcf.gz



fasterq-dump --split-files SRR17055867.sra

#геном Ochotona_princeps
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Ochotona_princeps/latest_assembly_versions/GCF_014633375.1_OchPri4.0/ ~/cow/oc 


#удаление адаптеров
export PATH="$HOME/miniconda/bin:$PATH"
conda install -c bioconda adapterremoval
AdapterRemoval --file1 1k_S1_R1_001.fastq  --file2 1k_S1_R2_001.fastq --basename output_paired --trimns --trimqualities



#индекс пищухи
bowtie2-build reference_sequence.fasta index_name

#выравнивание на геном пищухи bowtie. получение .sam
bowtie2 -p 6 -q --very-sensitive-local -x oc_index -U output_paired.pair1.truncated output_paired.pair2.truncated  -S lep.sam

# форматирование samtools
samtools view -bt oc_genome.fa.fai -o lep.bam lep.sam
samtools sort lep.bam -o lep_sort.bam
samtools index lep_sort.bam
bcftools mpileup -f oc_genome.fa lep_sort.bam | bcftools call -mv -Oz -o calls.vcf.gz

#фильтрация
bcftools filter calls.vcf.gz -s LowQual -e 'QUAL<20 && DP<10' > filtered.vcf

#статистика
samtools flagstat -@ 10 dna.bam > flagstat.txt 


#получаем формат .traw 
~/plink --vcf SRR10011655.vcf.gz --recode A-transpose --allow-extra-chr --out SRR10011655




