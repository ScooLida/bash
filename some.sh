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
tar xvzf archive.tar.gz


fasterq-dump --split-files SRR17055867.sra

#геном Ochotona_princeps
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Ochotona_princeps/latest_assembly_versions/GCF_014633375.1_OchPri4.0/ ~/cow/oc 

#удаление пцр-дупликатов(rd)
export JAVA_HOME=/mss_users/ltursunova/jdk/jdk-17.0.7/
export PATH="$JAVA_HOME/bin:$PATH"
java -jar ~/picard/picard.jar MarkDuplicatesWithMateCigar REMOVE_DUPLICATES=true       I=lep_sort.bam       O=rdlep.bam       M=mark_dups_w_mate_cig_metrics.txt


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

for var in SRR10011655	SRR11020300	SRR17908655	SRR5949623	SRR6485281	SRR17044867	SRR5949630	SRR6485284	SRR11020211	SRR17908654	SRR17908659	SRR17908658	SRR5949632 
do 
scp oc2* ./$var
cd $var
bowtie2 -p 11 -q --very-sensitive-local -x oc2 -U output_paired.pair1.truncated output_paired.pair2.truncated  -S lep.sam
samtools view -bt oc_genome.fa.fai -o $var.bam lep.sam
samtools flagstat -@ 10 $var.bam > flagstat.txt
cd ..
scp ./$var/$var.bam ./New
done

samtools sort SRR10011655.bam -o SRR10011655_sort.bam
samtools index SRR10011655_sort.bam
samtools sort SRR17044867.bam -o SRR17044867_sort.bam
samtools index SRR17044867_sort.bam
samtools sort SRR17908655.bam -o SRR17908655_sort.bam
samtools index SRR17908655_sort.bam
samtools sort SRR17908659.bam -o SRR17908659_sort.bam
samtools index SRR17908659_sort.bam
samtools sort SRR5949630.bam -o SRR5949630_sort.bam
samtools index SRR5949630_sort.bam
bcftools mpileup -f oc2_genome.fa SRR5949630_sort.bam SRR17908659_sort.bam SRR17908655_sort.bam SRR17044867_sort.bam SRR10011655_sort.bam | bcftools call -mv -Oz -o calls.vcf.gz


bcftools filter calls.vcf.gz -s LowQual -e 'QUAL<20 && DP<10' > filtered.vcf





