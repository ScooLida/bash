for var in SRR10011655	SRR11020300	SRR17908655	SRR5949623	SRR6485281	SRR17044867	SRR5949630	SRR6485284	SRR11020211	SRR17908654	SRR17908659	SRR17908658	SRR5949632 

export PATH="$HOME/miniconda/bin:$PATH"

#проверка на наличие 
if test -f "$var.vcf"; then
    echo "$var exists."
else echo "$var not"
fi  

#показать строку
sed '2221,2222!d' SRR17908658.vcf

#архив
bgzip -c Z_filtered.vcf > Z_filtered.vcf.gz
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
samtools faidx ref.fa
samtools view -bt ~/cow/oc2_genome.fa.fai  -o lep.bam lep.sam
samtools sort lep.bam -o lep_sort.bam
samtools index lep_sort.bam
bcftools mpileup -f oc2_genome.fa lep_sort.bam | bcftools call -mv -Oz -o calls.vcf.gz
vcf-merge A.vcf.gz B.vcf.gz C.vcf.gz | bgzip -c > out.vcf.gz
 tabix -p vcf zz1.vcf.gz
#фильтрация
bcftools filter out.vcf.gz -s LowQual -e 'QUAL<20 && DP<10' > filtered.vcf

#статистика
samtools flagstat -@ 10 dna.bam > flagstat.txt 


#получаем формат .traw 
~/plink --vcf out.vcf.gz --recode A-transpose --out Z_filtered

#################################################################################################



bcftools mpileup -f oc2_genome.fa SRR10011655.sort.bam  SRR11020300.sort.bam     SRR17908655.sort.bam     SRR5949623.sort.bam      SRR6485281.sort.bam      SRR17044867.sort.bam     SRR5949630.sort.bam      SRR6485284.sort.bam      SRR11020211.sort.bam     SRR17908654.sort.bam     SRR17908659.sort.bam     SRR17908658.sort.bam     SRR5949632.sort.bam  | bcftools call -mv -Oz -o calls2.vcf.gz

samtools view -bt ~/cow/oc2_genome.fa.fai  -o z3.bam z3.sam
samtools sort z3.bam -o z3_sort.bam
export JAVA_HOME=/mss_users/ltursunova/jdk/jdk-17.0.7/
export PATH="$JAVA_HOME/bin:$PATH"
java -jar ~/picard/picard.jar MarkDuplicatesWithMateCigar REMOVE_DUPLICATES=true       I=z3_sort.bam       O=rd_z3_sort.bam       M=mark_dups_w_mate_cig_metrics.txt
samtools index rd_z3_sort.bam 
samtools flagstat rd_z3_sort.bam  > rd_z3_sort.txt 
samtools flagstat rd_z3.bam > flagstat.txt 

cd ./TESTR
bcftools mpileup -f oc_genome.fa 1.bam 3.bam 4.bam 5.bam | bcftools call -mv -Oz -o calls.vcf.gz
bcftools filter calls.vcf.gz -s LowQual -e 'QUAL<20 && DP<10' > Z_filtered.vcf





