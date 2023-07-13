for var in SRR10011655	SRR11020300	SRR17908655	SRR5949623	SRR6485281	SRR17044867	SRR5949630	SRR6485284	SRR11020211	SRR17908654	SRR17908659	SRR17908658	SRR5949632 

export PATH="$HOME/miniconda/bin:$PATH"

screen -wipe #очистка замерших сессий

#копирование
scp /home/scoollida/picard.jar ltursunova@sy2.computing.kiae.ru:~/picard
scp ltursunova@sy2.computing.kiae.ru:~/don/TESTR/pca2.eigenval /home/scoollida/
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
samtools view -bt ~/cow/oc2_genome.fa.fai -@ 12 -o z1_1207.bam z1.sam
samtools flagstat -@ 10 z1_1207.bam > z1_1207.txt 
samtools sort lep.bam -o lep_sort.bam
samtools index lep_sort.bam
bcftools mpileup -f oc2_genome.fa lep_sort.bam | bcftools call -mv -Oz -o calls.vcf.gz

bcftools view -m2 -M2 -v snps input.vcf.gz
 tabix -p vcf zz1.vcf.gz
vcf-merge A.vcf.gz B.vcf.gz C.vcf.gz | bgzip -c > out.vcf.gz

#фильтрация
bcftools filter out.vcf.gz -s LowQual -e 'QUAL<20 && DP<10' > filtered.vcf

#корявый bam в норм sam
samtools view -@ 12 -h rd_z1_sort.bam > z11207_2.bam 
samtools view -@ 12 -bt ~/cow/oc2_genome.fa.fai  -o .bam file1.sam
#статистика
samtools flagstat -@ 10 dna.bam > flagstat.txt 


#получаем формат .traw 
~/plink --vcf out.vcf.gz --recode A-transpose --out Z_filtered

#################################################################################################

samtools quickcheck -qvvv z11207_2.bam 
z11207_2.bam cannot be checked for EOF block because its filetype does not contain one.
 z1_sort.bam
 z1_sort.bam has good EOF block.
rd_z1_sort.bam was missing EOF block when one should be present.

~/plink --vcf out.vcf  --allow-extra-chr --make-bed --out  1307


samtools view -bt ~/cow/oc2_genome.fa.fai  -o z3.bam z3.sam
samtools sort z3.bam -o z3_sort.bam
export JAVA_HOME=/mss_users/ltursunova/jdk/jdk-17.0.7/
export PATH="$JAVA_HOME/bin:$PATH"
java -jar ~/picard/picard.jar MarkDuplicatesWithMateCigar REMOVE_DUPLICATES=true       I=z5_sort.bam      O=rd_z5_sort.bam       M=mark_dups_w_mate_cig_metrics.txt
samtools flagstat rd_z5_sort.bam > rd_z5_sort.txt 

bcftools filter calls.vcf.gz -s LowQual -e 'QUAL<20 && DP<10' > Z_filtered.vcf
 
bcftools view -M2 -v snps out.vcf.gz #убирает больше 2х аллелей, тк bim не записывался
plink --vcf 2all --allow-extra-chr --make-bed --out  2all
~/plink --bfile 2all --pca 2 --allow-extra-chr --bad-freqs  --out pca2
