                                                                                    
fasterq-dump --split-files SRR17055867.sra

#геном Ochotona curzoniae
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/696/305/GCF_001696305.1_UCN72.1 my_dir/ 

#удаление адаптеров
export PATH="$HOME/miniconda/bin:$PATH"
conda install -c bioconda adapterremoval
AdapterRemoval --file1 1k_S1_R1_001.fastq  --file2 1k_S1_R2_001.fastq --basename output_paired --trimns --trimqualities

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

bgzip -c SRR10011655.vcf > SRR10011655.vcf.gz

#должен был получиться формат eigenstrat
vcf_format_conversions.py --vcf SRR10011655.vcf --out SRR10011655 --out-format eigenstrat

#получаем формат .traw 
~/plink --vcf SRR10011655.vcf.gz --recode A-transpose --allow-extra-chr --out SRR10011655


library (smartsnp, lib.loc = "~/_lib")

My <- c(rep("A",2), rep("B",2))
sm.pca <- smart_pca(snp_data = "1.traw", "2.traw", "3.traw", "SRR10011655.traw",
                    missing_value = NA,
                        sample_group = My)
pdf("my_plot.pdf")
plot(sm.pca$pca.sample_coordinates[, c(3,4)])
dev.off()
 
Checking argument options selected...
Argument options are correct...
Loading data...
Imported 11426973 SNP by 4 sample genotype matrix
Time elapsed: 0h 0m 32s
Filtering data...
No samples projected after PCA computation
Error in smart_pca(snp_data = "SRR10011655.traw", "1.traw", "2.traw",  : 
  length(sample_group) should be equal to number of samples in dataset: computation aborted


library (smartsnp, lib.loc = "~/_lib")
My <- c(rep("A",1), rep("B",2))
hh <- c("/mss_users/ltursunova/cow/R/PPP/1.traw", "/mss_users/ltursunova/cow/R/PPP/3.traw", "/mss_users/ltursunova/cow/R/PPP/2.traw",
sm.pca <- smart_pca(snp_data =hh,
                    missing_value = NA,
                        sample_group = My)
pdf("my_plot.pdf")
plot(sm.pca$pca.sample_coordinates[, c(3,4)])
dev.off()

Error in file(filepath, "rb", raw = TRUE) : 
  invalid 'description' argument
Calls: smart_pca -> is.binary -> file
Execution halted


#проверка на наличие 
if test -f "$var.vcf"; then
    echo "$var exists."
else echo "$var not"
fi  

#показать строку
sed '2221,2222!d' SRR17908658.vcf
