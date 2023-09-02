export PATH="$HOME/miniconda/bin:$PATH"
export PATH="$HOME/.local/bin:$PATH"

#копирование
scp /home/scoollida/picard.jar ltursunova@sy2.computing.kiae.ru:~/picard
scp ltursunova@sy2.computing.kiae.ru:~/don/TESTR/pca2.eigenval /home/scoollida/

#проверка на наличие 
do
    cd  $var
    if ls *.fastq >/dev/null 2>&1; then
        echo "$var exists."
    else
        echo "$var does not exist."
    fi
    cd ..
done
#Здесь использована команда ls *.fastq >/dev/null 2>&1, 
#которая проверяет наличие файлов с расширением .fastq в текущей директории. 
#Если файлы существуют, вывод команды перенаправляется в /dev/null, 
#чтобы не показывать список файлов. Если файлов нет, вывод перенаправляется в /dev/null и вместо этого выводится сообщение "does not exist."

#показать строку
sed '2221,2222!d' SRR17908658.vcf

#архив
bgzip -c Z_filtered.vcf > Z_filtered.vcf.gz
gzip -dk SRR10011655.vcf.gz
tar xvzf archive.tar.gz

#скачивание файлов с ncbi
fasterq-dump --split-files SRR17055867.sra

#геном Ochotona_princeps
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Ochotona_princeps/latest_assembly_versions/GCF_014633375.1_OchPri4.0/ ~/cow/oc 

#получение bam в paleomix
paleomix bam run makefile.yaml --max-threads 15

# форматирование samtools
samtools faidx ref.fa
samtools view -bt ~/cow/oc2_genome.fa.fai -@ 12 -o z1_1207.bam z1.sam
samtools flagstat -@ 10 z1_1207.bam > z1_1207.txt 
samtools sort lep.bam -o lep_sort.bam
samtools index lep_sort.bam
bcftools mpileup -f ~/cow/oc2_genome.fa rd_z3_sort.bam | bcftools call -mv -Oz -o rd_z3.vcf.gz

#фильтрация
bcftools filter out.vcf.gz -s LowQual -e 'QUAL<20 && DP<10' > filtered.vcf

#статистика
samtools flagstat -@ 10 dna.bam > flagstat.txt 

#соединяет vcf файлы
tabix -p vcf zz1.vcf.gz
vcf-merge zz1.vcf.gz zz4.vcf.gz zz5.vcf.gz | bgzip -c > Allzz.vcf.gz

#чистит от мультиаллельных локусов, оставляет 1-2
bcftools view  -M2 -v snps input.vcf.gz

#получаем нужные bim fam bed
~/plink --vcf out.vcf  --allow-extra-chr --make-bed --out  1307

#выполняем pca
~/plink --bfile 1407 --pca 4 --allow-extra-chr --bad-freqs  --out pcagen




############################################################################################################

#удаление адаптеров(есть в paleomix)
export PATH="$HOME/miniconda/bin:$PATH"
conda install -c bioconda adapterremoval
AdapterRemoval --file1 1k_S1_R1_001.fastq  --file2 1k_S1_R2_001.fastq --basename output_paired --trimns --trimqualities

#индекс пищухи(есть в paleomix)
bowtie2-build reference_sequence.fasta index_name

#выравнивание на геном пищухи bowtie. получение .sam(есть в paleomix)
bowtie2 -p 6 -q --very-sensitive-local -x oc_index -U output_paired.pair1.truncated output_paired.pair2.truncated  -S lep.sam


#удаление пцр-дуликатов(есть в paleomix)
export JAVA_HOME=/mss_users/ltursunova/jdk/jdk-17.0.7/
export PATH="$JAVA_HOME/bin:$PATH"
java -jar ~/picard/picard.jar MarkDuplicatesWithMateCigar REMOVE_DUPLICATES=true       I=z5_sort.bam      O=rd_z5_sort.bam       M=mark_dups_w_mate_cig_metrics.txt
samtools flagstat rd_z5_sort.bam > rd_z5_sort.txt


MyFilename:
    MySample:
      SRR10011655:
        Path:  /mss_users/ltursunova/cow/new2/SRR10011655*.fastq
    MySample2:
       SRR12518920:
        Path:  /mss_users/ltursunova/cow/new2/SRR12518920*.fastq
    MySample3:
       SRR10012547:
        Path:  /mss_users/ltursunova/cow/new2/SRR10012547*.fastq



    SRR10011655: /mss_users/ltursunova/cow/new2/SRR10011655*.fastq
    SRR12518920: /mss_users/ltursunova/cow/new2/SRR12518920*.fastq
    SRR10012547: /mss_users/ltursunova/cow/new2/SRR10012547*.fastq
    SRR10082097: /mss_users/ltursunova/cow/new2/SRR10082097*.fastq
    SRR10082089: /mss_users/ltursunova/cow/new2/SRR10082089*.fastq
    SRR17908653: /mss_users/ltursunova/cow/new2/SRR17908653*.fastq
    SRR5949621: /mss_users/ltursunova/cow/new2/SRR5949621*.fastq
    SRR5949624: /mss_users/ltursunova/cow/new2/SRR5949624*.fastq
    SRR5949634: /mss_users/ltursunova/cow/new2/SRR5949634*.fastq
    SRR6485240: /mss_users/ltursunova/cow/new2/SRR6485240*.fastq
    SRR6485265: /mss_users/ltursunova/cow/new2/SRR6485265*.fastq


