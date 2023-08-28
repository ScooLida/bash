for var in SRR10011655	SRR11020300	SRR17908655	SRR5949623	SRR6485281	SRR17044867	SRR5949630	SRR6485284	SRR11020211	SRR17908654	SRR17908659	SRR17908658	SRR5949632 

export PATH="$HOME/miniconda/bin:$PATH"
export PATH="$HOME/.local/bin:$PATH"

screen -wipe #очистка замерших сессий

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
bcftools mpileup -f ~/cow/oc2_genome.fa rd_z3_sort.bam | bcftools call -mv -Oz -o rd_z3.vcf.gz

#фильтрация
bcftools filter out.vcf.gz -s LowQual -e 'QUAL<20 && DP<10' > filtered.vcf

#статистика
samtools flagstat -@ 10 dna.bam > flagstat.txt 

#удаление пцр-дуликатов
export JAVA_HOME=/mss_users/ltursunova/jdk/jdk-17.0.7/
export PATH="$JAVA_HOME/bin:$PATH"
java -jar ~/picard/picard.jar MarkDuplicatesWithMateCigar REMOVE_DUPLICATES=true       I=z5_sort.bam      O=rd_z5_sort.bam       M=mark_dups_w_mate_cig_metrics.txt
samtools flagstat rd_z5_sort.bam > rd_z5_sort.txt

[E::bgzf_read_block] Failed to read BGZF block data at offset 2240401020 expected 13782 bytes; hread returned 0
[E::bgzf_read] Read block operation failed with error 4 after 30 of 39 bytes

samtools quickcheck -qvvv z11207_2.bam 
#rd_z1_sort.bam was missing EOF block when one should be present.
#rd_z1_sort.bam has 9351 targets in header.
#z11207_2.bam cannot be checked for EOF block because its filetype does not contain one.

#z1_sort.bam has good EOF block.

# bam > sam > bam
# bam  > bam

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
for var in SRR10011655	SRR11020300	SRR17908655	SRR5949623	SRR6485281	SRR17044867	SRR5949630	SRR6485284	SRR11020211	SRR17908654	SRR17908659	SRR17908658	SRR5949632
do 
tabix -p vcf $var.vcf.gz
done
vcf-merge SRR10011655.vcf.gz SRR11020300.vcf.gz SRR17908655.vcf.gz SRR5949623.vcf.gz SRR6485281.vcf.gz SRR17044867.vcf.gz SRR5949630.vcf.gz SRR6485284.vcf.gz SRR11020211.vcf.gz SRR17908654.vcf.gz SRR17908659.vcf.gz SRR17908658.vcf.gz SRR5949632.vcf.gz | bgzip -c > lepus.vcf.gz
                       

 WARNING: The script paleomix is installed in '/mss_users/ltursunova/.local/bin' which is not on PATH.
  Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.




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


