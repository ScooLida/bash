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




