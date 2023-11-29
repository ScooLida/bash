with open('new5.txt', 'r') as file_in, open ('new7.txt', 'w') as out:
 for line in file_in:
    # Разделяем строку на отдельные слова\n",
    line2 = line.replace("##bcftoolsCommand=mpileup -f /mss_users/ltursunova/cow/new2/och.fasta", '') 

     if '##' not in line2:
      line3=line2.replace('.sort.bam','')
      print(line2)
      out.write(line2)

