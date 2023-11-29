 "source": [
    "with open('new5.txt', 'r') as file_in, open ('new7.txt', 'w') as out:\n",
    "    for line in file_in:\n",
    "# Разделяем строку на отдельные слова\n",
    "        line2 = line.replace(\"##bcftoolsCommand=mpileup -f /mss_users/ltursunova/cow/new2/och.fasta \", '')        \n",
    "        \n",
    "        if '##' not in line2:\n",
    "            print(line2)\n",
    "            out.write(line2)\n",
    "            \n",
    "\n"
   ],
