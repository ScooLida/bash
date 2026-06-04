with open('list_11_12.csv', 'r') as file_in, open ('listmk.csv', 'w') as out:
 for line in file_in:
     line = line.strip()
     form = ("  "+line+":"+"\n"+"   "+line+":"+"\n"+"    "+"Path:  /home/users/"+line+"*fastq.gz"+"\n")
     print(form)
     out.write(form)
