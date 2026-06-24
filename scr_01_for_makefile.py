#This script takes a list of samples and generates a list of file paths for fastq or fastq.gz data for a PALEOMIX makefile
with open('list_11_12.csv', 'r') as file_in, open ('listmk.csv', 'w') as out:
 for line in file_in:
     line = line.strip()
     form = ("  "+line+":"+"\n"+"   "+line+":"+"\n"+"    "+"Path:  /home/users/"+line+"*fastq*"+"\n")
     print(form)
     out.write(form)
