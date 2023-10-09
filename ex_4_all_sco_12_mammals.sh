## data
#wget https://ftp.ensembl.org/pub/release-110/fasta/physeter_catodon/pep/Physeter_catodon.ASM283717v2.pep.all.fa.gz
#wget https://ftp.ensembl.org/pub/release-110/fasta/ovis_aries_rambouillet/pep/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.pep.all.fa.gz
#wget https://ftp.ensembl.org/pub/release-110/fasta/equus_caballus/pep/Equus_caballus.EquCab3.0.pep.all.fa.gz
#wget https://ftp.ensembl.org/pub/release-110/fasta/felis_catus/pep/Felis_catus.Felis_catus_9.0.pep.all.fa.gz
#wget https://ftp.ensembl.org/pub/release-110/fasta/canis_lupus_familiaris/pep/Canis_lupus_familiaris.ROS_Cfam_1.0.pep.all.fa.gz
#wget https://ftp.ensembl.org/pub/release-110/fasta/myotis_lucifugus/pep/Myotis_lucifugus.Myoluc2.0.pep.all.fa.gz
#wget https://ftp.ensembl.org/pub/release-110/fasta/mus_spretus/pep/Mus_spretus.SPRET_EiJ_v1.pep.all.fa.gz
#wget https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz
#wget https://ftp.ensembl.org/pub/release-110/fasta/marmota_marmota_marmota/pep/Marmota_marmota_marmota.marMar2.1.pep.all.fa.gz
#wget https://ftp.ensembl.org/pub/release-110/fasta/macaca_mulatta/pep/Macaca_mulatta.Mmul_10.pep.all.fa.gz
#wget https://ftp.ensembl.org/pub/release-110/fasta/ornithorhynchus_anatinus/pep/Ornithorhynchus_anatinus.mOrnAna1.p.v1.pep.all.fa.gz
#wget https://ftp.ensembl.org/pub/release-110/fasta/vombatus_ursinus/pep/Vombatus_ursinus.bare-nosed_wombat_genome_assembly.pep.all.fa.gz

mkdir 1_pep_all/
cd 1_pep_all

cp /home/zin_data/*all.fa.gz 1_pep_all/

#РАЗАРХИВИРОВАТЬ ФАСТУ


# full set of single copy orthologs
## proteinortho has trouble with long names
for file in ./1_pep_all/*fa; do sed -e's/\ .*//' $file >$file.clean.fa; done
mkdir 2_pep_renamed
mv ./1_pep_all/*clean* 2_pep_renamed/ 
## long command!!!
proteinortho 2_pep_renamed/*fa -t 12
## proteinortho runs for around 3h at one processor

grep -v  \* myproject.proteinortho.tsv  | grep -v -c "," 
#686 genes

cat 2_pep_renamed/*fa >all_pep.fa

## let's test bbmap ## had to update bbmap!
#filterbyname.sh in=../all_pep.fa out=family1.aa.fa include=t names=family1.names.txt overwrite=true ignorejunk=true
for list in family*txt
do filterbyname.sh in=../all_pep.fa out=$list.fa include=t names=$list overwrite=true ignorejunk=true
done

## now rename...
mkdir ../families_seqs
for file in *fa
do sed -e 's/ENSCAFP.*/Canis_lupus_familiaris/' $file | sed -e 's/ENSECAP.*/Equus_caballus/' | sed -e 's/ENSFCAP.*/Felis_catus/' | sed -e 's/ENSMMUP.*/Macaca_mulatta/' | sed -e 's/ENSMMMP.*/Marmota_marmota/' | sed -e 's/ENSMUSP.*/Mus_musculus/' | sed -e 's/MGP_SPRETEiJ.*/Mus_spretus/' | sed -e 's/ENSMLUP.*/Myotis_lucifugus/' | sed -e 's/ENSOANP.*/Ornithorhynchus_anatinus/' | sed -e 's/ENSOARP.*/Ovis_aries/' | sed -e 's/ENSPCTP.*/Physeter_catodon/' | sed -e 's/ENSVURP.*/Vombatus_ursinus/' > ../families_seqs/$file
done

cd ../
mkdir families_alns/
 for file in *fa; do mafft --auto $file >../families_alns/$file.aln; done
rm families_alns/family79.names.txt.fa.aln 
## this is for building the tree without trimming
#iqtree2 -s families_alns/ -alrt 1000 -abayes -o Ornithorhynchus_anatinus,Vombatus_ursinus --redo

## we can also trim 
mkdfir families_trimmed_alns/
cd families_alns/
for file in *aln; do in trimal -in $file -out ../families_trimmed_alns/$file -automated1 
## long command!!!
iqtree2 -s families_trimmed_alns/ -alrt 1000 -abayes -o Ornithorhynchus_anatinus,Vombatus_ursinus -nt 12
