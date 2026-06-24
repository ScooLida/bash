#!/bin/bash
DATA = "data"                                                                                            
#geno 0 - не будет вообще пропусков
#--keep samples_to_keep.txt, в котором образец должен быть написан на отдельной строке два раза через табуляцию
~/plink --vcf ${DATA}.vcf.gz --geno 0.2          --make-bed      --threads 8  --out ${DATA}_pca   --allow-extra-chr

~/plink --bfile ${DATA}_pca       --pca 4    --threads 8     --allow-extra-chr   --out ${DATA}_pca


#~/plink --bfile ${DATA}_pca  --make-bed --allow-extra-chr  --threads 8  --out ${DATA}_admixture

#k - число с минимальным количеством предполагаемых популяций до максимального
for K in $(seq 5 15)
do ~/admixture/dist/admixture_linux-1.3.0/admixture -j4 --cv ${DATA}_admixture.bed $K | tee log${K}.out; 
done
