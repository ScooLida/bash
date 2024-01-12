#https://devonderaad.github.io/aph.rad/admixture/plot.admixture.results.html
vcftools --vcf re.12.01.oc.miss02.vcf --plink-tped --out re.12.01.oc.miss02
~/plink --tped re.12.01.oc.miss02.tped --tfam re.12.01.oc.miss02.tfam --make-bed --allow-extra-chr --out re.12.01.oc.miss02
for K in echo $(seq 22) ; do ~/admixure/dist/admixture_linux-1.3.0/admixture --cv=10 -B2000 -j8 re.12.01.oc.miss02.bed $K | tee log${K}.out; done
