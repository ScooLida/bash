#!/bin/bash
# The script performs PCA and admixture analysis on a VCF file.
# It converts the VCF to PLINK binary format, runs PCA, and then runs
# admixture for a range of assumed ancestral populations (K).

# Prefix for input and output file names
DATA="data"

# PLINK binary (change the path if necessary)
PLINK="~/plink"

# Admixture binary (change the path if necessary)
ADMIXTURE="~/admixture/dist/admixture_linux-1.3.0/admixture"

# --geno 0.2 removes variants with >20% missing genotypes
# (use --geno 0 to remove all variants with any missing data)

# --keep samples_to_keep.txt can be added here;
# in that file each sample must be written on a separate line twice,
# separated by a tab.
$PLINK --vcf ${DATA}.vcf.gz --geno 0.2 --make-bed --threads 8 --out ${DATA}_pca --allow-extra-chr

$PLINK --bfile ${DATA}_pca --pca 4 --threads 8 --allow-extra-chr --out ${DATA}_pca

# Uncomment the line below to re-create the bed file specifically for admixture
# $PLINK --bfile ${DATA}_pca --make-bed --allow-extra-chr --threads 8 --out ${DATA}_admixture

# K is the assumed number of ancestral populations, from minimum to maximum
for K in $(seq 5 15)
do
  $ADMIXTURE -j4 --cv ${DATA}_admixture.bed $K | tee log${K}.out
done
