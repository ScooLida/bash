
######################################
# Principal Component Analysis - PCA #
######################################

## Genetic distances between individuals
system("plink --cow --allow-no-sex --nonfounders --file ADAPTmap_TOP --distance-matrix --out dataForPCA")

## Load data
dist_populations<-read.table("dataForPCA.mdist",header=F)
### Extract breed names
fam <- data.frame(famids=read.table("dataForPCA.mdist.id")[,1])
### Extract individual names 
famInd <- data.frame(IID=read.table("dataForPCA.mdist.id")[,2])
  
## Perform PCA using the cmdscale function 
# Time intensive step - takes a few minutes with the 4.5K animals
mds_populations <- cmdscale(dist_populations,eig=T,5)

## Extract the eigen vectors
eigenvec_populations <- cbind(fam,famInd,mds_populations$points)

## Proportion of variation captured by each eigen vector
eigen_percent <- round(((mds_populations$eig)/sum(mds_populations$eig))*100,2)
