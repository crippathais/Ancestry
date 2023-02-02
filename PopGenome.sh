###################################
# 1.6 Evaluating population structure 
# PCA
# ADMIXTURE
# FST
# AMOVA
# F4
###################################

###
# PCA
###

# estimates PCs Total
plink --bfile array_SABE_1KPG_Natives_MAF_LD_IBD --pca 3333 header --allow-no-sex --out array_SABE_1KPG_Natives_MAF_LD_IBD

# extract the 10 first PCs
cut -d ' ' -f 1-12 array_SABE_1KPG_Natives_MAF_LD_IBD.eigenvec > temp.eigenvec
mv temp.eigenvec array_SABE_1KPG_Natives_MAF_LD_IBD.eigenvec

# cut -d ' ' -f 1-12 alldatahg38_MAF_LD.eigenvec > temp.eigenvec
# mv temp.eigenvec alldatahg38_MAF_LD.eigenvec

#### Plot the PCA results

R CMD PlotPCA.R

###
# ADMIXTURE
###

# Running ADMIXTURE
for K in 1 2 3 4 5 6 7 8 9 10; do ~/Pharmaco/ancestry/array/BRA/popGenomic/scripts/dist/admixture_linux-1.3.0/admixture --cv /home/thais/Pharmaco/ancestry/arrayExome/array_SABE_1KPG_Natives_MAF_LD_IBD.bed $K | tee log${K}.out; done

## Seeing the best k for the population
grep -h CV log*.out

## Make a graph in R
R CMD PlotAdmix.R

###
# FST
###

## FST by method Hudson

~/bin/plink2 --bfile /home/thais/Pharmaco/ancestry/arrayExome/array_SABE_1KPG_Natives_MAF_LD_IBD --family --fst CATPHENO method=hudson --out FSTHudsonTotal

###
# AMOVA
###







