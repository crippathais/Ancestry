

###################################
# 1.5 Filtering with 1KGP
###################################

# Before estimate the Principal Components (PCs) from the total sample, we need to prune SNPs with Linkage Disequilibrium (LD), in order to avoid bias in PCs due to close SNPs with same frequency. We used the following parameters for pruning: window size=50 SNPs, shift step=5 SNPs, and r2=0.5.

# exclude variants with AF<0.01
plink --bfile array_SABE_1KPG_Natives --maf 0.01 --make-bed --out array_SABE_1KPG_Natives_MAF

# pruning SNPs by LD
plink --bfile array_SABE_1KPG_Natives_MAF --out array_SABE_1KPG_Natives_MAF_LD --indep-pairwise 50 5 0.5
  
# mantain pruned SNPs only and allele frequencies > 0.01
plink --bfile array_SABE_1KPG_Natives_MAF --extract array_SABE_1KPG_Natives_MAF_LD.prune.in --make-bed --out array_SABE_1KPG_Natives_MAF_LD

#Relatedness among samples is analyzed by estimating the proportion of identical-by-descent (IBD) alleles between pairs of individuals by $\hat{pi}$ and relatedness matrix (GRM) estimations. In this context, it is expected that third-degree relatives have values of $\hat{pi}$ or GRM = 0.125. Therefore, pairs of samples with $\hat{pi}$ or GRM > 0.125 are removed from the data.
# The expectation is that IBD = 1 for duplicates or monozygotic twins, IBD = 0.5 for first-degree relatives, IBD = 0.25 for second-degree relatives and IBD = 0.125 for third-degree relatives. Due to genotyping error, LD and population structure there is often some variation around these theoretical values and it is typical to remove one individual from each pair with an IBD > 0.1875, which is halfway between the expected IBD for third- and second-degree relatives. For these same reasons an IBD > 0.98 identifies duplicates.

# estimates IBD from the sample by PLINK
plink --bfile array_SABE_1KPG_Natives_MAF_LD --allow-extra-chr --allow-no-sex --genome --out array_SABE_1KPG_Natives_MAF_LD_temp

#Estimate Missigness
plink --bfile array_SABE_1KPG_Natives_MAF_LD --missing --out array_SABE_1KPG_Natives_MAF_LD_temp.Miss

cp array_SABE_1KPG_Natives_MAF_LD_temp.genome test.genome
cp array_SABE_1KPG_Natives_MAF_LD_temp.Miss.imiss test.imiss

# Run to identify all pairs of individuals with IBD > 0.125
# change directories, if necessary, in perl script
perl /home/thais/Pharmaco/ancestry/array/BRA_SABE/scripts/run-IBD-QC.pl test

plink --bfile array_SABE_1KPG_Natives_MAF_LD --double-id --biallelic-only strict --allow-extra-chr --allow-no-sex --remove /home/thais/Pharmaco/ancestry/array/BRA/ibd/fail-IBD-QC.txt --make-bed --out array_SABE_1KPG_Natives_MAF_LD_IBD

plink --bfile array_SABE_1KPG_Natives_MAF_LD_IBD --double-id --biallelic-only strict --allow-extra-chr --allow-no-sex --remove /home/thais/Pharmaco/ancestry/array/BRA/ibd/fail-IBD-QC.txt --recode vcf --out array_SABE_1KPG_Natives_MAF_LD_IBD

mv test* ibd/
rm array_SABE_1KPG_Natives_MAF_LD_temp*

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
for K in 1 2 3 4 5 6 7 8; \
do ~/Pharmaco/ancestry/array/BRA/popGenomic/scripts/dist/admixture_linux-1.3.0/admixture --cv /home/thais/Pharmaco/ancestry/array/BRA_SABE/array_SABE_1KPG_Natives_MAF_LD_IBD.bed $K  | tee log${K}.out; done

## Seeing the best k for the population
grep -h CV log*.out

## Make a graph in R
R CMD PlotAdmix.R

###
# FST
###

bcftools query -l array_SABE_1KPG_Natives_MAF_LD_IBD.vcf > samplesTotal.txt

grep "^BRA" samplesTotal.txt > BRA.txt
grep "^NAT_BRA" samplesTotal.txt > NAT_BRA.txt
grep "^EUR" samplesTotal.txt > EUR.txt
grep "^AFR" samplesTotal.txt > AFR.txt
grep "^NAT_PER" samplesTotal.txt > NAT_PER.txt
grep "^SAS" samplesTotal.txt > SAS.txt
grep "^EAS" samplesTotal.txt > EAS.txt
grep "^AMR" samplesTotal.txt > AMR.txt

rm samplesTotal.txt

for m in *.txt ; do \
for s in *.txt; do \
vcftools --vcf /home/thais/Pharmaco/ancestry/array/BRA_SABE/array_SABE_1KPG_Natives_MAF_LD_IBD.vcf --weir-fst-pop ${m} --weir-fst-pop ${s} --out ${m}_${s}
done
done

## Mean FST values
grep "Weir and Cockerham mean Fst estimate:" *.log > fst.txt

###
# AMOVA
###

###
# F4
###


