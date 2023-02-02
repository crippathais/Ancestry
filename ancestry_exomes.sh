############################################################################
#############################################################################
#
# Ancestry haplotypes pipeline
#
# AUTHOR: Thais de Oliveira (PostDoc at FMC - UNICAMP), Rodrigo Secolin (PostDoc at FMC - UNICAMP)
#
#
# INPUTS: .vcf file name and full path
#         file with sample status (mapping sample_id to status and gender)
#         full path of output directory
#
# REQUIREMENTS:
#          PLINK
#          BCFTOOLS
#          R
#          Perl
#
#
# Last modified: 14.Sep.2022
#
#############################################################################
#############################################################################

#Arrays lab
VCF="/home/thais/raw/array/liftover/snparrayTotal_Filter.vcf.gz"
#SABE genome
VCF="/home/thais/raw/SABE/SABE1172.1a1172.CONCAT.hg38.VQSR-AS_b.PASS.recode.vcf.gz"
###################################
# Steps:
# 1.1 Preparing Data
# 1.2 Filtering by Genotypes and Samples
# 1.3 Creating 1KGP SNP list -- if snparray or exomes as input 
# 1.4 Merging with 1000 Genome Project
# 1.4 - part 2 Merging with Natives from Brazil
# 1.5 Filtering with 1KGP
# 1.6 Evaluating population structure by Principal Component analysis (PCA)
# 2 Haplotype phasing by SHAPEIT4

# Scripts:
# FilterGenotypes.Rmd 

###################################


###################################
# 1.1 Preparing Data
###################################

# Convert the VCF files to PLINK FILES (BED/BIM/FAM). Here the missing SNP ids are replaced by chromosome and position. 

### For all VCF files
plink2 --vcf ${VCF} --allow-extra-chr --max-alleles 2 --const-fid --set-all-var-ids @:#:\$1:\$2 --snps-only --make-bed --out hg38

## If Natives
plink --vcf ${VCF} --vcf-half-call missing --allow-extra-chr --biallelic-only strict --const-fid --set-all-var-ids @:#:\$1:\$2 --new-id-max-allele-len 10 --snps-only --make-bed --out hg38


###################################
# 1.2 Filtering by Genotypes and Samples
###################################

#Create a list of ambiguous SNPs, which are those with A/T or G/C genotypes. Remove non-autosomal SNPs, ambigous SNPs, and HWE p-value < 0.000001 (parameters that have been used by the ENIGMA2 international collaboration for imputation). In this case I did not use missing data filtering, since when we merge individuals exomes, the different variants among them is high, and the missing data is more than 50%

R CMD FilterGenotypes.Rmd

###################################
# 1.3 Creating 1KGP SNP list -- if snparray or exomes as input 
# DO not use with genomics datasets
###################################

# Now, we merge our sampe with 1000 Genome Project data phase 3 in order to overlap SNPs and correct strand orientation. This reference data is also build in hg19. The 1000 genome data is located in PLINK format in /home/nfs/ref/1kgp/ folder.

#### Extract variants found in hg38 dataset from 1KGP dataset.

### Download data from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/

for i in {1..22}; do \
  plink \
    --bfile hg38_QC_het_rel \
    --chr $i \
    --make-just-bim --out hg38_temp; 
  awk '{print $1"\t"$4"\t"$4"\t"$2}' hg38_temp.bim > hg38_temp.txt;
  plink \
    --vcf /home/nfs/ref/1KGP_phase3_hg38/1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
    --extract range hg38_temp.txt \
    --set-missing-var-ids @:#\$1,\$2 \
    --make-bed --out thousand_genome_chr${i}_extract; 
done

# merging 1000 genomes chromosome files, 
ls thousand_genome_chr*extract*log | cut -d'.' -f1 > 1KGPmergelist.txt

plink --merge-list 1KGPmergelist.txt --out ThousandGenomeForHg38

grep "^Warning" ThousandGenomeForHg38.log | awk -F' ' ' {print $3 "\t" $3}' | sed "s/'//g" > RemoveWarnings.txt

# #if there is SNP errors, exclude them
for i in {1..22}; do \
  plink \
    --bfile thousand_genome_chr${i}_extract \
    --exclude RemoveWarnings.txt \
    --make-bed --out thousand_genome_chr${i};
  mv thousand_genome_chr${i}.bed thousand_genome_chr${i}_extract.bed;
  mv thousand_genome_chr${i}.bim thousand_genome_chr${i}_extract.bim;
  mv thousand_genome_chr${i}.fam thousand_genome_chr${i}_extract.fam;
done

# merge again
plink --merge-list 1KGPmergelist.txt --out ThousandGenomeForHg38_2

#remove extract files
rm thousand_genome_*

###################################
# 1.4 Merging with 1000 Genome Project
###################################
  
#  Copy 1KGP population info to ThousandGenomeForHg38_2.fam file
#  It's necessary to have in Family ID the correspondent SuperPopulation in each sample are in
# Available at: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/
# file: /home/thais/Pharmaco/ancestry/20130606_g1k_3202_samples_ped_population.txt

# change id's if not in the begining
plink2 --bfile hg38_QC_het_rel --set-all-var-ids @:#:\$1:\$2 --make-bed --out hg38_QC_het_rel_2

# create a snplist from our data
plink --bfile hg38_QC_het_rel_2 --write-snplist --out hg38_QC_het_rel_2

# extract the SNPs from 1000 Genome data
plink --bfile ThousandGenomeForHg38_2 --extract hg38_QC_het_rel_2.snplist --make-bed --out ThousandGenomeForHg38_extract

# create a list of SNPs from 1000genomes to overlap exactly the same SNPs from our sample
plink --bfile ThousandGenomeForHg38_extract --write-snplist --out ThousandGenomes_extract

# extract the 1000 genome SNPs from our sample
plink --bfile hg38_QC_het_rel_2 --extract ThousandGenomes_extract.snplist --make-bed --out hg38_QC_het_rel_extract

# merge our sample with 1000genomes data
plink --bfile hg38_QC_het_rel_extract --bmerge ThousandGenomeForHg38_extract --allow-no-sex --make-bed --out alldatahg38_1KGP

# If still appears some SNPs (often two or three), remove them from our sample and from 1000 genome data
# plink --bfile hg38_QC_het_rel_flip --exclude alldatahg38_1KGP-merge.missnp --make-bed --out hg38_QC_het_rel_exclude

# plink --bfile ThousandGenomeForHg38_extract --exclude alldata1KGP-merge.missnp --make-bed --out ThousandGenomeForHg38_exclude

# # merge again
# plink --bfile hg38_QC_het_rel_exclude --bmerge ThousandGenomeForHg38_exclude --make-bed --out alldatahg38_1KGP

echo "It's worked"

###################################
# 1.4 - part 2 Merging with Natives from Brazil
###################################

#Usinf snplist from step 1.4

NATIVES="/home/thais/raw/Natives/tupiniquim_guarani_Axiom_InCor_BB"
NAME="natives"

# Bed to VCF file and change id's
plink --bfile ${NATIVES} --double-id --allow-no-sex --recode vcf --out ${NAME}
#change to fiz max-alleles=2
vcftools --vcf ${NAME}.vcf --min-alleles 1 --max-alleles 2 --recode --recode-INFO-all --out ${NAME}
# Remove lines that contain "-""
grep -v  "-" natives.recode.vcf > natives_2.recode.vcf

#Liftover Natives Dataset : hg19 to hg38
java -Xms150g -jar /home/thais/Pharmaco/ancestry/array/scripts/picard.jar LiftoverVcf \
I=natives_2.recode.vcf \
O=/home/thais/Pharmaco/ancestry/array/${NAME}.hg19toHg38.vcf \
CHAIN=/home/nfs/ref/Liftover_Chain_Files/b37ToHg38.over.chain \
REJECT=/home/thais/Pharmaco/ancestry/array/${NAME}.hg19toHg38.reject.vcf \
R=/home/nfs/ref/hg38/Homo_sapiens_assembly38.fasta  \
MAX_RECORDS_IN_RAM=40000

#Processed 176003 variants.
#88 variants failed to liftover.
#542 variants lifted over but had mismatching reference alleles after lift over.
#0.3579% of variants were not successfully lifted over and written to the output.
#### Has Y

NATIVE="/home/thais/Pharmaco/ancestry/array/BRA/natives.hg19toHg38.vcf"

# Changes id's
plink2 --vcf ${NATIVE} --allow-extra-chr --double-id --new-id-max-allele-len 2 missing --set-all-var-ids @:#:\$1:\$2 --make-bed --out Natives

# extract the SNPs from Natives Genome data
plink --bfile Natives --allow-extra-chr --double-id --extract hg38_QC_het_rel_2.snplist --make-bed --out NativesForHg38_extract

# create a list of SNPs from Natives Genome data to overlap exactly the same SNPs from our sample
plink --bfile NativesForHg38_extract --allow-extra-chr --double-id --write-snplist --out Natives_extract

# extract the Natives Genome data SNPs from our sample
plink --bfile hg38_QC_het_rel_2 --allow-extra-chr --double-id --extract Natives_extract.snplist --make-bed --out hg38_QC_het_rel_extract_natives

# merge our sample with Natives Genome data
plink --bfile hg38_QC_het_rel_extract_natives --allow-extra-chr --double-id --bmerge NativesForHg38_extract --allow-no-sex --make-bed --out NativesHg38_extractMerge

#### Merging Natives with 1KGP + all snparrays
plink --bfile alldatahg38_1KGP --bmerge NativesHg38_extractMerge --allow-extra-chr --double-id --allow-no-sex --make-bed --out alldatahg38_1KGP_natives

###################################
# 1.5 Filtering with 1KGP
###################################

# Before estimate the Principal Components (PCs) from the total sample, we need to prune SNPs with Linkage Disequilibrium (LD), in order to avoid bias in PCs due to close SNPs with same frequency. We used the following parameters for pruning: window size=50 SNPs, shift step=5 SNPs, and r2=0.5.

# exclude variants with AF<0.01
plink --bfile alldatahg38_1KGP_natives --maf 0.01 --make-bed --out alldatahg38_1KGP_natives_MAF

# pruning SNPs by LD
plink --bfile alldatahg38_1KGP_natives_MAF --out alldatahg38_1KGP_natives_MAF_LD --indep-pairwise 50 5 0.5
  
# mantain pruned SNPs only and allele frequencies > 0.01
plink --bfile alldatahg38_1KGP_natives_MAF --extract alldatahg38_1KGP_natives_MAF_LD.prune.in --make-bed --out alldatahg38_1KGP_natives_MAF_LD

#Relatedness among samples is analyzed by estimating the proportion of identical-by-descent (IBD) alleles between pairs of individuals by $\hat{pi}$ and relatedness matrix (GRM) estimations. In this context, it is expected that third-degree relatives have values of $\hat{pi}$ or GRM = 0.125. Therefore, pairs of samples with $\hat{pi}$ or GRM > 0.125 are removed from the data.
# The expectation is that IBD = 1 for duplicates or monozygotic twins, IBD = 0.5 for first-degree relatives, IBD = 0.25 for second-degree relatives and IBD = 0.125 for third-degree relatives. Due to genotyping error, LD and population structure there is often some variation around these theoretical values and it is typical to remove one individual from each pair with an IBD > 0.1875, which is halfway between the expected IBD for third- and second-degree relatives. For these same reasons an IBD > 0.98 identifies duplicates.

# estimates IBD from the sample by PLINK
plink --bfile alldatahg38_1KGP_natives_MAF_LD --allow-extra-chr --allow-no-sex --genome --out alldatahg38_1KGP_natives_MAF_LD_temp

#Estimate Missigness
plink --bfile alldatahg38_1KGP_natives_MAF_LD --missing --out alldatahg38_1KGP_natives_MAF_LD_temp.Miss

cp alldatahg38_1KGP_natives_MAF_LD_temp.genome test.genome
cp alldatahg38_1KGP_natives_MAF_LD_temp.Miss.imiss test.imiss

# Run to identify all pairs of individuals with IBD > 0.125
# change directories, if necessary, in perl script
perl /home/thais/Pharmaco/ancestry/genome/SABE/scripts/run-IBD-QC.pl test

plink --bfile alldatahg38_1KGP_natives_MAF_LD --double-id --biallelic-only strict --allow-extra-chr --allow-no-sex --remove /home/thais/Pharmaco/ancestry/array/BRA/ibd/fail-IBD-QC.txt --make-bed --out alldatahg38_1KGP_natives_MAF_LD_IBD

plink --bfile alldatahg38_1KGP_natives_MAF_LD_IBD --double-id --biallelic-only strict --allow-extra-chr --allow-no-sex --keep 1kgp_Natives_BRA.txt --make-bed --out alldatahg38_1KGP_natives_MAF_LD_IBD_BRA

mv test* ibd/
rm alldatahg38_1KGP_natives_MAF_LD_temp*


###################################
# 2 Haplotype phasing by SHAPEIT4
###################################

# Input files
# It's imperative to know that for shapeit we do not use the LD and IBD filters
#plink="/home/thais/Pharmaco/ancestry/array/BRA/alldatahg38_1KGP_natives_MAF"
plink="/home/thais/Pharmaco/ancestry/genome/SABE/alldatahg38_1KGP_natives_MAF"
ref="/home/thais/Pharmaco/ref/1KGP_phase3_hg38/1kGP_high_coverage_Illumina"
mapfile="/home/thais/Pharmaco/ancestry/genome/SABE/shapeit/geneticmap"
data="/home/thais/Pharmaco/ancestry/genome/SABE/shapeit/hg38SNP"

# Split data by chromossome and extracting only our data from merged file
# Create VCF files, compress and index
for chr in {1..22}; do \
  plink \
    --bfile ${plink} \
    --chr ${chr} \
    --recode vcf --out ${plink}_chr${chr}_split;
  bgzip ${plink}_chr${chr}_split.vcf
  bcftools index ${plink}_chr${chr}_split.vcf.gz
done

## rm chr from reference hg38
for chr in {1..22}; do \
  zcat 1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz | sed 's/chr//g' > 1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel_2.vcf
  bgzip 1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel_2.vcf
  bcftools index 1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel_2.vcf.gz
done
### Parsing genetic_map into chrs
# for i in {1..22}; do
# awk '$1 ~ /^'${i}'$/' genetic_map_hg38_withX.txt > genetic_map_hg38_chr${i}.txt
# done

# Now, we run shapeit in phasing mode (software default). We recommend to parallelize across chromosomes by screen command in the terminal console (could be for chromosomes 1 to 11 and 12 to 22).

# Using Shapeit4
for chr in {1..5}; do \
   /bin/shapeit4-4.2.2/bin/shapeit4.2 \
    -I ${plink}_chr${chr}_split.vcf.gz \
    -H ${ref}.chr${chr}.filtered.SNV_INDEL_SV_phased_panel_2.vcf.gz \
    -M ${mapfile}/chr${chr}.b38.gmap.gz \
    --region ${chr} \
    -O ${data}.chr${chr}.phased.vcf.gz \
    -T 6 \
    --log  ${data}.chr${chr}.phased.log;
done

# Using Shapeit4
for chr in {6..12}; do \
   /bin/shapeit4-4.2.2/bin/shapeit4.2 \
    -I ${plink}_chr${chr}_split.vcf.gz \
    -H ${ref}.chr${chr}.filtered.SNV_INDEL_SV_phased_panel_2.vcf.gz \
    -M ${mapfile}/chr${chr}.b38.gmap.gz \
    --region ${chr} \
    -O ${data}.chr${chr}.phased.vcf.gz \
    -T 6 \
    --log  ${data}.chr${chr}.phased.log;
done

# Using Shapeit4
for chr in {13..18}; do \
   /bin/shapeit4-4.2.2/bin/shapeit4.2 \
    -I ${plink}_chr${chr}_split.vcf.gz \
    -H ${ref}.chr${chr}.filtered.SNV_INDEL_SV_phased_panel_2.vcf.gz \
    -M ${mapfile}/chr${chr}.b38.gmap.gz \
    --region ${chr} \
    -O ${data}.chr${chr}.phased.vcf.gz \
    -T 6 \
    --log  ${data}.chr${chr}.phased.log;
done

# Using Shapeit4
for chr in {19..22}; do \
   /bin/shapeit4-4.2.2/bin/shapeit4.2 \
    -I ${plink}_chr${chr}_split.vcf.gz \
    -H ${ref}.chr${chr}.filtered.SNV_INDEL_SV_phased_panel_2.vcf.gz \
    -M ${mapfile}/chr${chr}.b38.gmap.gz \
    --region ${chr} \
    -O ${data}.chr${chr}.phased.vcf.gz \
    -T 6 \
    --log  ${data}.chr${chr}.phased.log;
done

###################################
# 3 Local ancestry estimation usign RfMix2
###################################

################ Query and reference panel files
# Creating query file an reference panel: VCF with all chromossomes -- indexed with bcftools
bcftools concat ${data}.chr1.phased.vcf.gz ${data}.chr2.phased.vcf.gz ${data}.chr3.phased.vcf.gz ${data}.chr4.phased.vcf.gz ${data}.chr5.phased.vcf.gz ${data}.chr6.phased.vcf.gz ${data}.chr7.phased.vcf.gz ${data}.chr8.phased.vcf.gz ${data}.chr9.phased.vcf.gz ${data}.chr10.phased.vcf.gz ${data}.chr11.phased.vcf.gz ${data}.chr12.phased.vcf.gz ${data}.chr13.phased.vcf.gz ${data}.chr14.phased.vcf.gz ${data}.chr15.phased.vcf.gz ${data}.chr16.phased.vcf.gz ${data}.chr17.phased.vcf.gz ${data}.chr18.phased.vcf.gz ${data}.chr19.phased.vcf.gz ${data}.chr20.phased.vcf.gz ${data}.chr21.phased.vcf.gz ${data}.chr22.phased.vcf.gz -Oz -o allhg381kGPNatives_mergeChr.vcf.gz

MERGE="allhg381kGPNatives_mergeChr.vcf.gz"

#verify if all the chrms are present, and in order
zgrep -v "#" allhg381kGPNatives_mergeChr.vcf.gz | cut -f1 | uniq
#verify if all snps
bcftools view -H ${MERGE} | awk '{print $3}' | wc -l

# Creating sample list from merge Chr with all individuals
bcftools query -l ${MERGE} > samplesQueryRef.txt

# Creating sample and ref lists
grep '^0' samplesQueryRef.txt > samplesQuery.txt
grep -v '^0' samplesQueryRef.txt > samplesREF.txt

samplesQUERY=samplesQuery_BRA.txt
samplesREF=samplesREF.txt

###### Creating query VCF 
vcftools --gzvcf ${MERGE} --keep ${samplesQUERY} --recode --recode-INFO-all --out query_hg38_mergeChr

# Index
bgzip query_hg38_mergeChr.recode.vcf
bcftools index query_hg38_mergeChr.recode.vcf.gz

###### Creating sample file from reference panel: sample name \t population
# Using $sampleREF
awk -F'_' '{print $2 "\t" $1}' ${samplesREF} > samplefile_temp.txt

# Create sample file only with AFR, EUR and AMR (PERU - Nativos)
# Part did on excel

## Create samplefile with the correct names 
awk '{print $2"_"$1 "\t" $2}' samplefile_temp2.txt | sed 's/NAT/AMR/' > samplefile_temp3.txt

awk '{print $1}' samplefile_temp3.txt > samplesREF_2.txt

###### Creating reference panel VCF
samplesREF=samplesREF_2.txt
vcftools --gzvcf ${MERGE} --keep ${samplesREF} --recode --recode-INFO-all --out Ref_hg38_mergeChr

#Compress and index
bgzip Ref_hg38_mergeChr.recode.vcf
bcftools index Ref_hg38_mergeChr.recode.vcf.gz

################ Gene map file
# Creating map file: VCF with all chromossomes -- indexed with bcftools
zcat chr1.b38.gmap.gz chr2.b38.gmap.gz chr3.b38.gmap.gz chr4.b38.gmap.gz chr5.b38.gmap.gz chr6.b38.gmap.gz chr7.b38.gmap.gz chr8.b38.gmap.gz chr9.b38.gmap.gz chr10.b38.gmap.gz chr11.b38.gmap.gz chr12.b38.gmap.gz chr13.b38.gmap.gz chr14.b38.gmap.gz chr15.b38.gmap.gz chr16.b38.gmap.gz chr17.b38.gmap.gz chr18.b38.gmap.gz chr19.b38.gmap.gz chr20.b38.gmap.gz chr21.b38.gmap.gz chr22.b38.gmap.gz > hg38.gmap

#Chr \t pos \t cM
awk '{print $2 "\t" $1 "\t" $3}' hg38.gmap > hg38_2.gmap

#Do not compress gmap

################ RfMix

QUERY="/home/thais/Pharmaco/ancestry/genome/SABE/shapeit/query_hg38_mergeChr.recode.vcf.gz" #VCF with all chrs
#Ref with NAT from BRA
REF="/home/thais/Pharmaco/ancestry/array/BRA/shapeit/Ref_hg38_mergeChr.recode.vcf.gz" #VCF reference
MAP="/home/thais/Pharmaco/ancestry/array/BRA/shapeit/geneticmap/hg38_2.gmap" #genetic map
SAMPLE="/home/thais/Pharmaco/ancestry/array/BRA/shapeit/samplefile_temp3.txt" #sample map

## Run RfMix

### BRA
for chr in {1..22}
do
/home/thais/bin/rfmix/rfmix \
-f ${QUERY} \
-r ${REF} \
-m ${SAMPLE} \
-g ${MAP} \
-o /home/thais/Pharmaco/ancestry/genome/SABE/shapeit/RFMix/BRA/BRA_allhg38SABE_rfmix_${chr} \
--n-threads=12 \
--chromosome=${chr}
done

### BRA - RR
for chr in {1..22}
do
/home/thais/bin/rfmix/rfmix \
-f ${QUERY} \
-r ${REF} \
-m ${SAMPLE} \
-g ${MAP} \
-o /home/thais/Pharmaco/ancestry/genome/SABE/shapeit/RFMix/BRA_RR/BRA_allhg38SABE_rfmix_RR_${chr} \
--n-threads=12 \
--chromosome=${chr} \
--reanalyze-reference
done

#VCF with all chrs
#QUERY="/home/thais/Pharmaco/ancestry/array/BRA/shapeit/query_hg38_mergeChr.recode.vcf.gz" 
QUERY="/home/thais/Pharmaco/ancestry/genome/SABE/shapeit/query_hg38_mergeChr.recode.vcf.gz"
#Ref with NAT from PER
REF="/home/thais/Pharmaco/ancestry/array/PER/shapeit/Ref_hg38_mergeChr.recode.vcf.gz" #VCF reference
MAP="/home/thais/Pharmaco/ancestry/array/BRA/shapeit/geneticmap/hg38_2.gmap" #genetic map
SAMPLE="/home/thais/Pharmaco/ancestry/array/PER/shapeit/samplefile_temp3.txt" #sample map


### PER
for chr in {1..22}
do
/home/thais/bin/rfmix/rfmix \
-f ${QUERY} \
-r ${REF} \
-m ${SAMPLE} \
-g ${MAP} \
-o /home/thais/Pharmaco/ancestry/genome/SABE/shapeit/RFMix/PER/PER_allhg38SABE_rfmix_${chr} \
--n-threads=12 \
--chromosome=${chr}
done


### PER - RR
for chr in {1..22}
do
/home/thais/bin/rfmix/rfmix \
-f ${QUERY} \
-r ${REF} \
-m ${SAMPLE} \
-g ${MAP} \
-o /home/thais/Pharmaco/ancestry/genome/SABE/shapeit/RFMix/PER_RR/PER_allhg38SABE_rfmix_RR_${chr} \
--n-threads=12 \
--chromosome=${chr} \
--reanalyze-reference
done

for m in BRA PER
do
grep "#" ${m}_allhg38_rfmix_1.sis.tsv > header.sis.tsv
cat ${m}_allhg38_rfmix_*.sis.tsv | grep -v "#" > ${m}_allhg38_rfmix_total.sis.tsv
cat header.sis.tsv ${m}_allhg38_rfmix_total.sis.tsv | sort -k1 -n > ${m}_allhg38_rfmix_total_h.sis.tsv

awk 'NR==2' ${m}_allhg38_rfmix_1.fb.tsv > header.fb.tsv
cat ${m}_allhg38_rfmix_*.fb.tsv | grep -v "chromosome" | grep -v "#" > ${m}_allhg38_rfmix_total.fb.tsv
cat header.fb.tsv ${m}_allhg38_rfmix_total.fb.tsv | sort -k1 -n > ${m}_allhg38_rfmix_total_h.fb.tsv

grep "#" ${m}_allhg38_rfmix_1.msp.tsv > header.msp.tsv
cat ${m}_allhg38_rfmix_*.msp.tsv | grep -v "#" > ${m}_allhg38_rfmix_total.msp.tsv
cat header.msp.tsv ${m}_allhg38_rfmix_total.msp.tsv | sort -k1 -n > ${m}_allhg38_rfmix_total_h.msp.tsv

grep "#" ${m}_allhg38_rfmix_1.rfmix.Q > header.rfmix.Q
cat ${m}_allhg38_rfmix_*.rfmix.Q | grep -v "#" > ${m}_allhg38_rfmix_total.rfmix.Q
cat header.rfmix.Q ${m}_allhg38_rfmix_total.rfmix.Q | sort -k1 -n  > ${m}_allhg38_rfmix_total_h.rfmix.Q
done

### BRA_RR

grep "#" BRA_allhg38_rfmix_RR_1.sis.tsv > header.sis.tsv
cat BRA_allhg38_rfmix_RR_*.sis.tsv | grep -v "#" > BRA_allhg38_rfmix_RR_total.sis.tsv
cat header.sis.tsv BRA_allhg38_rfmix_RR_total.sis.tsv | sort -k1 -n > BRA_allhg38_rfmix_RR_total_h.sis.tsv

awk 'NR==2' BRA_allhg38_rfmix_RR_1.fb.tsv > header.fb.tsv
cat BRA_allhg38_rfmix_RR_*.fb.tsv | grep -v "chromosome" | grep -v "#" > BRA_allhg38_rfmix_RR_total.fb.tsv
cat header.fb.tsv BRA_allhg38_rfmix_RR_total.fb.tsv | sort -k1 -n > BRA_allhg38_rfmix_RR_total_h.fb.tsv

grep "#" BRA_allhg38_rfmix_RR_1.msp.tsv > header.msp.tsv
cat BRA_allhg38_rfmix_RR_*.msp.tsv | grep -v "#" > BRA_allhg38_rfmix_RR_total.msp.tsv
cat header.msp.tsv BRA_allhg38_rfmix_RR_total.msp.tsv | sort -k1 -n > BRA_allhg38_rfmix_RR_total_h.msp.tsv

grep "#" BRA_allhg38_rfmix_RR_1.rfmix.Q > header.rfmix.Q
cat BRA_allhg38_rfmix_RR_*.rfmix.Q | grep -v "#" > BRA_allhg38_rfmix_RR_total.rfmix.Q
cat header.rfmix.Q BRA_allhg38_rfmix_RR_total.rfmix.Q | sort -k1 -n  > BRA_allhg38_rfmix_RR_total_h.rfmix.Q



##PER
grep "#" PER_allhg38_rfmix_1.sis.tsv > header.sis.tsv
cat PER_allhg38_rfmix_*.sis.tsv | grep -v "#" > PER_allhg38_rfmix_total.sis.tsv
cat header.sis.tsv PER_allhg38_rfmix_total.sis.tsv | sort -k1 -n > PER_allhg38_rfmix_total_h.sis.tsv

awk 'NR==2' PER_allhg38_rfmix_1.fb.tsv > header.fb.tsv
cat PER_allhg38_rfmix_*.fb.tsv | grep -v "chromosome" | grep -v "#" > PER_allhg38_rfmix_total.fb.tsv
cat header.fb.tsv PER_allhg38_rfmix_total.fb.tsv | sort -k1 -n > PER_allhg38_rfmix_total_h.fb.tsv

grep "#" PER_allhg38_rfmix_1.msp.tsv > header.msp.tsv
cat PER_allhg38_rfmix_*.msp.tsv | grep -v "#" > PER_allhg38_rfmix_total.msp.tsv
cat header.msp.tsv PER_allhg38_rfmix_total.msp.tsv | sort -k1 -n > PER_allhg38_rfmix_total_h.msp.tsv

grep "#" PER_allhg38_rfmix_1.rfmix.Q > header.rfmix.Q
cat PER_allhg38_rfmix_*.rfmix.Q | grep -v "#" > PER_allhg38_rfmix_total.rfmix.Q
cat header.rfmix.Q PER_allhg38_rfmix_total.rfmix.Q | sort -k1 -n  > PER_allhg38_rfmix_total_h.rfmix.Q

### PER_RR

grep "#" PER_allhg38_rfmix_RR_1.sis.tsv > header.sis.tsv
cat PER_allhg38_rfmix_RR_*.sis.tsv | grep -v "#" > PER_allhg38_rfmix_RR_total.sis.tsv
cat header.sis.tsv PER_allhg38_rfmix_RR_total.sis.tsv | sort -k1 -n > PER_allhg38_rfmix_RR_total_h.sis.tsv

awk 'NR==2' PER_allhg38_rfmix_RR_1.fb.tsv > header.fb.tsv
cat PER_allhg38_rfmix_RR_*.fb.tsv | grep -v "chromosome" | grep -v "#" > PER_allhg38_rfmix_RR_total.fb.tsv
cat header.fb.tsv PER_allhg38_rfmix_RR_total.fb.tsv | sort -k1 -n > PER_allhg38_rfmix_RR_total_h.fb.tsv

grep "#" PER_allhg38_rfmix_RR_1.msp.tsv > header.msp.tsv
cat PER_allhg38_rfmix_RR_*.msp.tsv | grep -v "#" > PER_allhg38_rfmix_RR_total.msp.tsv
cat header.msp.tsv PER_allhg38_rfmix_RR_total.msp.tsv | sort -k1 -n > PER_allhg38_rfmix_RR_total_h.msp.tsv

grep "#" PER_allhg38_rfmix_RR_1.rfmix.Q > header.rfmix.Q
cat PER_allhg38_rfmix_RR_*.rfmix.Q | grep -v "#" > PER_allhg38_rfmix_RR_total.rfmix.Q
cat header.rfmix.Q PER_allhg38_rfmix_RR_total.rfmix.Q | sort -k1 -n  > PER_allhg38_rfmix_RR_total_h.rfmix.Q


