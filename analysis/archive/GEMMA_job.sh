### ======================== RUNNING GEMMA GWAS ANALYSIS ========================

mamba create -n mtnuc
mamba activate mtnuc

# Install relevant software
mamba install -c bioconda gemma
mamba install plink

# Check that GEMMA is working:
gemma -h # v.0.98.5

# https://github.com/genetics-statistics/GEMMA
# https://github.com/genetics-statistics/GEMMA/blob/master/doc/manual.pdf
# https://adegenet.r-forge.r-project.org/files/Glasgow2015/practical-GWAS_day4.pdf
# Following the GitHub tutorial: https://github.com/rcc-uchicago/genetic-data-analysis-2.git

# ===================================== INPUT FILES =====================================

### GEMMA INPUT FILES
# GEMMA requires four main input files:
#		1) Genotypes, provided in Plink ped format (bed, bim, fam) including phenotype info > cz_snps
#		2) SNP annotations > map.txt [optional?]
#		2) Phenotypes/traits -- MITOTYPE
#		3) Relatedness matrix -- ANGSD IBS MATRIX
#   	4) Covariates (optional) -- POPULATION STRUCTURE

# How many SNPs are there in cz_snps vcf?
bcftools query -f '%POS\n' cz_snps.vcf.gz | wc -l # 26,988,597 sites

#		1. Generate bed files (which will be used as input into GEMMA) and run PCA for pop 
#		   structure correction:
plink --vcf cz_snps.vcf --pca --make-bed --out cz_snps --allow-extra-chr --const-fid

#		2. Add phenotypes (i.e., mitotypes) to the fam file. You should be able to do this
#		   by doing the following, but it didn't work for me:
# plink --vcf cz_snps.vcf.gz --pheno plink_pheno.txt --make-bed --out cz_snps --allow-extra-chr --const-fid
#		   Because it didn't work, I ended up just manually editing the last column of the fam
#		   file to match the phenotypes (i.e., mitotypes) for each individual.

#		3. Run the first two steps of `GEMMA_analysis.R` to generate the properly formatted
#		   GEMMA input files for the covariates and phenotype files.

# ===================================== PERFORM GWAS =====================================

# -k is the path to the relatedness matrix, -lmm specifies which kind of mixed effect model (4 runs all 3 kinds) and -o is the output. 
gemma -bfile cz_snps -k czsnps.ibsMat -c cz_snps_covar.txt -lmm 4 -o cz_snps_gemma

# ===================================== OUTPUT FILES =====================================

.cXX.txt.gz kinship matrix
.log.txt summary of association analysis
.assoc.txt full table of association results (p-values, effect size estimates, etc.)

# ========================== REMAINING QUESTIONS / THINGS TO TRY =========================

# Do we need to do LD-pruning on the data prior to running the GWAS? If so:
bcftools +prune -m 0.6 -w 50 cz_snps.vcf.gz -O v -o cz_snps_ldp.vcf
bcftools query -f '%POS\n' cz_snps_ldp.vcf | wc -l # now only 10,128,602 sites
plink --vcf cz_snps_ldp.vcf --pca --make-bed --out cz_snps_ldp --allow-extra-chr --const-fid

# Rather than using the ANGSD IBS matrix, you can also estimate relatedness matrix using GEMMA which is commonly done
# I am not sure if there's a benefit to doing so but it might be worth seeing how much the results change:
# -bfile is prefix for plink .ped formatted files, which also includes the phenotype info, -gk is the relatedness matrix it can be either centered gk 1 or standardized gk2, -o output file name
gemma -bfile cz_snps -gk 1 -o cz_snps_relate
# Then you should be able to run GEMMA the same way but with this new relatedness file, like so:
gemma -bfile cz_snps -k cz_snps_relate.txt -c cz_snps_covar.txt -lmm 4 -o cz_snps_gemma

# ================================== RUNNING THIS ON TACC ================================

mamba create -n mtnuc
mamba activate mtnuc
mamba install -c r r-essentials

