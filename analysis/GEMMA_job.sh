git clone https://github.com/rcc-uchicago/genetic-data-analysis-2.git

wget https://github.com/genetics-statistics/GEMMA/archive/refs/tags/v0.98.5.tar.gz

tar -xzf v0.98.5.tar.gz # will create GEMMA directory

conda create -n mtnuc
conda activate mtnuc

conda install -c bioconda gemma
conda install plink

# Check that GEMMA is working:
gemma -h 

# https://github.com/genetics-statistics/GEMMA
# https://github.com/genetics-statistics/GEMMA/blob/master/doc/manual.pdf

# Generate bed file and run PCA on vcf
plink --vcf cz_snps.vcf.gz --pca --make-bed --out cz_snps --allow-extra-chr --const-fid # Total genotyping rate is 0.878658. 26988597 variants and 30 people pass filters and QC.


plink --file mydata --pheno pheno.raw --assoc --maf 0.05 --out run1


# Estimate relatedness matrix
gemma -bfile JuncoForGWAS_Feb2025 -gk 1 -o JuncoRelate # -bfile is prefix for plink .ped formatted files, which also includes the phenotype info, -gk is the relatedness matrix it can be either centered gk 1 or standardized gk2, -o output file name
# Going to use the IBGS matrix from ANGSD for relatedness

# Perform GWAS
# -k is the path to the relatedness matrix, -lmm specifies which kind of mixed effect model (4 runs all 3 kinds) and -o is the output. 
gemma -bfile JuncoForGWAS_Feb2025 -k output/JuncoRelate.sXX.txt -lmm 4 -o Junco_BasicGemma6Feb25

gemma -bfile XXX -k czsnps.ibsMat -lmm 4 -o cz_snps_gemma

gemma -p pheno.txt -c covar.txt -a map.txt -g geno.txt -notsnp -lm 2 -outdir . -o tibia


wget https://www.well.ox.ac.uk/~gav/resources/qctool_v2.0.8-osx.tgz
tar -zxvf qctool_v2.0.8-osx.tgz
cd qctool_v2.0.8-osx
# Check that it works:
./qctool -help

./qctool -g ../../data/cz_snps.vcf.gz -ofiletype bimbam_dosage -og ../../data/cz_snps.txt

### GEMMA INPUT FILES
# GEMMA requires four main input files:
#		1) Genotypes > geno.txt
#		2) SNP annotations > map.txt
#		2) Phenotypes/traits -- MITOTYPE
#		3) Relatedness matrix -- ANGSD IBS MATRIX
#   	4) Covariates (optional) -- POPULATION STRUCTURE

### Following the GitHub tutorial:

