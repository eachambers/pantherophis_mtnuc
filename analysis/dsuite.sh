### Use Dsuite package to run ABBA-BABA analysis and calculate Fstats  (using SNPs called with FreeBayes on high coverage sequencing dataset)

## create input files for analyses

# create sets.txt file:
PanObs	Outgroup
H16057	guttatus
H21189	slowinskii
O42736	emoryi 
TJH3395	xxx

# create newick tree file:
(((emoryi,slowinskii),guttatus),Outgroup);

# create test_trios.txt:
emoryi	slowinskii	guttatus

## run ABBA-BABA analysis on genome-wide SNPs
Dsuite Dtrios HCSamples_PanObs_fb_snps_filtered.recode.vcf.gz sets.txt

## calculate fbranch statistic 
Dsuite Fbranch treefile.nwk sets_tree.txt > fbranch_matrix.txt

## calculate Fstats to identify regions of the genome with signatures of introgression
Dsuite Dinvestigate HCSamples_PanObs_fb_snps_filtered.recode.vcf.gz sets.txt test_trios.txt


## run ABBA-BABA analysis on N-mt SNPs only (run in separate directory to avoid output file confusion)

# use N-mt gene coordinates (saved in nmtcoords.txt) to make a vcf containing only N-mt SNPs
bcftools view -R nmtcoords.txt -Oz -o HCSamples_PanObs_fb_nmt_snps_filtered.recode.vcf.gz HCSamples_PanObs_fb_snps_filtered.recode.vcf.gz

#run ABBA-BABA analysis with nmt SNPs only
Dsuite Dtrios HCSamples_PanObs_fb_nmt_snps_filtered.recode.vcf.gz sets.txt -t treefile.nwk

#calculate fbranch statistic for nmt dataset
Dsuite Fbranch treefile.nwk sets_tree.txt > fbranch_matrix_nmt.txt


## run ABBA-BABA analysis on control SNPs only

# use control gene coordinates (saved in control_coords.txt) to make a vcf containing only control SNPs
bcftools view -R control_coords.txt -Oz -o HCSamples_PanObs_fb_cont_snps_filtered.recode.vcf.gz HCSamples_PanObs_fb_snps_filtered.recode.vcf.gz 

# run ABBA-BABA analysis with control SNPs only
Dsuite Dtrios HCSamples_PanObs_fb_cont_snps_filtered.recode.vcf.gz sets.txt -t treefile.nwk

#calculate fbranch statistic for control dataset
/work/05104/thomm/ls6/Dsuite/Build/Dsuite Fbranch treefile.nwk sets_tree.txt > fbranch_matrix_cont.txt


## subsample SNPs from control and genome-wide datasets to get an equal number as N-mt dataset (so that p-values are comparable)

# randomly subsample control SNPs to obtain a number comparable to N-mt SNPs
bcftools view HCSamples_PanObs_fb_cont_snps_filtered.recode.vcf.gz | vcfrandomsample -r 0.23215 > cont_snps_subset.vcf

# randomly subsample genome-wide SNPs to obtain a number comparable to N-mt SNPs
bcftools view HCSamples_PanObs_fb_snps_filtered.recode.vcf.gz | vcfrandomsample -r 0.00149 > all_snps_subset.vcf

# run ABBA-BABA analysis on subsampled control SNPs
Dsuite Dtrios cont_snps_subset.vcf sets.txt

# run ABBA-BABA analysis on subsampled genome-wide SNPs
Dsuite Dtrios all_snps_subset.vcf sets.txt