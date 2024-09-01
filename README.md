# Introgression and coevolution in the *Pantherophis guttatus* complex

The following is code to run analyses and generate figures from [Marshall et al. (XXX)](REFER). Raw data files are provided on Dryad [here](https://datadryad.org/XXX). The structure of this repository is that scripts are divided into those that [perform analyses](https://github.com/eachambers/pantherophis_mtnuc/tree/main/analysis) and those that create [data visualizations](https://github.com/eachambers/pantherophis_mtnuc/tree/main/data_viz).

## Scripts in repository

- **Bioinformatics pipeline:**
  - Trimming, mapping, and calling variants for high coverage samples with freebayes
  - Extracting random 5-kb blocks for species tree analysis; extracting N-mt and control gene datasets
  - Extracting mitogenomes with Mitofinder
  - Call variants for contact zone (i.e., lower coverage) samples using ANGSD
  - Download outgroup sequences from SRA
  
- **Performing phylogenomic, mitochondrial and gene tree analyses:**
  - Estimating species tree with Starbeast3
  - Estimating mitochondrial tree
  - Estimating N-mt and control gene trees

- **Examining topologies of gene trees:**
    - Calculating Robinson-Foulds distances among gene trees ([functions](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/RFdists.R) and [analysis script](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/RFdists_analysis.R))
    - Visualizing results ([Fig. 2C](https://github.com/eachambers/pantherophis_mtnuc/blob/main/data_viz/RFdists_figure.R))

- **GWAS analysis:**
    - Running the GWAS itself ([script](REFER))
    - Processing GWAS results and getting summary statistics ([functions](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/GWAS.R) and [analysis script](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/GWAS_analysis.R))
    - Visualizing results ([Fig. 3A & S4 & S5](https://github.com/eachambers/pantherophis_mtnuc/blob/main/data_viz/GWAS_figures.R))

- **Diagnostic differences:**
    - Calculating diagnostic differences and running DAPC ([functions](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/Diagnosticdiffs.R) and [analysis script](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/Diagnosticdiffs_analysis.R))
    - Visualizing diagnostic differences results ([Figs. 3C & 3D](https://github.com/eachambers/pantherophis_mtnuc/blob/main/data_viz/Diagnosticdiffs_figures.R))
    - Visualizing DAPC results ([Fig. S6](https://github.com/eachambers/pantherophis_mtnuc/blob/main/data_viz/Diagnosticdiffs_figures.R))
 
- **ABBA-BABA analysis:**
    - Running ABBA-BABA analyses with Dtrios and Fstat analysis with Dinvestigate
    - Calculating mean fdM in sliding windows from fstat results [analysis script]
    - (https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/ABBABABA.R))
    - Visualizing results ([Figs. 2B & S3](https://github.com/eachambers/pantherophis_mtnuc/blob/main/data_viz/ABBABABA_figure.R))
