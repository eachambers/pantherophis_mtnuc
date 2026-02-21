# Introgression and coevolution in the *Pantherophis guttatus* complex

The following is code to run analyses and generate figures from [Marshall et al. (XXX)](REFER). Raw data files are provided on Dryad [here](https://datadryad.org/XXX). The structure of this repository is that scripts are divided into those that [perform analyses](https://github.com/eachambers/pantherophis_mtnuc/tree/main/analysis) and those that create [data visualizations](https://github.com/eachambers/pantherophis_mtnuc/tree/main/data_viz).

## Scripts in repository

**1. Bioinformatics pipeline:**
  - [Bioinformatics walkthrough](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/1a_bioinformatics_pipeline.sh) contains code for: trimming, mapping, calling variants for high coverage samples with freebayes, and calling variants for lower coverage samples with ANGSD
  - [Extracting alignments](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/1b_extract_alignments.sh) extracts random 5-kb blocks for species tree analysis; extracts N-mt and control gene datasets
  - [Extracting mitogenomes](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/1c_mitofinder.sh) extracts mitogenomes with MitoFinder
  - [Downloading outgroup sequences](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/1d_outgroup_sra_download.sh) from SRA
  - [Principal coordinates analysis](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/1e_PCoA.R) and data visualization (Fig. S2)
  
**2. Performing phylogenomic, mitochondrial and gene tree analyses:**
  - Estimating species tree with Starbeast3
  - Estimating mitochondrial tree
  - Estimating N-mt and control gene trees

**3. ABBA-BABA analysis:**
  - Running ABBA-BABA analyses with Dtrios and Fstat analysis with Dinvestigate [script](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/3a_dsuite.sh)
  - Calculating mean fdM in sliding windows from fstat results [analysis script](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/3b_ABBABABA.R)
  - Visualizing results ([Figs. 2B & S3](https://github.com/eachambers/pantherophis_mtnuc/blob/main/data_viz/ABBABABA_figure.R))

**4. Examining topologies of gene trees:**
  - Calculating Robinson-Foulds distances among gene trees ([functions](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/4a_RFdists.R) and [analysis script](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/4b_RFdists_analysis.R))
  - Visualizing results ([Fig. 2C](https://github.com/eachambers/pantherophis_mtnuc/blob/main/data_viz/RFdists_figure.R))

**5. Association analysis:**
  - [Walkthrough](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/5a_Association_analysis_walkthrough.txt) for running the association analysis
  - Running the GWAS itself ([script](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/5b_RDA_GWAS_marshall_2022.R))
  - Processing GWAS results and getting summary statistics ([analysis script](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/5c_GWAS_processing.R))
  - Visualizing results ([Figs. 3A & S4 & S5](https://github.com/eachambers/pantherophis_mtnuc/blob/main/data_viz/GWAS_figures.R))
  - Calculating SNP scores based on PCA on non-LD pruned data [script](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/5d_SNP_scores.R)

**6. Diagnostic differences:**
  - Performing LD-pruning on input data [script](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/6a_ld_pruning.sh)
  - Calculating diagnostic differences and running DAPC ([functions](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/6b_Diagnosticdiffs.R) and [analysis script](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/6c_Diagnosticdiffs_analysis.R))
  - Visualizing diagnostic differences and DAPC results ([Figs. 3C, 3D, and S6](https://github.com/eachambers/pantherophis_mtnuc/blob/main/data_viz/Diagnosticdiffs_figures.R))
