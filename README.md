# pantherophis_mtnuc

The following is code to run analyses and generate figures from [Marshall et al. (XXX)](REFER). Raw data files are provided on Dryad [here](https://datadryad.org/XXX). The structure of this repository is that scripts are divided into those that perform analyses ([here](https://github.com/eachambers/pantherophis_mtnuc/tree/main/analysis)) and data visualizations ([here](https://github.com/eachambers/pantherophis_mtnuc/tree/main/data_viz)).

## Scripts in repository

- **Examining topologies of gene trees:**
    - Calculating Robinson-Foulds distances among gene trees ([functions](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/RFdists.R) and [analysis script](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/RFdists_analysis.R))
    - Visualizing results ([Fig. 2C](https://github.com/eachambers/pantherophis_mtnuc/blob/main/data_viz/RFdists_figure.R))

- **GWAS analysis:**
    - Running the GWAS itself ([script](REFER))
    - Processing GWAS results and getting summary statistics ([functions](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/GWAS.R) and [analysis script](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/GWAS_analysis.R))
    - Visualizing results ([Fig. 3A & S3](https://github.com/eachambers/pantherophis_mtnuc/blob/main/data_viz/GWAS_figures.R))

- **Diagnostic differences:**
    - Calculating diagnostic differences ([functions](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/Diagnosticdiffs.R) and [analysis script](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/Diagnosticdiffs_analysis.R))
    - Visualizing results ([Figs. 3C & 3D](https://github.com/eachambers/pantherophis_mtnuc/blob/main/data_viz/Diagnosticdiffs_figures.R))
 
- **ABBA-BABA analysis:**
    - Calculating mean fdM in sliding windows from ABBA-BABA results ([analysis script](https://github.com/eachambers/pantherophis_mtnuc/blob/main/analysis/ABBABABA.R))
    - Visualizing results ([Figs. 2B & S2](https://github.com/eachambers/pantherophis_mtnuc/blob/main/data_viz/ABBABABA_figure.R))
