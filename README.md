# phosphoproteomics_src

This repository contains the scripts used for SRC phosphoproteomic analysis in HEK293T cells. If you use this repository, please cite:

Jain, P. et al. (2025), “Minimal Perturbation of Activation Loop Dynamics Rewires Kinase Signaling,” bioRxiv, p. 2025.10.15.682502. doi: 10.1101/2025.10.15.682502

## Scripts and Data

The R scripts are numbered in the order of execution within the "scripts" folder. Note that some parts may include redundancies for additional cosmetic improvements. The "data" folder in this repository includes data files required for analysis.

**IMPORTANT:** The directory paths and initial input filenames used in the scripts are configured for local settings on Ubuntu 22.04 and 24.04 LTS. Ensure that you update these paths as needed for your specific computing environment.

## Software and Dependencies

The analysis was conducted using R version 4.5.1 with the following packages:

| Package                | Version |  
|------------------------|---------|  
| AnnotationDbi          | 1.70.0  |  
| clusterProfiler        | 4.16.0  |  
| ComplexHeatmap         | 2.24.1  |  
| dplyr                  | 1.1.4   |  
| enrichplot             | 1.28.2  |  
| FragPipeAnalystR       | 1.1.0   |  
| ggplot2                | 4.0.0   |  
| ggrepel                | 0.9.6   |  
| magick                 | 2.9.0   |  
| matrixStats            | 1.5.0   |  
| msigdbr                | 25.1.1  |  
| org.Hs.eg.db           | 3.21.0  |  
| patchwork              | 1.3.2   |  
| RColorBrewer           | 1.1.3   |  
| stringr                | 1.5.2   |  
| SummarizedExperiment   | 1.38.1  |  
| tibble                 | 3.3.0   |  
| tidyr                  | 1.3.1   |  
| tidyverse              | 2.0.0   |  

## Additional Notes

Please be aware that the scripts labeled **A0-B2** were not used to generate the results shown in the manuscript.

