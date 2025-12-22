# phosphoproteomics_src
Scripts used for SRC phosphoproteomic analysis in HEK293T cells as reported in the manuscript. Please cite:
Jain, P. et al. (2025) “Minimal Perturbation of Activation Loop Dynamics Rewires Kinase Signaling,” bioRxiv, p. 2025.10.15.682502. doi: 10.1101/2025.10.15.682502

The R Scripts are numbered in the order of execution in the folder "scripts", though some parts may be redundant to add cosmetic improvements. The folder 'data' in this repository includes the data files used to start the first script "1_LFQpeptideFragPipeAnalystR_filt-vsn-bpca_or_man_251014.R". The scripts that follow the first script will then use these data files and the output files to run. 

IMPORTANT: The directory and initial input names in the scripts are given for local settings in Ubuntu 22.04 and 24.04 LTR. The directory may need to be updated for your specific computing environment. 

R version 4.5.1 was used with the following packages:

              Package Version
        AnnotationDbi  1.70.0
      clusterProfiler  4.16.0
       ComplexHeatmap  2.24.1
                dplyr   1.1.4
           enrichplot  1.28.2
     FragPipeAnalystR   1.1.0
              ggplot2   4.0.0
              ggrepel   0.9.6
               magick   2.9.0
          matrixStats   1.5.0
              msigdbr  25.1.1
         org.Hs.eg.db  3.21.0
            patchwork   1.3.2
         RColorBrewer   1.1.3
              stringr   1.5.2
    SummarizedExperiment  1.38.1
               tibble   3.3.0
                tidyr   1.3.1
            tidyverse   2.0.0

Please note that A0-B2 scripts are not used to generate results shown in the manuscript.

