# IFNγ cancer project
R and Matlab code used for the analysis performed in the paper "Bystander IFNγ activity promotes widespread and sustained cytokine signaling altering the tumor microenvironment".

### Recommendations

Before running any of the *R* scripts please install the following R libraries :
- **Pagoda2** (https://github.com/hms-dbmi/pagoda2/tree/master/R)
- **fifer** (https://github.com/dustinfife/fifer)
- **uwot** (https://github.com/jlmelville/uwot)
- **igraph** (https://igraph.org/r/)
- **ggplot2** (https://ggplot2.tidyverse.org)
- **DESeq2** (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- **FactoMineR** (https://cran.r-project.org/web/packages/FactoMineR/index.html)
- **pheatmap** (https://cran.r-project.org/web/packages/pheatmap/index.html)

Concerning the *Matlab®* scripts, please use at least *Matlab* 2017 and have the **Image Processing Toolbox** available.

### List of the scripts
Each script corresponds to a specific part of the paper :

- **In_vitro_RNAseq_analysis.R** : R script used to analyze the RNAseq data generated from B cell lymphoma cell line stimulated by IFNγ during various time.
- **Melanoma_data_analysis.R** : R script used to re-analyze the melanoma MARS-seq data from Li et al. using Pagoda2 pipeline.
- **HNSCC_data_analysis.R** : R script use to re-analyze the Head and Neck Squamous Cell Carcinoma from Puram et al. using Pagoda2 pipeline.
- **Visualize_Image_analysis.R** : R script that gathers results of image analysis performed by *Matlab* and create the figures used in the paper.
- **Microscopy_Image_analysis_main.m** : Matlab script that performs co-localisation analysis between STAT1-GFP and nuclear mCherry signal both for in-vitro and intra-vital images shown in the main figures. A detailed description of the image processing pipeline is available in as a supplementary figure of the paper. 
- **Microscopy_Image_analysis_sup.m** : Matlab script that performs co-localisation analysis between STAT1-GFP and nuclear mCherry signal for intra-vital images shown in supplementary figures.
- **Kinetic_analysis.m**: Matlab script used to perform the co-localisation analysis of IFNγ in-vitro stimulated cells across different times.
- **LoadImage.m** : Basic Matlab script that loads multiple-stack .tiff files into Matlab.
- **Pre_processing.m** : Basic Matlab script used to clean images before the analysis (background removal, intensity adjustment).
- **Nuclei_identification.m** : Matlab script used to segment cells through image binarisation and watershed transform.
