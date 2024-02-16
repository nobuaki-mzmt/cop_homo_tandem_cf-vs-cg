# README
## Article Information
This repository provides access to the data and source code used for the manuscript    
### **The strength of sexual signals predicts same-sex paring in termites**  
Author names are commented out for DBR.
Nobuaki Mizumoto, Sang-Bin Lee, Thomas Chouvenc  

Preprint will be available at bioRxiv. [![DOI:XXX](http://img.shields.io/badge/DOI-10.1101/XXX.svg)]  
The all data will be uploaded in Zenodo upon acceptance: [![DOI](https://zenodo.org/badge/DOI/XXXDOIXXX.svg)](https://doi.org/XXXDOIXXX)
  
This study compared the same-sex tandem running behavior between two termtie species (Coptotermes formosanus and Coptotermes gestroi) that use the same chemicals for tandem runs but have it in different quantities.  

This includes tracking data, R codes to analyze it, and Python code for video analysis.  

## Table of Contents
* [README](./README.md) - this file
* [code](./analysis/code) - folder containing scripts for the analysis
  * [format_trajectories.R](./analysis/code/format_trajectories.R)
  * [output.R](./analysis/code/output.R)
  * [video_scale_BL.py](./analysis/code/video_scale_BL.py)
* [data_raw](./analysis/data_raw) - folder containing raw data
* [data_fmt](./analysis/data_fmt) - folder containing data converted from raw data
* [output](./analysis/output) - folder containing outputs
  
## Session information
```
R version 4.3.1 (2023-06-16 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 22621)

Matrix products: default


locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

time zone: America/Chicago
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] rstatix_0.7.2     multcomp_1.4-25   TH.data_1.1-2     MASS_7.3-60      
 [5] mvtnorm_1.2-3     car_3.1-2         carData_3.0-5     lme4_1.1-34      
 [9] Matrix_1.6-1      coxme_2.2-18.1    survival_3.5-5    patchwork_1.2.0  
[13] Hmisc_5.1-1       survminer_0.4.9   ggpubr_0.6.0      tidyr_1.3.0      
[17] viridis_0.6.3     viridisLite_0.4.2 ggplot2_3.4.2     dplyr_1.1.2      
[21] stringr_1.5.0     data.table_1.14.8 bdsmatrix_1.3-6  

loaded via a namespace (and not attached):
 [1] gtable_0.3.3      xfun_0.41         htmlwidgets_1.6.2 lattice_0.21-8   
 [5] vctrs_0.6.3       tools_4.3.1       generics_0.1.3    sandwich_3.0-2   
 [9] tibble_3.2.1      fansi_1.0.4       cluster_2.1.4     pkgconfig_2.0.3  
[13] checkmate_2.3.1   lifecycle_1.0.3   compiler_4.3.1    munsell_0.5.0    
[17] codetools_0.2-19  htmltools_0.5.7   yaml_2.3.7        htmlTable_2.4.2  
[21] Formula_1.2-5     pillar_1.9.0      nloptr_2.0.3      rpart_4.1.19     
[25] boot_1.3-28.1     abind_1.4-5       nlme_3.1-162      km.ci_0.5-6      
[29] tidyselect_1.2.0  digest_0.6.33     stringi_1.7.12    purrr_1.0.2      
[33] splines_4.3.1     fastmap_1.1.1     grid_4.3.1        colorspace_2.1-0 
[37] cli_3.6.1         magrittr_2.0.3    base64enc_0.1-3   utf8_1.2.3       
[41] broom_1.0.5       foreign_0.8-84    withr_2.5.0       scales_1.2.1     
[45] backports_1.4.1   rmarkdown_2.25    nnet_7.3-19       gridExtra_2.3    
[49] ggsignif_0.6.4    zoo_1.8-12        evaluate_0.21     knitr_1.45       
[53] KMsurv_0.1-5      survMisc_0.5.6    rlang_1.1.1       Rcpp_1.0.10      
[57] xtable_1.8-4      glue_1.6.2        rstudioapi_0.15.0 minqa_1.2.6      
[61] R6_2.5.1        
```
