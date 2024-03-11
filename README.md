# README
## Article Information
This repository provides access to the data and source code used for the manuscript    
### **The strength of sexual signals predicts same-sex paring in termites**  
<!--Author names are commented out for DBR.-->
#### **Nobuaki Mizumoto, Sang-Bin Lee, Thomas Chouvenc**  
  
Preprint will be available at bioRxiv. [![DOI:10.1101/2024.03.07.583902](http://img.shields.io/badge/10.1101/2024.03.07.583902.svg)]  
The all data will be uploaded in Zenodo upon acceptance: [![DOI](https://zenodo.org/badge/DOI/XXXDOIXXX.svg)](https://doi.org/XXXDOIXXX)
  
This study compared the same-sex tandem running behavior between two termtie species (Coptotermes formosanus and Coptotermes gestroi) that use the same chemicals for tandem runs but have it in different quantities.  
We recorded movement patterns of termite pairs, either in combination of Female-Male, Female-Female, Male-Male, by using a tracking software, [UMATracker](https://ymnk13.github.io/UMATracker/).  
Also, measured the body length of termites and the size of arena, using a [python program](./analysis/code/video_scale_BL.py). All these raw data were stored at [data_raw](./analysis/data_raw).
First, run the [format_trajectories.R](./analysis/code/format_trajectories.R) to format data for further statistical analysis and data visualization. Formatted data will be stored at [data_fmt](./analysis/data_fmt).  
Then, run the [output.R](./analysis/code/output.R) to obtain all outputs.  

## Table of Contents
This repository includes tracking data, R codes to analyze it, and Python code for video analysis.  
* [README](./README.md)
* [analysis](./analysis)
  * [code](./analysis/code)
    * [format_trajectories.R](./analysis/code/format_trajectories.R)
    * [output.R](./analysis/code/output.R)
    * [video_scale_BL.py](./analysis/code/video_scale_BL.py)
  * [data_raw](./analysis/data_raw) - folder containing raw data
  * [data_fmt](./analysis/data_fmt) - folder containing data converted from raw data
  * [output](./analysis/output) - folder containing outputs
* [draft](./draft)

## Session information
```
R version 4.3.1 (2023-06-16 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

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
 [1] dplyr_1.1.4       stringr_1.5.0     data.table_1.14.8 rstatix_0.7.2     multcomp_1.4-23  
 [6] TH.data_1.1-2     MASS_7.3-60       mvtnorm_1.1-3     car_3.1-2         carData_3.0-5    
[11] lme4_1.1-33       Matrix_1.5-4.1    coxme_2.2-18.1    survival_3.5-5    magick_2.8.3     
[16] patchwork_1.2.0   Hmisc_5.1-1       viridis_0.6.2     viridisLite_0.4.1 survminer_0.4.9  
[21] ggpubr_0.6.0      ggplot2_3.4.4     tidyr_1.3.0       bdsmatrix_1.3-6  

loaded via a namespace (and not attached):
 [1] gtable_0.3.3      xfun_0.39         htmlwidgets_1.6.2 lattice_0.21-8    vctrs_0.6.4      
 [6] tools_4.3.1       generics_0.1.3    sandwich_3.0-2    tibble_3.2.1      fansi_1.0.4      
[11] cluster_2.1.4     pkgconfig_2.0.3   checkmate_2.3.0   lifecycle_1.0.3   compiler_4.3.1   
[16] munsell_0.5.0     codetools_0.2-19  htmltools_0.5.5   yaml_2.3.7        htmlTable_2.4.2  
[21] Formula_1.2-5     nloptr_2.0.3      pillar_1.9.0      boot_1.3-28.1     rpart_4.1.19     
[26] abind_1.4-5       nlme_3.1-162      km.ci_0.5-6       tidyselect_1.2.0  digest_0.6.31    
[31] stringi_1.7.12    purrr_1.0.1       splines_4.3.1     fastmap_1.1.1     grid_4.3.1       
[36] colorspace_2.1-0  cli_3.6.1         magrittr_2.0.3    base64enc_0.1-3   utf8_1.2.3       
[41] broom_1.0.4       foreign_0.8-84    withr_2.5.0       scales_1.2.1      backports_1.4.1  
[46] rmarkdown_2.22    nnet_7.3-19       gridExtra_2.3     ggsignif_0.6.4    zoo_1.8-12       
[51] evaluate_0.20     knitr_1.42        KMsurv_0.1-5      survMisc_0.5.6    rlang_1.1.0      
[56] Rcpp_1.0.10       xtable_1.8-4      glue_1.6.2        minqa_1.2.5       rstudioapi_0.14  
[61] R6_2.5.1       
```
