## Task 9  

This purpose of this task is for the applicant to perform differential expression analysis while controlling for age and gender with DESEQ2. Present the results in an R Shiny application having:
+ An ordered table with the results.
+ And a way to let the user pick any gene and have the system present a chart with the expression of cases and controls.  


### Input files:  
mayo.path_aging.con.salmon.gene.counts.txt  
mayo.path_aging.con.phenotype.csv  

NOTE: I needed to rename some of the genes in the counts data which were converted to 1-Mar, 2-Mar.. ect. due to opening/saving the provided count data file in excel.  

### Tools required:  
R  

### R packages required:  
DESeq2 (1.24.0)  
shiny (1.5.0)  

### Output files:  
diffExp_DESeq2_complete.txt  
ShinyApp_Input.RData  

### Execution:  
`Rscript ./DESeq2_DEA_script.R`  

INSTRUCTION: Please copy the above-mentioned input files in the current directory (Task9_DESeq2) before launching the script.  

NOTE: If you run this R script interactively then upon runing the shiny serve (line 102) it will automatically open a new window in your browser having the presentable/required results (DE table and boxplot of gene-of-interest). If you run it via `Rscript`, then you will need to copy the server address into your browser window manually.  

### R session_info  

`sessioninfo::session_info()`  

`─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 setting  value                       
 version  R version 3.6.1 (2019-07-05)
 os       macOS Catalina 10.15.7      
 system   x86_64, darwin15.6.0        
 ui       RStudio                     
 language (EN)                        
 collate  en_US.UTF-8                 
 ctype    en_US.UTF-8                 
 tz       Europe/Luxembourg           
 date     2020-10-16                  

`─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 package              * version    date       lib source                            
 acepack                1.4.1      2016-10-29 [1] CRAN (R 3.6.0)                    
 annotate               1.62.0     2019-05-02 [1] Bioconductor                      
 AnnotationDbi          1.46.1     2019-08-20 [1] Bioconductor                      
 assertthat             0.2.1      2019-03-21 [1] CRAN (R 3.6.0)                    
 backports              1.1.5      2019-10-02 [1] CRAN (R 3.6.0)                    
 base64enc              0.1-3      2015-07-28 [1] CRAN (R 3.6.0)                    
 Biobase              * 2.44.0     2019-05-02 [1] Bioconductor                      
 BiocGenerics         * 0.30.0     2019-05-02 [1] Bioconductor                      
 BiocParallel         * 1.18.1     2019-08-06 [1] Bioconductor                      
 bit                    1.1-14     2018-05-29 [1] CRAN (R 3.6.0)                    
 bit64                  0.9-7      2017-05-08 [1] CRAN (R 3.6.0)                    
 bitops                 1.0-6      2013-08-17 [1] CRAN (R 3.6.0)                    
 blob                   1.2.0      2019-07-09 [1] CRAN (R 3.6.0)                    
 checkmate              1.9.4      2019-07-04 [1] CRAN (R 3.6.0)                    
 cli                    2.0.2      2020-02-28 [1] CRAN (R 3.6.0)                    
 cluster                2.1.0      2019-06-19 [1] CRAN (R 3.6.1)                    
 colorspace             1.4-1      2019-03-18 [1] CRAN (R 3.6.0)                    
 crayon                 1.3.4      2017-09-16 [1] CRAN (R 3.6.0)                    
 data.table             1.12.8     2019-12-09 [1] CRAN (R 3.6.0)                    
 DBI                    1.0.0      2018-05-02 [1] CRAN (R 3.6.0)                    
 DelayedArray         * 0.10.0     2019-05-02 [1] Bioconductor                      
 DESeq2               * 1.24.0     2019-05-02 [1] Bioconductor                      
 digest                 0.6.25     2020-02-23 [1] CRAN (R 3.6.0)                    
 dplyr                  0.8.3      2019-07-04 [1] CRAN (R 3.6.0)                    
 fansi                  0.4.1      2020-01-08 [1] CRAN (R 3.6.0)                    
 fastmap                1.0.1      2019-10-08 [1] CRAN (R 3.6.0)                    
 foreign                0.8-72     2019-08-02 [1] CRAN (R 3.6.0)                    
 Formula                1.2-3      2018-05-03 [1] CRAN (R 3.6.0)                    
 genefilter             1.66.0     2019-05-02 [1] Bioconductor                      
 geneplotter            1.62.0     2019-05-02 [1] Bioconductor                      
 GenomeInfoDb         * 1.20.0     2019-05-02 [1] Bioconductor                      
 GenomeInfoDbData       1.2.1      2019-09-10 [1] Bioconductor                      
 GenomicRanges        * 1.36.1     2019-09-06 [1] Bioconductor                      
 ggplot2                3.3.0.9000 2020-05-14 [1] Github (tidyverse/ggplot2@f1422ea)
 glue                   1.4.1      2020-05-13 [1] CRAN (R 3.6.1)                    
 gridExtra              2.3        2017-09-09 [1] CRAN (R 3.6.0)                    
 gtable                 0.3.0      2019-03-25 [1] CRAN (R 3.6.0)                    
 Hmisc                  4.2-0      2019-01-26 [1] CRAN (R 3.6.0)                    
 htmlTable              1.13.1     2019-01-07 [1] CRAN (R 3.6.0)                    
 htmltools              0.5.0      2020-06-16 [1] CRAN (R 3.6.2)                    
 htmlwidgets            1.5.2      2020-10-03 [1] CRAN (R 3.6.2)                    
 httpuv                 1.5.4      2020-06-06 [1] CRAN (R 3.6.2)                    
 IRanges              * 2.18.2     2019-08-24 [1] Bioconductor                      
 knitr                  1.24       2019-08-08 [1] CRAN (R 3.6.0)                    
 later                  1.1.0.1    2020-06-05 [1] CRAN (R 3.6.2)                    
 lattice                0.20-38    2018-11-04 [1] CRAN (R 3.6.1)                    
 latticeExtra           0.6-28     2016-02-09 [1] CRAN (R 3.6.0)                    
 locfit                 1.5-9.4    2020-03-25 [1] CRAN (R 3.6.0)                    
 magrittr               1.5        2014-11-22 [1] CRAN (R 3.6.0)                    
 Matrix                 1.2-17     2019-03-22 [1] CRAN (R 3.6.1)                    
 matrixStats          * 0.55.0     2019-09-07 [1] CRAN (R 3.6.0)                    
 memoise                1.1.0      2017-04-21 [1] CRAN (R 3.6.0)                    
 mime                   0.7        2019-06-11 [1] CRAN (R 3.6.0)                    
 munsell                0.5.0      2018-06-12 [1] CRAN (R 3.6.0)                    
 nnet                   7.3-12     2016-02-02 [1] CRAN (R 3.6.1)                    
 pillar                 1.4.2      2019-06-29 [1] CRAN (R 3.6.0)                    
 pkgconfig              2.0.2      2018-08-16 [1] CRAN (R 3.6.0)                    
 promises               1.1.1      2020-06-09 [1] CRAN (R 3.6.2)                    
 purrr                  0.3.3      2019-10-18 [1] CRAN (R 3.6.0)                    
 R6                     2.4.1      2019-11-12 [1] CRAN (R 3.6.0)                    
 RColorBrewer           1.1-2      2014-12-07 [1] CRAN (R 3.6.0)                    
 Rcpp                   1.0.5      2020-07-06 [1] CRAN (R 3.6.2)                    
 RCurl                  1.95-4.12  2019-03-04 [1] CRAN (R 3.6.0)                    
 rlang                  0.4.5      2020-03-01 [1] CRAN (R 3.6.0)                    
 rpart                  4.1-15     2019-04-12 [1] CRAN (R 3.6.1)                    
 RSQLite                2.1.2      2019-07-24 [1] CRAN (R 3.6.0)                    
 rstudioapi             0.11       2020-02-07 [1] CRAN (R 3.6.0)                    
 S4Vectors            * 0.22.0     2019-05-02 [1] Bioconductor                      
 scales                 1.0.0      2018-08-09 [1] CRAN (R 3.6.0)                    
 sessioninfo            1.1.1      2018-11-05 [1] CRAN (R 3.6.0)                    
 shiny                * 1.5.0      2020-06-23 [1] CRAN (R 3.6.2)                    
 stringi                1.4.6      2020-02-17 [1] CRAN (R 3.6.0)                    
 stringr                1.4.0      2019-02-10 [1] CRAN (R 3.6.0)                    
 SummarizedExperiment * 1.14.1     2019-07-31 [1] Bioconductor                      
 survival               2.44-1.1   2019-04-01 [1] CRAN (R 3.6.1)                    
 tibble                 2.1.3      2019-06-06 [1] CRAN (R 3.6.0)                    
 tidyselect             0.2.5      2018-10-11 [1] CRAN (R 3.6.0)                    
 vctrs                  0.3.0      2020-05-11 [1] CRAN (R 3.6.2)                    
 withr                  2.1.2      2018-03-15 [1] CRAN (R 3.6.0)                    
 xfun                   0.9        2019-08-21 [1] CRAN (R 3.6.0)                    
 XML                    3.98-1.20  2019-06-06 [1] CRAN (R 3.6.0)                    
 xtable                 1.8-4      2019-04-21 [1] CRAN (R 3.6.0)                    
 XVector                0.24.0     2019-05-02 [1] Bioconductor                      
 zlibbioc               1.30.0     2019-05-02 [1] Bioconductor                      
`