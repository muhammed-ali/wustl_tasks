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
