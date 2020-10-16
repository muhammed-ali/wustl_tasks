## Task 13  

Perform a linear regression adjusted by age, gender and the first two principal components. Report lowest p value for the analyses with the effect and the 95% confident intervals.  


### Input files:  
CHR1_Task_Renamed.bim  
CHR1_Task_Renamed.bed  
CHR1_Task_Renamed.fam  

Phenotype.txt  

### Tools required:  
Plink  
R  

### Output files:  
LR.assoc.linear  
LR.log  

### Execution:  
`sh ./linear_regression_script.sh`  

INSTRUCTION: Please copy the above-mentioned input files in the current directory (Task13_linearRegression) before launching the script.  

NOTE: At the end of execution the program reports the lowest p value for the analyses with the effect and the 95% confident intervals.  