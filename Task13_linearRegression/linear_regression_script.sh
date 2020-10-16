#!/bin/bash -l

./plink --bfile CHR1_Task_Renamed --pheno Phenotype.txt --pheno-name Phenotype --covar Phenotype.txt --covar-name age, gender, PC1, PC2 --out LR --linear --ci 0.95

echo "PLINK linear regression computed, launching R script now for getting the top SNP: "

Rscript ./linear_regression_script.R