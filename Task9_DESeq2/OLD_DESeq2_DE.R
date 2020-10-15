# NOTE: First rename some of the genes which have been converted to 1-Mar, 2-Mar.. ect. due to opening/saving the provided count data file in excel.
# Tutorial: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# ssh iris-cluster
# cd /home/users/mali/Tasks/WUSTL_Tasks/Task9_DESeq2
# srun -N 1 -t 480 --ntasks-per-node=6 -p gpu --pty bash
# module load swenv/default-env/devel
# module load lang/R/3.6.0-foss-2019a-bare

setwd("/Users/Muhammad_Ali/Downloads/WUSTL_Tasks/Task9_DESeq2")
library(DESeq2)

counts <- read.table("mayo.path_aging.con.salmon.gene.counts.txt", sep="\t", header=T, row.names=1, check.names = FALSE, stringsAsFactors=F)
pheno <- read.csv("mayo.path_aging.con.phenotype.csv", sep=",", header=T, stringsAsFactors=F)

# sanity check
identical(colnames(counts), pheno$UID) # TRUE

# pheno file correction
pheno$AgeAtDeath <- as.factor(pheno$AgeAtDeath)
pheno$Sex <- as.factor(pheno$Sex)
pheno$Diagnosis <- as.factor(gsub(" ", "_", pheno$Diagnosis))
row.names(pheno) <- pheno$UID

# convert count matrix to integer
for (i in 1:ncol(counts))
 {
  counts[,i] <- as.integer(counts[,i])
 }

# make a DESeq object (https://sites.tufts.edu/biotools/files/2019/04/bioinformatics_for_rnaseq_day2.pdf check slide 25)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = pheno, design = ~ AgeAtDeath + Sex + Diagnosis)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="Diagnosis_Pathologic_Aging_vs_Control") # get the statistics
dim(res) # 55857     6
summary(res)

# out of 42936 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1732, 4%
# LFC < 0 (down)     : 2919, 6.8%
# outliers [1]       : 34, 0.079%
# low counts [2]     : 15567, 36%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

sum(res$padj < 0.05, na.rm=TRUE) # 3224
sum(res$pvalue < 0.05, na.rm=TRUE) # 7029
res <- res[order(res$padj),]
write.table(res, file="DiffExp_DESeq2.txt", sep="\t", quote=F, row.names=F)
save(dds, res, counts, pheno, file="dds.RData")

# directory: /mnt/irisgpfs/users/mali/Tasks/WUSTL_Tasks/Task9_DESeq2

head(res)
log2 fold change (MLE): Diagnosis Pathologic Aging vs Control 
Wald test p-value: Diagnosis Pathologic Aging vs Control 
DataFrame with 6 rows and 6 columns
                 baseMean     log2FoldChange              lfcSE
                <numeric>          <numeric>          <numeric>
TMSB4XP8 136.151158952484  -1.69480630773748  0.153301073031935
YBX1P10  10.8575186411274  -3.03142073102026  0.322047366285225
SPCS2P4  53.9110959484258  -2.11881688134672  0.226189921658405
TMSB4XP4 62.1209219895609  -1.75338883367661  0.191871469470421
DYNLL2      15250.4222931 -0.713643105699521 0.0880789615123444
RPL35P5  8.81066416977477  -1.92073920764157  0.246191771619061
                      stat               pvalue                 padj
                 <numeric>            <numeric>            <numeric>
TMSB4XP8 -11.0554105996664 2.06390263789487e-28 5.64188425094942e-24
YBX1P10  -9.41296544662766 4.82332796734826e-21 6.59252466577161e-17
SPCS2P4  -9.36742391443324 7.43243598680164e-21 6.77243567117365e-17
TMSB4XP4  -9.1383509935901 6.34118313115712e-20 4.33356455183277e-16
DYNLL2   -8.10231062498964  5.3925044217513e-16 2.94819001745987e-12
RPL35P5  -7.80180099038233 6.10298388783729e-15 2.78051945929867e-11




##################


setwd("/Users/Muhammad_Ali/Downloads/WUSTL_Tasks/Task9_DESeq2")
library(DESeq2)

counts <- read.table("mayo.path_aging.con.salmon.gene.counts.txt", sep="\t", header=T, row.names=1, check.names = FALSE, stringsAsFactors=F)
pheno <- read.csv("mayo.path_aging.con.phenotype.csv", sep=",", header=T, stringsAsFactors=F)

# sanity check
identical(colnames(counts), pheno$UID) # TRUE

# pheno file correction
pheno$AgeAtDeath <- as.factor(pheno$AgeAtDeath)
pheno$Sex <- as.factor(pheno$Sex)
pheno$Diagnosis <- as.factor(gsub(" ", "_", pheno$Diagnosis))
row.names(pheno) <- pheno$UID

# convert count matrix to integer
for (i in 1:ncol(counts))
 {
  counts[,i] <- as.integer(counts[,i])
 }

# make a DESeq object (https://sites.tufts.edu/biotools/files/2019/04/bioinformatics_for_rnaseq_day2.pdf check slide 25)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = pheno, design = ~ AgeAtDeath + Sex + Diagnosis)
dds <- DESeq(dds, test="LRT", reduced = ~ AgeAtDeath + Sex) # https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/DESeq # https://support.bioconductor.org/p/106649/
resultsNames(dds) # lists the coefficients

res <- results(dds, name="Diagnosis_Pathologic_Aging_vs_Control") # get the statistics
dim(res) # 55857     6
summary(res)

# out of 42936 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1667, 3.9%
# LFC < 0 (down)     : 2872, 6.7%
# outliers [1]       : 34, 0.079%
# low counts [2]     : 15567, 36%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

sum(res$padj < 0.05, na.rm=TRUE) # 3175
sum(res$pvalue < 0.05, na.rm=TRUE) # 6959
res <- res[order(res$padj),]
write.table(res, file="DiffExp_DESeq2_LRT.txt", sep="\t", quote=F, row.names=F)
save(dds, res, counts, pheno, file="dds_LRT.RData")

# save only the significant DEGs (and their counts) for R Shiny App table and boxplots
de.result <- as.data.frame(res)
de.result <- na.omit(de.result)
de.result <- de.result[de.result$padj < 0.05,]
de.result[,1:4] <- round(de.result[,1:4], 2)
de.result[,5:6] <- format(de.result[,5:6], digits=3)
expression.values <- counts
expression.values <- expression.values[row.names(de.result),]
expression.values <- as.matrix(expression.values)
condition <- pheno$Diagnosis
save(expression.values, condition, de.result, file="ShinyApp_Input.RData")


head(res)
log2 fold change (MLE): Diagnosis Pathologic Aging vs Control 
LRT p-value: '~ AgeAtDeath + Sex + Diagnosis' vs '~ AgeAtDeath + Sex' 
DataFrame with 6 rows and 6 columns
                 baseMean     log2FoldChange              lfcSE
                <numeric>          <numeric>          <numeric>
TMSB4XP8 136.151158952484  -1.69480630774407  0.153301073180988
YBX1P10  10.8575186411274  -3.03142073102252  0.322047366261412
SPCS2P4  53.9110959484258  -2.11881688135142  0.226189921554572
TMSB4XP4 62.1209219895609  -1.75338883368047  0.191871469518722
DYNLL2      15250.4222931 -0.713643105700687 0.0880789617864137
RPL35P5  8.81066416977477  -1.92073920775202   0.24619177231926
                     stat               pvalue                 padj
                <numeric>            <numeric>            <numeric>
TMSB4XP8 130.601430242898 3.02664778167179e-30   8.273644375978e-26
YBX1P10  96.4168915072769 9.30697204071798e-23 1.27207693852533e-18
SPCS2P4   90.208191266692 2.14372722968982e-21 1.95336425169337e-17
TMSB4XP4 87.6445369312147 7.83406287049757e-21 5.35379856569804e-17
DYNLL2   66.2669014243006 3.93816395002849e-16 2.15307299475957e-12
RPL35P5  63.9997348949099 1.24435954897199e-15  5.6693021051164e-12



###############################

R shiny App: https://bioinformatics-core-shared-training.github.io/shiny-bioinformatics/use-case


# load("dds_LRT.RData")
# expression.values <- counts
# er.status <- pheno$Diagnosis
# expression.values <- as.matrix(expression.values)
# boxplot(expression.values[2,] ~ er.status)


library(shiny)

ui <- fluidPage(
   
   # Application title
   titlePanel("Interrogating the Aging vs. Control dataset"),
   
   sidebarLayout(
      sidebarPanel(
         selectInput("thegene","Gene to Analyse",
                     choices=c("TMSB4XP8","YBX1P10","SPCS2P4"),
                       selected  = "YBX1P10")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("boxplot")
      )
   )
)

server <- function(input, output) {

 load("ShinyApp_Input.RData")

 output$boxplot <- renderPlot({
  gene <- input$thegene
  counts <- expression.values[gene,]
  boxplot(counts ~ condition)
  })
}


# Run the application 
shinyApp(ui = ui, server = server)



###############################

R shiny App: https://bioinformatics-core-shared-training.github.io/shiny-bioinformatics/use-case


# load("dds_LRT.RData")
# expression.values <- counts
# er.status <- pheno$Diagnosis
# expression.values <- as.matrix(expression.values)
# boxplot(expression.values[2,] ~ er.status)


library(shiny)

load("ShinyApp_Input.RData")
DEGs <- row.names(de.result)

ui <- fluidPage(
   
   # Application title
   titlePanel("Interrogating the Aging vs. Control dataset"),
   
   sidebarLayout(
      sidebarPanel(
         selectInput("thegene","Gene to Analyse",
                     choices=DEGs,
                       selected  = "YBX1P10")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        plotOutput("boxplot"),
        h2("Summary of the DE"),
        verbatimTextOutput("DE_Results")
        
      )
   )

# fluidRow(dataTableOutput('DE_Results'))

)

server <- function(input, output) {

# load("ShinyApp_Input.RData")

  output$boxplot <- renderPlot({
  gene <- input$thegene
  counts <- expression.values[gene,]
  boxplot(counts ~ condition)
  })
  
  output$DE_Results <- renderPrint({
    de.result
  })
}


# Run the application 
shinyApp(ui = ui, server = server)

