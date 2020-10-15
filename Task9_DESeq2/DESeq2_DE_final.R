# NOTE: First rename some of the genes which have been converted to 1-Mar, 2-Mar.. ect. due to opening/saving the provided count data file in excel.

setwd("/Users/Muhammad_Ali/Downloads/WUSTL_Tasks/Task9_DESeq2")

suppressMessages(library(DESeq2))
suppressMessages(library(shiny))

counts <- read.table("mayo.path_aging.con.salmon.gene.counts.txt", sep="\t", header=T, row.names=1, check.names = FALSE, stringsAsFactors=F)
pheno <- read.csv("mayo.path_aging.con.phenotype.csv", sep=",", header=T, stringsAsFactors=F)
padj_threshold <- 0.05

writeLines(paste("Input files read for conducting DEA"))

# sanity check
identical(colnames(counts), pheno$UID) # TRUE

# pheno file correction
pheno$AgeAtDeath <- as.factor(pheno$AgeAtDeath)
pheno$Sex <- as.factor(pheno$Sex)
pheno$Diagnosis <- as.factor(gsub(" ", "_", pheno$Diagnosis))
row.names(pheno) <- pheno$UID

# convert count matrix from numeric to integer
counts <- as.matrix(counts)
mode(counts) <- "integer"

# make a DESeq object (https://sites.tufts.edu/biotools/files/2019/04/bioinformatics_for_rnaseq_day2.pdf check slide 25)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = pheno, design = ~ AgeAtDeath + Sex + Diagnosis)
dds <- DESeq(dds)
# resultsNames(dds) # lists all the coefficients
res <- results(dds, name="Diagnosis_Pathologic_Aging_vs_Control") # get the statistics
writeLines(paste("Dimensions of the DESeq2 DEA results table: "))
dim(res) # 55857     6

writeLines(paste("Summary of DEA: "))
summary(res)

writeLines(paste("No.of significant (padj < 0.05) DEGs: "))
sum(res$padj < padj_threshold, na.rm=TRUE) # 3224
writeLines(paste("No.of nominal (pvalue < 0.05) DEGs: "))
sum(res$pvalue < padj_threshold, na.rm=TRUE) # 7029
res <- res[order(res$padj),]
write.table(res, file="diffExp_DESeq2_complete.txt", sep="\t", quote=F, row.names=F)

# save only the significant DEGs (and their counts) for R shiny app's summary table and boxplots
de_result <- as.data.frame(res)
de_result <- na.omit(de_result)
de_result <- de_result[de_result$padj < padj_threshold,]
de_result[,1:4] <- round(de_result[,1:4], 2)
de_result[,5:6] <- format(de_result[,5:6], digits=3)
expression_values <- counts
expression_values <- expression_values[row.names(de_result),]
expression_values <- as.matrix(expression_values)
condition <- pheno$Diagnosis
save(expression_values, condition, de_result, file="ShinyApp_Input.RData")

writeLines(paste("DEA is complete, now running the shinny app script"))

###############################

# R Shiny app script to show the boxplot of differentially expressed genes (DEG) and table of DE result

suppressMessages(library(shiny))
load("ShinyApp_Input.RData")
DEGs <- row.names(de_result)

ui <- fluidPage(
   
   titlePanel("Interrogating the Aging vs. Control dataset"),
   sidebarLayout(
      sidebarPanel(
         selectInput("thegene","Gene to Analyse",
                     choices=DEGs,
                       selected  = "TMSB4XP8")
      ),
      
      # Show boxplot and DE table in the panels
      mainPanel(
        plotOutput("boxplot"),
        h2("Summary of the DE"),
        verbatimTextOutput("DE_Results")   
      )
   )
)

server <- function(input, output) {

  output$boxplot <- renderPlot({
  gene <- input$thegene
  counts <- expression_values[gene,]
  boxplot(counts ~ condition)
  })
  
  output$DE_Results <- renderPrint({
    de_result
  })
}

writeLines(paste("Open the below URL in your browser to view the results in R shiny app: "))

# Run the application 
shinyApp(ui = ui, server = server)

