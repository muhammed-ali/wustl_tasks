
unadjusted <- read.table("LR.assoc.linear", sep="", header=T, stringsAsFactors=F)
# dim(unadjusted) # 1223350      12

unadjusted <- unadjusted[order(unadjusted$P),]
unadjusted <- unadjusted[unadjusted$TEST == "ADD",]
top_SNP <- unadjusted[unadjusted$P == min(unadjusted$P),]
top_SNP <- top_SNP[,c("SNP", "BETA", "L95", "P")]
writeLines(paste("Top SNP of the analysis: "))
top_SNP

