library(dplyr)
library(qqman)
library(data.table)
library(readr)

args <- commandArgs(trailingOnly = TRUE)

gwasFile <- args[1]

gwasResults <- read.table(gwasFile)
gwasResults <- gwasResults[complete.cases(gwasResults),]
colnames(gwasResults) <- c('CHR','SNP','A1','P')
gwasResults <- gwasResults[gwasResults$CHR != "X",]
ID_split <- strsplit(gwasResults$SNP,split=":")
gwasResults$BP <- as.numeric(sapply(ID_split, "[[", 2))
gwasResults$CHR <- as.numeric(gwasResults$CHR)

chi_squared_stats <- qchisq(1 - gwasResults$P, df=1)
median_observed <- median(chi_squared_stats)
median_expected <- qchisq(0.5, df=1)

lambda <- median_observed / median_expected

tiff(paste(gwasFile,".tiff", sep=""), units="mm", width=190, height=142, compression="lzw", res=1000)
manhattan(gwasResults,main = paste("Manhattan Plot (lambda =", round(lambda, 3), ")"),cex=0.5, cex.axis=1, genomewideline = -log10(5e-8), suggestiveline = -log10(5e-6), col=c("aquamarine4", "darkgoldenrod4", "blue", "blueviolet", "chartreuse4", "firebrick2", "deepskyblue", "hotpink", "darkgoldenrod4", "darkolivegreen", "darkorange", "darkviolet"))
dev.off()

jpeg(paste(gwasFile,".jpeg", sep="")) 
qq(gwasResults$P, main = paste("Q-Q plot of GWAS p-values : log   lambda =",round(lambda,3)),ylim = c(0, 8),xlim = c(0, 8))
dev.off()