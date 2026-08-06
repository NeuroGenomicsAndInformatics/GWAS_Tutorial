# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign command-line arguments to variables
RefBim_file <- args[1]
MydataBim_file <- args[2]
OutputDir <- args[3]
# Read input files
RefBim <- read.table(RefBim_file, header = FALSE)
MydataBim <- read.table(MydataBim_file, header = FALSE)

RefBim_unique<- RefBim[!duplicated(paste(RefBim$V1, RefBim$V4, RefBim$V5, RefBim$V6)), ]
RefBim_unique$SNP_ID <- paste0("chr", RefBim_unique$V1, ":", RefBim_unique$V4)

head(RefBim_unique)

MydataBim_unique<- MydataBim[!duplicated(paste(MydataBim$V1, MydataBim$V4, MydataBim$V5, MydataBim$V6)), ]
removed_duplicates <- MydataBim[duplicated(paste(MydataBim$V1, MydataBim$V4, MydataBim$V5, MydataBim$V6)), ]
cat("Duplicates in MydataBim per chr, pos, ref, alt:\n")
dim(removed_duplicates)
MydataBim_unique$SNP_ID <- paste0("chr", MydataBim_unique$V1, ":", MydataBim_unique$V4)

#If there are no duplicate SNPs, merge bim files.
commonSNPsMerge = merge.data.frame(RefBim_unique, MydataBim_unique, by = "SNP_ID")
#V2 is SNPID chr:bp
cat("Common snps between mydata and ref data:\n")
dim(commonSNPsMerge)

#Function Definition
allele.qc = function(a1,a2,ref1,ref2) {
  # a1 and a2 are the first data-set
  # ref1 and ref2 are the 2nd data-set
  # Make all the alleles into upper-case, as A,T,C,G:
  a1 = toupper(a1)
  a2 = toupper(a2)
  ref1 = toupper(ref1)
  ref2 = toupper(ref2)
  # Strand flip, to change the allele representation in the 2nd data-set
  strand_flip = function(ref) {
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip
  }
  flip1 = strand_flip(ref1)
  flip2 = strand_flip(ref2)
  snp = list()
  # Remove strand ambiguous SNPs (scenario 3)
  snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
  # Remove non-ATCG coding
  snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
  snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
  # as long as scenario 1 is involved, sign_flip will return TRUE
  snp[["sign_flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
  # as long as scenario 2 is involved, strand_flip will return TRUE
  snp[["strand_flip"]] = (a1 == flip1 & a2 == flip2) | (a1 == flip2 & a2 == flip1)
  # remove other cases, eg, tri-allelic, one dataset is A C, the other is A G, for example.
  exact_match = (a1 == ref1 & a2 == ref2) 
  snp[["keep"]][!(exact_match | snp[["sign_flip"]] | snp[["strand_flip"]])] = F
  return(snp)
}
head(commonSNPsMerge)
#Call function allele.qc
qc = allele.qc(commonSNPsMerge[,'V5.x'] , commonSNPsMerge[,'V6.x'] , commonSNPsMerge[,'V5.y'], commonSNPsMerge[,'V6.y'])

keep = unlist(qc[1])
#Shows the number of SNPs kept and dropped
cat("Number of SNPs kept and dropped:\n")
table(keep)

signFlip = unlist(qc[2])
#Shows the number of SNPs with alleles flipped
cat("Number of SNPs with alleles flipped:\n")
table(signFlip)

strandFlip = unlist(qc[3])
#Shows the number of SNPs that are strand flipped
cat("Number of SNPs that are strand flipped:\n")
table(strandFlip)

qcMerge = cbind.data.frame(keep, signFlip, strandFlip)

qcCommonSNPsMerge = cbind.data.frame(qcMerge, commonSNPsMerge)

#Create data frame of SNPs kept, dropped, and strand flipped
qcCommonSNPsMergeKeep = subset(qcCommonSNPsMerge, qcCommonSNPsMerge$keep == "TRUE")
qcCommonSNPsMergeDrop = subset(qcCommonSNPsMerge, qcCommonSNPsMerge$keep == "FALSE")
qcCommonSNPsMergeStrandFlip = subset(qcCommonSNPsMergeKeep, qcCommonSNPsMergeKeep$strandFlip == "TRUE")
qcCommonSNPsMergeSignFlip = subset(qcCommonSNPsMergeKeep, qcCommonSNPsMergeKeep$signFlip == "TRUE")

#Creates a list of SNPs to extract from each data set; 
cat("Writing SNPs to use in the pseudo-GWAS:\n")

refdata<-qcCommonSNPsMergeKeep[ , c(1:10)]
refdata<-refdata[!duplicated(paste(refdata$SNP_ID)), ]

extractref = refdata[ , 6]
write.table(extractref, paste0(OutputDir, "/extractrefSnps.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)

extractmydata = qcCommonSNPsMergeKeep[ , 12]
write.table(extractmydata, paste0(OutputDir, "/extractmydataSnps.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)

#Make list of SNPs to be strand flipped
flipmydata = qcCommonSNPsMergeStrandFlip[ , 12]
write.table(flipmydata, paste0(OutputDir, "/mydataStrandFlip.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)

#Make list of snps to be used for snp update with sign(allele) flip
SignFlipmydata = qcCommonSNPsMergeSignFlip[ ,c(12,15,16,16,15)]
write.table(SignFlipmydata, paste0(OutputDir, "/mydataalleleFlip.txt"),sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


cat("Done!!\n")