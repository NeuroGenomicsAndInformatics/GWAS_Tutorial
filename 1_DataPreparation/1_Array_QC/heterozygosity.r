#!/bin/Rscript
# This script generates heterozygosity plots and a list of individuals outside 3 SDs of the mean.
args <- commandArgs(trailingOnly = TRUE)
OUTDIR <- args[1]
HET_FILE <- args[2]

het <- read.table(HET_FILE, head=TRUE)
het$meanHet = (het$N.NM. - het$O.HOM.)/het$N.NM.
het$HET_RATE = (het$N.NM. - het$O.HOM.)/het$N.NM.

# Plot heterozygosity distribution
png(paste(OUTDIR, "heterozygosity_distribution.png", sep="/"))
hist(het$meanHet, xlab="Heterozygosity Rate", ylab="Frequency", main= "Heterozygosity Rate")
abline(v=mean(het$meanHet)-(3*sd(het$meanHet)),col="RED",lty=2)
abline(v=mean(het$meanHet)+(3*sd(het$meanHet)),col="RED",lty=2)
dev.off()

# The following code generates a list of individuals who deviate more than 3 standard deviations 
#from the heterozygosity rate mean.

het_fail = subset(het, (het$HET_RATE < mean(het$HET_RATE)-3*sd(het$HET_RATE)) | (het$HET_RATE > mean(het$HET_RATE)+3*sd(het$HET_RATE)));
het_fail$HET_DST = (het_fail$HET_RATE-mean(het$HET_RATE))/sd(het$HET_RATE);

write.table(het_fail, paste(OUTDIR, "fail-het-qc.txt", sep="/"), row.names=FALSE)

het_fail_F = subset(het, (het$F > 0.25) | (het$F < -0.25))
write.table(het_fail_F, paste(OUTDIR, "fail-het-qc-F0.25.txt", sep="/"), row.names=FALSE)