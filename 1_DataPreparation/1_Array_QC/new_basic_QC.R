################################################################################
## Libraries

library(openxlsx)
library(plyr)
library(dplyr)
library(data.table) 

args <- commandArgs(trailingOnly = TRUE )

################################################################################
## Main input

pl <- "plink1.9"
pl2 <- "plink2"

in_bfile <- args[1]
out_fold <- args[2]
namebase <- args[3]
subj_list_EUR <- args[4]
if (length(args)==5) {
	ph_file <- args[5]
	dx_filt <- ".dx_filt"
	} else {
	dx_filt <- ""
}

out_file <- namebase

################################################################################
## subset data to EUR and apply basic QC

## extract subjects
command <- paste(pl2, 
                 "--bfile", in_bfile,
                 "--keep", subj_list_EUR,
                 "--make-bed --out", paste(out_fold, out_file, ".EUR.0.A", sep=""),
                 sep= " ")
system(command)

# Differential missingness cases/controls
if (dx_filt == ".dx_filt") {
## get pheno file
ph  <- fread(ph_file, header = TRUE)

command <- paste(pl,
                 "--bfile", paste(out_fold, out_file, ".EUR.0.A", sep=""),
				 "--pheno", ph_file, "--pheno-name", "STATUS",
                 "--test-missing --keep-allele-order",
                 "--out", paste(out_fold, out_file, ".EUR.0.A.dx", sep=""),
                 sep=" ")
system(command)
mis <- fread(paste(out_fold, out_file, ".EUR.0.A.dx.missing", sep=""), header = T)
mis <- mis[which(mis$P < 1e-5),"SNP"]
fwrite(mis, paste(out_fold, out_file, ".EUR.0.A.dx_specific_diff_P1e-5_vars", sep=""),col.names = F)
}
# Differential missingness males/females
command <- paste(pl,
                 "--bfile", paste(out_fold, out_file, ".EUR.0.A", sep=""),
                 "--keep-allele-order",
                 "--make-bed --out", paste(out_fold, out_file, ".EUR.0.A.for_sex_missingness", sep=""),
                 sep=" ")
system(command)
fam <- fread(paste(out_fold, out_file, ".EUR.0.A.for_sex_missingness.fam", sep=""))
#fam[which(fam$V2 %in% unlist(ph[which(ph$SEX==1),"IID"])),"V6"] <- 1
#fam[which(fam$V2 %in% unlist(ph[which(ph$SEX==2),"IID"])),"V6"] <- 2
fam$V6 <- fam$V5
fwrite(fam,paste(out_fold, out_file, ".EUR.0.A.for_sex_missingness.fam", sep=""),sep="\t",col.names=F)
command <- paste(pl,
                 "--bfile", paste(out_fold, out_file, ".EUR.0.A.for_sex_missingness", sep=""),
                 "--test-missing --keep-allele-order",
                 "--out", paste(out_fold, out_file, ".EUR.0.A.sex", sep=""),
                 sep=" ")
system(command)
mis <- fread(paste(out_fold, out_file, ".EUR.0.A.sex.missing", sep=""), header = T)
mis <- mis[which(mis$P < 1e-5),"SNP"]
fwrite(mis, paste(out_fold, out_file, ".EUR.0.A.sex_specific_diff_P1e-5_vars", sep=""),col.names = F)

## Run HWE in females and males and extract deviation
command <- paste(pl2, 
                 "--bfile", paste(out_fold, out_file, ".EUR.0.A", sep=""),
                 "--hardy",
                 "--out", paste(out_fold, out_file, ".EUR.0.B", sep=""),
                 sep= " ")
system(command)
if (file.exists(paste(out_fold, out_file, ".EUR.0.B.hardy", sep=""))) {
ttt <- fread(paste(out_fold, out_file, ".EUR.0.B.hardy", sep=""))
fwrite(ttt[which(ttt$P<1e-6),"ID"],
       paste(out_fold, out_file, ".EUR.0.B.hardy.var_list", sep=""), col.names = F)
}

if (file.exists(paste(out_fold, out_file, ".EUR.0.B.hardy.x", sep=""))) {
tt <- fread(paste(out_fold, out_file, ".EUR.0.B.hardy.x", sep=""))
fwrite(tt[which(tt$P<1e-6),"ID"],
       paste(out_fold, out_file, ".EUR.0.B.hardy.x.var_list", sep=""), col.names = F)
}

# apply the filters
if (dx_filt == ".dx_filt") {
command <- paste(pl2, 
                 "--bfile", in_bfile,
                 "--exclude", paste(out_fold, out_file, ".EUR.0.A.dx_specific_diff_P1e-5_vars", sep=""),
                 "--make-bed --out", paste(out_fold, out_file, ".dx_filt", sep=""),
                 sep= " ")
system(command)

command <- paste(pl2, 
                 "--bfile", paste(out_fold, out_file, ".dx_filt", sep=""),
                 "--exclude", paste(out_fold, out_file, ".EUR.0.A.sex_specific_diff_P1e-5_vars", sep=""),
                 "--make-bed --out", paste(out_fold, out_file, dx_filt, ".sex_filt", sep=""),
                 sep= " ")
system(command)
} else {
command <- paste(pl2, 
                 "--bfile", in_bfile,
                 "--exclude", paste(out_fold, out_file, ".EUR.0.A.sex_specific_diff_P1e-5_vars", sep=""),
                 "--make-bed --out", paste(out_fold, out_file, dx_filt, ".sex_filt", sep=""),
                 sep= " ")
system(command)
}
hwe <- ""
if (file.exists(paste(out_fold, out_file, ".EUR.0.B.hardy.var_list", sep=""))) {
hwe <- ".hwe"
command <- paste(pl2, 
                 "--bfile", paste(out_fold, out_file, dx_filt, ".sex_filt", sep=""),
				 "--exclude", paste(out_fold, out_file, ".EUR.0.B.hardy.var_list", sep=""),
                 "--make-bed --out", paste(out_fold, out_file, dx_filt, ".sex_filt.hwe", sep=""),
                 sep= " ")
system(command)

}

hwex <- ""

if (file.exists(paste(out_fold, out_file, ".EUR.0.B.hardy.x.var_list", sep=""))) {
hwex <- ".hwex"

command <- paste(pl2, 
                 "--bfile", paste(out_fold, out_file, dx_filt, ".sex_filt",hwe, sep=""),
                 "--exclude", paste(out_fold, out_file, ".EUR.0.B.hardy.x.var_list", sep=""),
                 "--make-bed --out", paste(out_fold, out_file, dx_filt, ".sex_filt",hwe,hwex, sep=""),
                 sep= " ")
system(command)

}

## set heterozygote hard calls to missing
command <- paste(pl2,
                 "--bfile", paste(out_fold, out_file, dx_filt, ".sex_filt",hwe,hwex, sep=""),
                 "--set-hh-missing",
                 "--make-bed --out", paste(out_fold, out_file, ".filts.hwe.hh", sep=""),
                 sep= " ")
system(command)

