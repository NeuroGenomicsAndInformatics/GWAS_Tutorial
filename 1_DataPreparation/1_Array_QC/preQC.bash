#!/bin/bash
# This script is intended to be a combination of the steps typically done in pre-QC. 
# It finds and removes variants with missing alleles, missing chromosomes, and D/I coded Indels and duplicates

# provide full path to plink files. NO FILE EXTENSION
IN_BFILE=$1
OUTDIR=$2
NAMEBASE=$3
REFALT_FILE=$4

awk '$5 == "0" && $6 == "0" {print $2}' ${IN_BFILE}.bim > ${OUTDIR}/double0s.txt
awk '$5 == "-" && $6 == "-" {print $2}' ${IN_BFILE}.bim >> ${OUTDIR}/double0s.txt
awk '$5 == "." && $6 == "." {print $2}' ${IN_BFILE}.bim >> ${OUTDIR}/double0s.txt

plink1.9 --bfile ${IN_BFILE} \
	--keep-allele-order \
	--exclude ${OUTDIR}/double0s.txt \
	--make-bed \
	--out ${OUTDIR}/00-${NAMEBASE}

# Get SNPs with missing alleles in 5th column
awk '$5 == "0" && $6 != "0" {print $1,$2,$3,$4,$5,$6}' ${OUTDIR}/00-${NAMEBASE}.bim > ${OUTDIR}/SNPs_tofix.txt
awk '$5 == "-" && $6 != "-" {print $1,$2,$3,$4,$5,$6}' ${OUTDIR}/00-${NAMEBASE}.bim >> ${OUTDIR}/SNPs_tofix.txt
awk '$5 == "." && $6 != "." {print $1,$2,$3,$4,$5,$6}' ${OUTDIR}/00-${NAMEBASE}.bim >> ${OUTDIR}/SNPs_tofix.txt

# Used SNPmap file to get the REF/ALT SNPs and made the missing_alleles_refalt.txt file with SNPname, REF, and ALT columns
awk 'FNR==NR{a[$2]=$0; next}{if(b=a[$1]){print b, $0;}}' ${OUTDIR}/SNPs_tofix.txt ${REFALT_FILE}  >  ${OUTDIR}/SNPs_info.txt

# Get alleles in case of flip
awk '{print $8,$9}' ${OUTDIR}/SNPs_info.txt \
	| sed  -e 's/A/Y/g; s/T/A/g; s/Y/T/g; s/G/W/g; s/C/G/g; s/W/C/g' > ${OUTDIR}/SNP_flip.txt

# Add flips to end of file
paste ${OUTDIR}/SNPs_info.txt ${OUTDIR}/SNP_flip.txt \
	| awk '{print $0}' > ${OUTDIR}/SNPs_infoflip.txt

# Make list for each case (correct orientation, flip, swap, flip+swap)
awk '$5 == "0" && $6 == $8 {print $2,$5,$6,$9,$8}' ${OUTDIR}/SNPs_infoflip.txt > ${OUTDIR}/SNPs_0allelesfix.txt
awk '$5 == "-" && $6 == $8 {print $2,$5,$6,$9,$8}' ${OUTDIR}/SNPs_infoflip.txt >> ${OUTDIR}/SNPs_0allelesfix.txt
awk '$5 == "." && $6 == $8 {print $2,$5,$6,$9,$8}' ${OUTDIR}/SNPs_infoflip.txt >> ${OUTDIR}/SNPs_0allelesfix.txt

awk '$5 == "0" && $6 == $9 {print $2,$5,$6,$8,$9}' ${OUTDIR}/SNPs_infoflip.txt > ${OUTDIR}/SNPs_0allelesfix2.txt
awk '$5 == "-" && $6 == $9 {print $2,$5,$6,$8,$9}' ${OUTDIR}/SNPs_infoflip.txt >> ${OUTDIR}/SNPs_0allelesfix2.txt
awk '$5 == "." && $6 == $9 {print $2,$5,$6,$8,$9}' ${OUTDIR}/SNPs_infoflip.txt >> ${OUTDIR}/SNPs_0allelesfix2.txt

awk '$5 == "0" && $6 == $10 {print $2,$5,$6,$11,$10}' ${OUTDIR}/SNPs_infoflip.txt > ${OUTDIR}/SNPs_0allelesfix3.txt
awk '$5 == "-" && $6 == $10 {print $2,$5,$6,$11,$10}' ${OUTDIR}/SNPs_infoflip.txt >> ${OUTDIR}/SNPs_0allelesfix3.txt
awk '$5 == "." && $6 == $10 {print $2,$5,$6,$11,$10}' ${OUTDIR}/SNPs_infoflip.txt >> ${OUTDIR}/SNPs_0allelesfix3.txt

awk '$5 == "0" && $6 == $11 {print $2,$5,$6,$10,$11}' ${OUTDIR}/SNPs_infoflip.txt > ${OUTDIR}/SNPs_0allelesfix4.txt
awk '$5 == "-" && $6 == $11 {print $2,$5,$6,$10,$11}' ${OUTDIR}/SNPs_infoflip.txt >> ${OUTDIR}/SNPs_0allelesfix4.txt
awk '$5 == "." && $6 == $11 {print $2,$5,$6,$10,$11}' ${OUTDIR}/SNPs_infoflip.txt >> ${OUTDIR}/SNPs_0allelesfix4.txt

# combine
cat ${OUTDIR}/SNPs_0allelesfix.txt ${OUTDIR}/SNPs_0allelesfix2.txt ${OUTDIR}/SNPs_0allelesfix3.txt ${OUTDIR}/SNPs_0allelesfix4.txt  > ${OUTDIR}/SNPs_0allelesfix_all.txt

# remove duplicates
awk '!seen[$1]++' ${OUTDIR}/SNPs_0allelesfix_all.txt > ${OUTDIR}/allele_change.txt

# update alleles
plink1.9 --bfile ${OUTDIR}/00-${NAMEBASE} \
	--update-alleles ${OUTDIR}/allele_change.txt \
	--keep-allele-order \
	--make-bed \
	--out ${OUTDIR}/00-${NAMEBASE}.alleles

# Find variants with at least one zero allele missing or . missing allele
awk '{if($5 == "0" || $6 == "0" || $5 == "." || $6 == "." || $5 == "-" || $6 == "-") print $0}' ${OUTDIR}/00-${NAMEBASE}.alleles.bim > ${OUTDIR}/vars_allele_missing.txt
echo "variants with missing allele = $(wc -l ${OUTDIR}/vars_allele_missing.txt)" &> ${OUTDIR}/preQC.log

# Find variants with no chr info (chr0)
awk '{if($1 == "0" ) print $2}' ${OUTDIR}/00-${NAMEBASE}.alleles.bim > ${OUTDIR}/vars_chr_missing.txt
echo "variants with missing chr = $(wc -l ${OUTDIR}/vars_chr_missing.txt)" &> ${OUTDIR}/preQC.log

# Check for Indels coded as D/I
awk '{if($5 =="I" || $5 == "D" || $6 == "I" || $6 == "D") print $2}' ${OUTDIR}/00-${NAMEBASE}.alleles.bim > ${OUTDIR}/vars_indels.txt
echo "indel variants coded as D/I = $(wc -l ${OUTDIR}/vars_indels.txt)" &> ${OUTDIR}/preQC.log

plink2 --bfile ${OUTDIR}/00-${NAMEBASE}.alleles \
	--rm-dup list \
	--make-bed \
	--out ${OUTDIR}/00-${NAMEBASE}.rmdups

# remove variants from above steps
cat ${OUTDIR}/vars_allele_missing.txt ${OUTDIR}/vars_chr_missing.txt ${OUTDIR}/vars_indels.txt > ${OUTDIR}/vars_toremove.txt

plink2 --bfile ${OUTDIR}/00-${NAMEBASE}.rmdups \
	--exclude ${OUTDIR}/vars_toremove.txt \
	--make-bed \
	--out ${OUTDIR}/00-${NAMEBASE}.preQC

# Remove zeroed individuals from input plink files if present
plink2 --bfile ${OUTDIR}/00-${NAMEBASE}.preQC \
	--mind 0.1 \
	--make-bed \
	--out ${OUTDIR}/00-${NAMEBASE}.preQC.0inds

# Hide XY to keep them
plink1.9 --bfile ${OUTDIR}/00-${NAMEBASE}.preQC.0inds \
	--merge-x no-fail \
	--keep-allele-order \
	--allow-no-sex \
	--make-bed \
	--out ${OUTDIR}/01-${NAMEBASE}
	
# Clean up
[[ -f ${OUTDIR}/00-${NAMEBASE}.nosex ]] && rm ${OUTDIR}/00-${NAMEBASE}.nosex
rm ${OUTDIR}/00-${NAMEBASE}.alleles.{bed,bim,fam}
[[ -f ${OUTDIR}/00-${NAMEBASE}.alleles.nosex ]] && rm ${OUTDIR}/00-${NAMEBASE}.alleles.nosex
rm ${OUTDIR}/00-${NAMEBASE}.rmdups.{bed,bim,fam}
rm ${OUTDIR}/00-${NAMEBASE}.preQC.{bed,bim,fam}
rm ${OUTDIR}/00-${NAMEBASE}.preQC.0inds.{bed,bim,fam}
