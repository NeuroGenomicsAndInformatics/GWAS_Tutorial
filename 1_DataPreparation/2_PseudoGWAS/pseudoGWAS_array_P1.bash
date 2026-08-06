#!/bin/bash
# This script is for performing a pseudoGWAS with 2 batches of our data. This script prepares the data.

OUTDIR=$1
IN_BFILEA=$2
NAMEBASEA=$3
EURCO_FILEA=$4
IN_BFILEB=$5
NAMEBASEB=$6
EURCO_FILEB=$7

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Get just EUR Controls
plink2 --bfile ${IN_BFILEA} \
	--keep ${EURCO_FILEA} \
	--make-bed \
	--out ${OUTDIR}/${NAMEBASEA}_EURCO

plink2 --bfile ${IN_BFILEB} \
	--keep ${EURCO_FILEB} \
	--make-bed \
	--out ${OUTDIR}/${NAMEBASEB}_EURCO

# Get common SNPs and flips
Rscript ${SCRIPT_DIR}/allele_qc.R \
	${OUTDIR}/${NAMEBASEA}_EURCO.bim \
	${OUTDIR}/${NAMEBASEB}_EURCO.bim \
	${OUTDIR} \
	> ${OUTDIR}/allele_qc_report.txt

plink2 --bfile ${OUTDIR}/${NAMEBASEA}_EURCO \
	--extract ${OUTDIR}/extractrefSnps.txt \
	--make-bed \
	--out ${OUTDIR}/${NAMEBASEA}_with_${NAMEBASEB}_snps

plink2 --bfile ${OUTDIR}/${NAMEBASEA}_with_${NAMEBASEB}_snps \
	--set-all-var-ids @:#:\$r:\$a \
	--output-chr chrM \
	--make-bed \
	--out ${OUTDIR}/${NAMEBASEA}_with_${NAMEBASEB}_snps_ids

plink2 --bfile ${OUTDIR}/${NAMEBASEB}_EURCO \
	--extract ${OUTDIR}/extractmydataSnps.txt \
	--make-bed \
	--out ${OUTDIR}/${NAMEBASEB}_with_${NAMEBASEA}_snps

#plink1.9 --bfile ${OUTDIR}/${NAMEBASEB}_with_${NAMEBASEA}_snps \
#	--flip ${OUTDIR}/mydataStrandFlip.txt \
#	--keep-allele-order \
#	--allow-no-sex \
#	--make-bed \
#	--out ${OUTDIR}/${NAMEBASEB}_with_${NAMEBASEA}_snps_flip

plink2 --bfile ${OUTDIR}/${NAMEBASEB}_with_${NAMEBASEA}_snps \
	--set-all-var-ids @:#:\$r:\$a \
	--output-chr chrM \
	--make-bed \
	--out ${OUTDIR}/${NAMEBASEB}_with_${NAMEBASEA}_snps_ids

# Merge
plink1.9 --bfile ${OUTDIR}/${NAMEBASEA}_with_${NAMEBASEB}_snps_ids \
	--bmerge ${OUTDIR}/${NAMEBASEB}_with_${NAMEBASEA}_snps_ids \
	--chr 1-22 \
	--allow-extra-chr \
	--allow-no-sex \
	--keep-allele-order \
	--make-bed \
	--out ${OUTDIR}/${NAMEBASEA}__${NAMEBASEB}_merged

bash ${SCRIPT_DIR}/get2columns_log.bash ${OUTDIR}/${NAMEBASEA}__${NAMEBASEB}_merged.log

plink1.9 --bfile ${OUTDIR}/${NAMEBASEA}__${NAMEBASEB}_merged \
	--exclude ${OUTDIR}/${NAMEBASEA}__${NAMEBASEB}_merged_Col2.txt \
	--allow-no-sex \
	--keep-allele-order \
	--make-bed \
	--out ${OUTDIR}/${NAMEBASEA}__${NAMEBASEB}_merged_clean

# Generate PCs
plink2 --bfile ${OUTDIR}/${NAMEBASEA}__${NAMEBASEB}_merged_clean --geno 0.02 --maf 0.02 --hwe 0.00005 --pca 10 --out ${OUTDIR}/${NAMEBASEA}__${NAMEBASEB}_merged_clean_PCs

rm ${OUTDIR}/${NAMEBASEA}_EURCO.{bed,bim,fam}
rm ${OUTDIR}/${NAMEBASEB}_EURCO.{bed,bim,fam}
rm ${OUTDIR}/${NAMEBASEA}_with_${NAMEBASEB}_snps.{bed,bim,fam}
rm ${OUTDIR}/${NAMEBASEA}_with_${NAMEBASEB}_snps_ids.{bed,bim,fam}
rm ${OUTDIR}/${NAMEBASEB}_with_${NAMEBASEA}_snps.{bed,bim,fam}
rm ${OUTDIR}/${NAMEBASEB}_with_${NAMEBASEA}_snps_ids.{bed,bim,fam}
rm ${OUTDIR}/${NAMEBASEA}__${NAMEBASEB}_merged.{bed,bim,fam}
