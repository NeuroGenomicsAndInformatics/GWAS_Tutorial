#!/bin/bash
# This script is for performing a pseudoGWAS with 2 batches of our data. This script runs the GWAS.

OUTDIR=$1
IN_BFILE=$2
NAMEBASEA=$3
NAMEBASEB=$4
PHENO_FILE=$5

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Run GWAS
plink1.9 --bfile ${IN_BFILE} \
	--logistic hide-covar \
	--adjust \
	--ci 0.95 \
	--pheno ${PHENO_FILE} \
	--pheno-name STATUS \
	--covar ${PHENO_FILE} \
	--covar-name SEX,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
	--keep-allele-order \
	--allow-no-sex \
	--out ${OUTDIR}/${NAMEBASEA}__${NAMEBASEB}_GWAS

awk '{if ($3 < 1e-4) print $2}' ${OUTDIR}/${NAMEBASEA}__${NAMEBASEB}_GWAS.assoc.logistic.adjusted > ${OUTDIR}/${NAMEBASEA}__${NAMEBASEB}.blacklist_SNPs.txt

# Make Manhattan
Rscript ${SCRIPT_DIR}/make_manhattan_update.R ${OUTDIR}/${NAMEBASEA}__${NAMEBASEB}_GWAS.assoc.logistic.adjusted

# Clean up
rm ${OUTDIR}/*.nosex