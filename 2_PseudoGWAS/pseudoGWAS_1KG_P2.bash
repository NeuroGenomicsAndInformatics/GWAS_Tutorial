#!/bin/bash
# This script is for performing a pseudoGWAS with a batch of our data with 1000 Genomes. This script runs the GWAS.

OUTDIR=$1
IN_BFILE=$2
NAMEBASE=$3
PHENO_FILE=$4
INDEX_NUM=$(echo ${IN_BFILE##*/} | cut -d- -f1)
printf -v TEMP "%s%02d%s" ${OUTDIR}/ $((10#${INDEX_NUM} + 1)) -${NAMEBASE}

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
	--out ${TEMP}_1KG

awk '{if ($3 < 1e-4) print $2}' ${TEMP}_1KG.assoc.logistic.adjusted > ${OUTDIR}/${NAMEBASE}.blacklist_SNPs_1KG.txt

# Make Manhattan
Rscript ${SCRIPT_DIR}/make_manhattan_update.R ${TEMP}_1KG.assoc.logistic.adjusted

# Clean up
[[ -f ${TEMP}_1KG.nosex ]] && rm ${TEMP}_1KG.nosex
