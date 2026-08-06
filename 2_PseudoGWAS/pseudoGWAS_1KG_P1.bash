#!/bin/bash
# This script is for performing a pseudoGWAS with a batch of our data with 1000 Genomes. This script prepares the data.

OUTDIR=$1
IN_BFILE=$2
NAMEBASE=$3
EURCO_FILE=$4
KG_EUR=$5
INDEX_NUM=$(echo ${IN_BFILE##*/} | cut -d- -f1)
printf -v TEMP "%s%02d%s" ${OUTDIR}/ $((10#${INDEX_NUM} + 1)) -${NAMEBASE}
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Get just EUR Controls
plink2 --bfile ${IN_BFILE} \
	--keep ${EURCO_FILE} \
	--make-bed \
	--out ${OUTDIR}/${NAMEBASE}_EURCO

# Get common SNPs and flips
Rscript ${SCRIPT_DIR}/allele_qc.R \
	$KG_EUR.bim \
	${OUTDIR}/${NAMEBASE}_EURCO.bim \
	${OUTDIR} \
	> ${OUTDIR}/allele_qc_report.txt

plink2 --bfile $KG_EUR \
	--extract ${OUTDIR}/extractrefSnps.txt \
	--make-bed \
	--out ${OUTDIR}/1KG_with_${NAMEBASE}_snps

plink2 --bfile ${OUTDIR}/1KG_with_${NAMEBASE}_snps \
	--set-all-var-ids @:#:\$r:\$a \
	--output-chr chrM \
	--make-bed \
	--out ${OUTDIR}/1KG_with_${NAMEBASE}_snps_ids

plink2 --bfile ${OUTDIR}/${NAMEBASE}_EURCO \
	--extract ${OUTDIR}/extractmydataSnps.txt \
	--make-bed \
	--out ${OUTDIR}/${NAMEBASE}_with_1KG_snps

plink2 --bfile ${OUTDIR}/${NAMEBASE}_with_1KG_snps \
	--rm-dup force-first \
	--make-bed \
	--out ${OUTDIR}/${NAMEBASE}_with_1KG_snps_dups

plink1.9 --bfile ${OUTDIR}/${NAMEBASE}_with_1KG_snps_dups \
	--flip ${OUTDIR}/mydataStrandFlip.txt \
	--keep-allele-order \
	--allow-no-sex \
	--make-bed \
	--out ${OUTDIR}/${NAMEBASE}_with_1KG_snps_flip

plink2 --bfile ${OUTDIR}/${NAMEBASE}_with_1KG_snps_flip \
	--set-all-var-ids @:#:\$r:\$a \
	--output-chr chrM \
	--make-bed \
	--out ${OUTDIR}/${NAMEBASE}_with_1KG_snps_flip_ids

# Merge
plink1.9 --bfile ${OUTDIR}/1KG_with_${NAMEBASE}_snps_ids \
	--bmerge ${OUTDIR}/${NAMEBASE}_with_1KG_snps_flip_ids \
	--chr 1-22 \
	--maf 0.01 \
	--allow-no-sex \
	--keep-allele-order \
	--make-bed \
	--out ${TEMP}_1KG_merged

bash ${SCRIPT_DIR}/get2columns_log.bash ${TEMP}_1KG_merged.log

if [ -s ${TEMP}_1KG_merged_Col2.txt ]; then

plink1.9 --bfile ${TEMP}_1KG_merged \
	--exclude ${TEMP}_1KG_merged_Col2.txt \
	--allow-no-sex \
	--keep-allele-order \
	--make-bed \
	--out ${TEMP}_1KG_merged_clean

else
	cp ${TEMP}_1KG_merged.bed ${TEMP}_1KG_merged_clean.bed
	cp ${TEMP}_1KG_merged.bim ${TEMP}_1KG_merged_clean.bim
	cp ${TEMP}_1KG_merged.fam ${TEMP}_1KG_merged_clean.fam
fi


# Generate PCs
#plink2 --bfile ${TEMP}_1KG_merged_clean --geno 0.02 --maf 0.02 --hwe 0.00005 --pca 20 --out ${TEMP}_1KG_merged_clean_PCs

# Clean up
rm ${OUTDIR}/${NAMEBASE}_EURCO.{bed,bim,fam}
rm ${OUTDIR}/1KG_with_${NAMEBASE}_snps.{bed,bim,fam}
rm ${OUTDIR}/1KG_with_${NAMEBASE}_snps_ids.{bed,bim,fam}
rm ${OUTDIR}/${NAMEBASE}_with_1KG_snps.{bed,bim,fam}
rm ${OUTDIR}/${NAMEBASE}_with_1KG_snps_dups.{bed,bim,fam}
rm ${OUTDIR}/${NAMEBASE}_with_1KG_snps_flip.{bed,bim,fam}
rm ${OUTDIR}/${NAMEBASE}_with_1KG_snps_flip_ids.{bed,bim,fam}
rm ${TEMP}_1KG_merged.{bed,bim,fam}
[[ -f ${TEMP}_1KG_merged.nosex ]] && rm ${TEMP}_1KG_merged.nosex
