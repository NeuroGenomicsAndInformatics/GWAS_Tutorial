#!/bin/bash
# This script is for doing the basic QC that each batch needs before imputation.
IN_BFILE=$1
OUTDIR=$2
NAMEBASE=$3
REF_FASTA=$4
INDEX_NUM=$(echo ${IN_BFILE##*/} | cut -d- -f1)
printf -v TEMP "%s%02d%s" ${OUTDIR}/ $((10#${INDEX_NUM} + 1)) -${NAMEBASE}
printf -v TEMP2 "%s%02d%s" ${OUTDIR}/ $((10#${INDEX_NUM} + 2)) -${NAMEBASE}
printf -v TEMP3 "%s%02d%s" ${OUTDIR}/ $((10#${INDEX_NUM} + 3)) -${NAMEBASE}
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# fix flips and swaps
plink2 --bfile ${IN_BFILE} \
	--recode vcf id-paste=iid bgz \
	--out ${TEMP}

bcftools annotate \
	--rename-chrs ${SCRIPT_DIR}/chr_change.txt \
	-Oz -o ${TEMP}.chrs.vcf.gz \
	${TEMP}.vcf.gz

bcftools +fixref ${TEMP}.chrs.vcf.gz \
	-Oz -o ${TEMP}.chrs.fixref.vcf.gz \
	-- -f ${REF_FASTA} \
	-m flip

plink2 --vcf ${TEMP}.chrs.fixref.vcf.gz \
	--double-id \
	--snps-only just-acgt \
	--vcf-half-call missing \
	--output-chr chrM \
	--make-bed \
	--out ${TEMP}.fixref

# Missingness filters
# SNPs genotyped in 95% of samples
plink2 --bfile ${TEMP}.fixref \
	--geno 0.05 \
	--make-bed \
	--out ${TEMP2}.geno05

# Individuals that don't have 95% of SNPs
plink2 --bfile ${TEMP2}.geno05 \
	--mind 0.05 \
	--make-bed \
	--out ${TEMP2}.geno05.mind05

# SNPs genotyped in 98% of samples
plink2 --bfile ${TEMP2}.geno05.mind05 \
	--geno 0.02 \
	--make-bed \
	--out ${TEMP2}.geno02.mind05

# Individuals that don't have 98% of SNPs
plink2 --bfile ${TEMP2}.geno02.mind05 \
	--mind 0.02 \
	--make-bed \
	--out ${TEMP2}.geno02.mind02

# check and remove palindromic SNPs
awk '($5=="A" && $6=="T") || ($5=="T" && $6=="A") || ($5=="G" && $6=="C") || ($5=="C" && $6=="G") {print $2}' ${TEMP2}.geno02.mind02.bim > ${OUTDIR}/palindromic_SNPs.txt

plink2 --bfile ${TEMP2}.geno02.mind02 \
	--exclude ${OUTDIR}/palindromic_SNPs.txt \
	--make-bed \
	--out ${TEMP2}.geno02.mind02.nopal

# Update IDs to CHR:POS:REF:ALT to match with GWAS master file and maf filter
plink2 --bfile ${TEMP2}.geno02.mind02.nopal \
	--set-all-var-ids @:#:\$r:\$a \
	--maf 0.01 \
	--output-chr chrM \
	--make-bed \
	--out ${TEMP3}

# Clean up
rm ${TEMP}.vcf.gz
rm ${TEMP}.chrs.vcf.gz
rm ${TEMP}.chrs.fixref.vcf.gz
rm ${TEMP}.fixref.{bed,bim,fam}
rm ${TEMP2}.geno05.{bed,bim,fam}
rm ${TEMP2}.geno05.mind05.{bed,bim,fam}
rm ${TEMP2}.geno02.mind05.{bed,bim,fam}
rm ${TEMP2}.geno02.mind02.{bed,bim,fam}
rm ${TEMP2}.geno02.mind02.nopal.{bed,bim,fam}
