#!/bin/bash
# This script is for doing the basic QC that each batch needs before imputation.
IN_BFILE=$1
OUTDIR=$2
NAMEBASE=$3
SEX_FILE=$4
EUR_LIST=$5
PHENO_FILE=$6
INDEX_NUM=$(echo ${IN_BFILE##*/} | cut -d- -f1)
printf -v TEMP "%s%02d%s" ${OUTDIR}/ $((10#${INDEX_NUM} + 1)) -${NAMEBASE}
printf -v TEMP2 "%s%02d%s" ${OUTDIR}/ $((10#${INDEX_NUM} + 2)) -${NAMEBASE}
printf -v TEMP3 "%s%02d%s" ${OUTDIR}/ $((10#${INDEX_NUM} + 3)) -${NAMEBASE}
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Update sex and perform a sex check
plink2 --bfile ${IN_BFILE} \
	--update-sex ${SEX_FILE} \
	--output-chr chrM \
	--make-bed \
	--out ${TEMP}.wsex

plink1.9 --bfile ${TEMP}.wsex \
	--split-x b38 no-fail \
	--keep-allele-order \
	--allow-no-sex \
	--make-bed \
	--out ${TEMP}.wsex.split-x

plink1.9 --bfile ${TEMP}.wsex.split-x \
	--check-sex \
	--keep-allele-order \
	--allow-no-sex \
	--make-bed \
	--out ${TEMP}.wsex.split-x.check

grep PROBLEM ${TEMP}.wsex.split-x.check.sexcheck > ${OUTDIR}/${NAMEBASE}_PROBLEM_SEX.txt

# Heterozygosity
plink2 --bfile ${TEMP}.wsex \
	--chr 1-22 \
	--make-bed \
	--out ${TEMP2}_auto

plink2 --bfile ${TEMP2}_auto \
	--set-all-var-ids @:#:\$r:\$a \
	--output-chr chrM \
	--rm-dup force-first \
	--make-bed \
	--out ${TEMP2}.ids

plink2 --bfile ${TEMP2}.ids \
	--indep-pairwise 50 5 0.2 \
	--out ${OUTDIR}/indepSNP

plink1.9 --bfile ${TEMP2}.ids \
	--extract ${OUTDIR}/indepSNP.prune.in \
	--het \
	--keep-allele-order \
	--allow-no-sex \
	--make-bed \
	--out ${TEMP2}.ids_pruned

Rscript ${SCRIPT_DIR}/heterozygosity.r \
	${OUTDIR} \
	${TEMP2}.ids_pruned.het

# Differential missingness, HWE, and setting male het calls to missing
Rscript ${SCRIPT_DIR}/new_basic_QC.R \
	${TEMP}.wsex \
	$OUTDIR/ \
	$NAMEBASE \
	$EUR_LIST \
	$PHENO_FILE
	
plink2 --bfile ${OUTDIR}/${NAMEBASE}.filts.hwe.hh \
	--update-sex ${SEX_FILE} \
	--set-all-var-ids @:#:\$r:\$a \
	--output-chr chrM \
	--make-bed \
	--out ${TEMP3}

# Clean up
# Note: have to provide filenames for both with and without diagnosis data
rm ${TEMP}.wsex.{bed,bim,fam}
rm ${TEMP}.wsex.split-x.{bed,bim,fam}
[[ -f ${TEMP}.wsex.split-x.nosex ]] && rm ${TEMP}.wsex.split-x.nosex
rm ${TEMP}.wsex.split-x.check.{bed,bim,fam}
[[ -f ${TEMP}.wsex.split-x.check.nosex ]] && rm ${TEMP}.wsex.split-x.check.nosex
rm ${TEMP2}_auto.{bed,bim,fam}
rm ${TEMP2}.ids.{bed,bim,fam}
rm ${TEMP2}.ids_pruned.{bed,bim,fam}
[[ -f ${TEMP2}.ids_pruned.nosex ]] && rm ${TEMP2}.ids_pruned.nosex
rm ${OUTDIR}/${NAMEBASE}.EUR.0.A.{bed,bim,fam}
[[ -f ${OUTDIR}/${NAMEBASE}.EUR.0.A.dx.nosex ]] && rm ${OUTDIR}/${NAMEBASE}.EUR.0.A.dx.nosex
[[ -f ${OUTDIR}/${NAMEBASE}.EUR.0.A.sex.nosex ]] && rm ${OUTDIR}/${NAMEBASE}.EUR.0.A.sex.nosex
rm ${OUTDIR}/${NAMEBASE}.EUR.0.A.for_sex_missingness.{bed,bim,fam}
[[ -f ${OUTDIR}/${NAMEBASE}.EUR.0.A.for_sex_missingness.nosex ]] && rm ${OUTDIR}/${NAMEBASE}.EUR.0.A.for_sex_missingness.nosex
if [[ -f ${OUTDIR}/${NAMEBASE}.dx_filt.bed ]]; then
	rm ${OUTDIR}/${NAMEBASE}.dx_filt.{bed,bim,fam}
	rm ${OUTDIR}/${NAMEBASE}.dx_filt.sex_filt.{bed,bim,fam}
	rm ${OUTDIR}/${NAMEBASE}.dx_filt.sex_filt.hwe.{bed,bim,fam}
	rm ${OUTDIR}/${NAMEBASE}.dx_filt.sex_filt.hwe.hwex.{bed,bim,fam}
	rm ${OUTDIR}/${NAMEBASE}.dx_filt.sex_filt.hwe.hwex.par.{bed,bim,fam}
else
	rm ${OUTDIR}/${NAMEBASE}.sex_filt.{bed,bim,fam}
	rm ${OUTDIR}/${NAMEBASE}.sex_filt.hwe.{bed,bim,fam}
	rm ${OUTDIR}/${NAMEBASE}.sex_filt.hwe.hwex.{bed,bim,fam}
	rm ${OUTDIR}/${NAMEBASE}.sex_filt.hwe.hwex.par.{bed,bim,fam}
fi
rm ${OUTDIR}/${NAMEBASE}.filts.hwe.hh.{bed,bim,fam}
