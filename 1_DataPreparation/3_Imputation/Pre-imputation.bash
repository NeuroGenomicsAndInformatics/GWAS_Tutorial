#!/bin/bash
# This script is for preparing data to be uploaded to the TOPMed server for imputation.
IN_BFILE=$1
OUTDIR=$2
NAMEBASE=$3
EXCLUDE_LIST=$4

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

FILT_DIR="${SCRIPT_DIR}/gnomAD_QC"

if [[ -n $EXCLUDE_LIST ]]; then
	plink2 --bfile ${IN_BFILE} \
		--exclude $EXCLUDE_LIST \
		--make-bed \
		--out ${OUTDIR}/${NAMEBASE}_excluded
else
	cp ${IN_BFILE}.bed ${OUTDIR}/${NAMEBASE}_excluded.bed
	cp ${IN_BFILE}.bim ${OUTDIR}/${NAMEBASE}_excluded.bim
	cp ${IN_BFILE}.fam ${OUTDIR}/${NAMEBASE}_excluded.fam
fi

for CHR in {1..22}; do

plink2 --bfile ${OUTDIR}/${NAMEBASE}_excluded \
	--chr $CHR \
	--exclude range \
	${FILT_DIR}/LCR_filter_gnomAD_v3/filtered_position_hg38/${CHR}_hg38_LCR_gnomAD_v3.txt \
	${FILT_DIR}/MultiAllelic_and_InsDel_filter_gnomAD_v3/filtered_position_hg38/${CHR}_hg38_MultiAllelic_and_InsDel_gnomAD_v3.txt \
	${FILT_DIR}/PrimerPoly_5perc_filter_gnomAD_v3/filtered_position_hg38/${CHR}_hg38_PrimerPoly_5perc_gnomAD_v3.txt \
	${FILT_DIR}/SV_filter_gnomAD_v2/filtered_position_hg38/hg38_SV_gnomAD_v2.txt \
	--output-chr chrM \
	--recode vcf bgz id-paste=iid \
	--out ${OUTDIR}/${NAMEBASE}_excluded.chr${CHR}
	
bcftools annotate \
	--set-id '%CHROM\:%POS\:%REF\:%ALT' \
	${OUTDIR}/${NAMEBASE}_excluded.chr${CHR}.vcf.gz \
	| bgzip -c > ${OUTDIR}/${NAMEBASE}_excluded.chr${CHR}.SNPID.vcf.gz

rm ${OUTDIR}/${NAMEBASE}_excluded.chr${CHR}.vcf.gz

done

plink2 --bfile ${OUTDIR}/${NAMEBASE}_excluded \
	--chr X \
	--exclude range \
	${FILT_DIR}/LCR_filter_gnomAD_v3/filtered_position_hg38/X_hg38_LCR_gnomAD_v3.txt \
	${FILT_DIR}/MultiAllelic_and_InsDel_filter_gnomAD_v3/filtered_position_hg38/X_hg38_MultiAllelic_and_InsDel_gnomAD_v3.txt \
	${FILT_DIR}/PrimerPoly_5perc_filter_gnomAD_v3/filtered_position_hg38/X_hg38_PrimerPoly_5perc_gnomAD_v3.txt \
	${FILT_DIR}/SV_filter_gnomAD_v2/filtered_position_hg38/hg38_SV_gnomAD_v2.txt \
	--output-chr chrM \
	--make-bed \
	--out ${OUTDIR}/${NAMEBASE}_excluded.chrX
	
plink1.9 --bfile ${OUTDIR}/${NAMEBASE}_excluded.chrX \
	--set-hh-missing \
	--keep-allele-order \
	--allow-no-sex \
	--output-chr chrM \
	--recode vcf-iid bgz \
	--out ${OUTDIR}/${NAMEBASE}_excluded.chrX.hh
	
bcftools annotate \
	--set-id '%CHROM\:%POS\:%REF\:%ALT' \
	${OUTDIR}/${NAMEBASE}_excluded.chrX.hh.vcf.gz \
	| bgzip -c > ${OUTDIR}/${NAMEBASE}_excluded.chrX.hh.SNPID.vcf.gz
	
rm ${OUTDIR}/${NAMEBASE}_excluded.chrX.hh.vcf.gz
rm ${OUTDIR}/${NAMEBASE}_excluded.chrX.{bed,bim,fam}
rm ${OUTDIR}/${NAMEBASE}_excluded.{bed,bim,fam}