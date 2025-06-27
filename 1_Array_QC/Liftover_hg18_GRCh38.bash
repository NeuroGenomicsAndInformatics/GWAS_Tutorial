#!/bin/bash
# This script lifts over plink files from hg19 to GRCh38 using the UCSC liftover tool

IN_BFILE=$1
OUTDIR=$2
NAMEBASE=$3
INDEX_NUM=$(echo ${IN_BFILE##*/} | cut -d- -f1)
printf -v TEMP "%s%02d%s" ${OUTDIR}/ $((10#${INDEX_NUM} + 1)) -${NAMEBASE}
printf -v TEMP2 "%s%02d%s" ${OUTDIR}/ $((10#${INDEX_NUM} + 2)) -${NAMEBASE}
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Update chromosome codes
plink1.9 --bfile ${IN_BFILE} \
	--output-chr chrM \
	--keep-allele-order \
	--allow-no-sex \
	--make-bed \
	--out ${TEMP}.chrs

# Turn bim into bed file
awk '{OFS="\t"}; {print $1,$4-1,$4,$2}' ${TEMP}.chrs.bim > ${TEMP}.chrs.tolift.bed

# Create lifted files
${SCRIPT_DIR}/liftOver \
	${TEMP}.chrs.tolift.bed \
	${SCRIPT_DIR}/hg18ToHg38.over.chain.gz \
	${TEMP}.chrs.tolift.bed-hg38 ${TEMP}.chrs.tolift.bed-unmapped

# Remove variants on weird chromosomes
grep -v  '#Deleted' ${TEMP}.chrs.tolift.bed-unmapped \
	| awk '{print $4}' > ${OUTDIR}/to_remove.txt

grep '_alt'  ${TEMP}.chrs.tolift.bed-hg38 \
	| awk '{print $2}'  >> ${OUTDIR}/to_remove.txt

grep '_random'  ${TEMP}.chrs.tolift.bed-hg38 \
	| awk '{print $2}' >> ${OUTDIR}/to_remove.txt

grep 'chrUn_'  ${TEMP}.chrs.tolift.bed-hg38 \
	| awk '{print $2}' >> ${OUTDIR}/to_remove.txt

#Creating files to update plink file to hg38
grep -Ev '_alt|_random|chrUn_' ${TEMP}.chrs.tolift.bed-hg38 > ${TEMP}.chrs.tolift.bed-hg38-clean

awk '{print $4, $3}'  ${TEMP}.chrs.tolift.bed-hg38-clean > ${OUTDIR}/update_pos.txt

awk '{print $4, $1}'  ${TEMP}.chrs.tolift.bed-hg38-clean > ${OUTDIR}/update_chr.txt

plink1.9 --bfile ${TEMP}.chrs \
	--exclude ${OUTDIR}/to_remove.txt \
	--keep-allele-order \
	--allow-no-sex \
	--update-chr ${OUTDIR}/update_chr.txt \
	--update-map ${OUTDIR}/update_pos.txt \
	--make-bed \
	--out ${TEMP2}

# Clean up
rm ${TEMP}.chrs.{bed,bim,fam}
