#!/bin/bash
# This script is for processing the data from the TOPMed server after imputation.
OUTDIR=$1
NAMEBASE=$2

# Decompress imputation output
parallel -j 23 "unzip \
	-P $3 \
	${OUTDIR}/chr_{}.zip" \
	::: {1..22} X

# Filter for Rsq > 0.7
parallel -j 23 "bcftools view \
	-i 'INFO/R2 > 0.7' \
	-Oz \
	-o ${OUTDIR}/chr{1}_r20.7.dose.vcf.gz \
	${OUTDIR}/chr{1}.dose.vcf.gz \
	&& rm ${OUTDIR}/chr{1}.dose.vcf.gz" \
	::: {1..22} X

# Add chr:pos:ref:alt SNPIDs
parallel -j 23 "bcftools annotate \
	--set-id '%CHROM:%POS:%REF:%ALT' \
	-Oz \
	-o ${OUTDIR}/chr{1}_r20.7_annSNPID.vcf.gz \
	${OUTDIR}/chr{1}_r20.7.dose.vcf.gz \
	&& rm ${OUTDIR}/chr{1}_r20.7.dose.vcf.gz" \
	::: {1..22} X

parallel -j 23 "bcftools head \
	${OUTDIR}/chr{1}_r20.7_annSNPID.vcf.gz \
	| grep -v contig \
	> ${OUTDIR}/chr{1}_r20.7_annSNPID.hdr" \
	::: {1..22} X
	
parallel -j 23 "bcftools reheader \
	--header ${OUTDIR}/chr{1}_r20.7_annSNPID.hdr \
	-o ${OUTDIR}/chr{1}_r20.7_annSNPID.nocontig.vcf.gz \
	${OUTDIR}/chr{1}_r20.7_annSNPID.vcf.gz \
	&& rm ${OUTDIR}/chr{1}_r20.7_annSNPID.vcf.gz \
	&& rm ${OUTDIR}/chr{1}_r20.7_annSNPID.hdr" \
	::: {1..22} X

# Add contigs to make headers compatible
parallel -j 23 "bcftools reheader \
	--fai /03-DryLab/02-Data/01-Array/02-Processed/2024_imputation/03_reference_files/GRCh38.fasta.fai \
	-o ${OUTDIR}/chr{1}_r20.7_annSNPID.hdr.vcf.gz \
	${OUTDIR}/chr{1}_r20.7_annSNPID.nocontig.vcf.gz \
	&& rm ${OUTDIR}/chr{1}_r20.7_annSNPID.nocontig.vcf.gz" \
	::: {1..22} X

# Gather vcfs
INPUTS=()
for CHR in {1..22} X; do
	INPUTS+=" -I ${OUTDIR}/chr${CHR}_r20.7_annSNPID.hdr.vcf.gz"
done

gatk GatherVcfs \
	${INPUTS[@]} \
	-O ${OUTDIR}/${NAMEBASE}_gathered.vcf.gz \
	&& rm ${OUTDIR}/chr*_r20.7_annSNPID.hdr.vcf.gz

# Turn into plink files
plink2 --vcf ${OUTDIR}/${NAMEBASE}_gathered.vcf.gz \
	--double-id \
	--make-bed \
	--out ${OUTDIR}/00-${NAMEBASE}_annSNPID_r20.7 \
&& rm ${OUTDIR}/${NAMEBASE}_gathered.vcf.gz
	
# Clean up files
mkdir ${OUTDIR}/imputation_logs
mv ${OUTDIR}/chr*.log ${OUTDIR}/imputation_logs
mkdir ${OUTDIR}/imputation_zips
mv ${OUTDIR}/chr*.zip ${OUTDIR}/imputation_zips
mkdir ${OUTDIR}/info_files
mv ${OUTDIR}/chr*.info.gz ${OUTDIR}/info_files