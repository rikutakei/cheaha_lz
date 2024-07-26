#!/bin/bash

module load Anaconda3

# Load in conda environment with bgenix:
source activate tools

# Make sure the directories to put LZs and LD data are there:
mkdir -p /data/scratch/$(whoami)/lz_outputs/{lz_plots,ld_files}

# Where to put the LD results:
OUT_DIR=/data/scratch/$(whoami)/lz_outputs/ld_files

# Function to calculate the LD information of a region from either UKBB or 1KGP
CHR=$1
START=$2
STOP=$3
SNP=$4
ANC=$5
DAT=$6

UKB_PATH=/data/project/merrimanlab/reference_files/ukbiobank
UKB_GENO=${UKB_PATH}/genotypes/imputed/
KGP_PATH=/data/project/merrimanlab/reference_files/1kgp_data/per_ancestry/${ANC}/vcf

# If UKBB, run UKBB code. Otherwise, 1KGP
if [[ $DAT == "UKB"* ]]
then
# UKBB:
# First adjust the chromosome so it has leading 0's for single-digit chromosome:
CHR_UKB=$(echo ${CHR} | awk '{printf "%02.0f", $1}' )
bgenix -g ${UKB_GENO}/imputed_chr${CHR}.bgen -i ${UKB_GENO}/imputed_chr${CHR}.bgen.bgi -incl-range ${CHR_UKB}:${START}-${STOP} -vcf | bcftools reheader -h ${UKB_PATH}/helper_files/new_header.txt | bcftools annotate --rename-chrs ${UKB_PATH}/helper_files/rename_contigs.txt | bgzip -c > ${OUT_DIR}/${SNP}_tmp.vcf.gz
plink --vcf ${OUT_DIR}/${SNP}_tmp.vcf.gz --make-bed --out ${OUT_DIR}/${SNP}_tmp
awk '{gsub(";.*", "", $2); print}' ${OUT_DIR}/${SNP}_tmp.bim | tr -s ' ' '\t' > ${OUT_DIR}/${SNP}_tmp.bim_tmp && mv ${OUT_DIR}/${SNP}_tmp.bim_tmp ${OUT_DIR}/${SNP}_tmp.bim
plink --bfile ${OUT_DIR}/${SNP}_tmp --allow-no-sex --snps-only --r2 --inter-chr --ld-snp ${SNP} --ld-window-r2 0 --out ${OUT_DIR}/UKBB_region_${CHR}.${START}-${STOP}_${SNP}
else
# 1KGP:
bcftools view --regions ${CHR}:${START}-${STOP} --output-type z --output-file ${OUT_DIR}/${SNP}_tmp.vcf.gz ${KGP_PATH}/1KGP_${ANC}_chr${CHR}.no_relatives.rsid.vcf.gz
plink --vcf ${OUT_DIR}/${SNP}_tmp.vcf.gz --allow-no-sex --snps-only --r2 --inter-chr --ld-snp ${SNP} --ld-window-r2 0 --out ${OUT_DIR}/1KGP_${ANC}_region_${CHR}.${START}-${STOP}_${SNP}
fi

rm ${OUT_DIR}/${SNP}_tmp.* ${OUT_DIR}/*_region_${CHR}.${START}-${STOP}_${SNP}.nosex

