#!/bin/bash

mkdir -p METASKAT
mkdir -p SKAT
mkdir -p logs
mkdir -p annotation

### COHORT: UKBB (!!! Contains a lot of variants and takes a while to process !!!)
# STEP 1: MULTI-ALLELIC VARIANT REMOVAL
bcftools view \
  # -r "chr5:110739007-110765157" \
  --threads 8 \
  /lustre03/project/6032411/neurobioinfo_ziv/liulang/SLC25A46/UKB.WGS.SLC25A46.vcf.gz \
  --max-alleles 2 \
  -l 1 \
  -Oz -o UKBB.vcf.gz


## STEP 2: QC SCRIPT FOR GQ, DP 95% MISSINGNESS
bash UKBB_QC.sh

bcftools view -i 'MAF<0.01' UKBB_GQ25_DP25_MISS95_filtered.vcf.gz -Oz -o UKBB_MAF_0.01.vcf.gz

## STEP 3: VEP ANNOTATION OF VARIANTS
#Reference
FASTA=/cvmfs/ref.mugqic/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa
DIR_CACHE=/lustre03/project/6004655/COMMUN/data/VEP
DIR_PLUGINS=${DIR_CACHE}/Plugins

#ClinVar custom annotation
CLINVAR_VCF=${DIR_CACHE}/clinvar/clinvar_20240407.GRCh38.vcf.gz
CLINVAR_FIELDS="ALLELEID,CLNDN,CLNREVSTAT,CLNSIG,CLNSIGCONF"

#Plugin data files
ALPHAMISSENSE_DATA=${DIR_PLUGINS}/AlphaMissense/AlphaMissense_hg38.tsv.gz
CADD_FILE=~/runs/sajanth/2025/CADD/GRCh38/whole_genome_SNVs.tsv.gz
CADD_INDEL_FILE=~/runs/sajanth/2025/CADD/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz

srun --account=def-grouleau -c 4 --mem=120g -t "3:00:00" vep -i UKBB_MAF_0.01.vcf.gz -o UKBB.VEPannotated --force_overwrite --safe \
    --cache --cache_version 111 --buffer_size 100 --dir_cache ${DIR_CACHE} \
    --format vcf --offline  --fasta ${FASTA} --xref_refseq --assembly GRCh38 \
    --pick_allele \
    --show_ref_allele --af_gnomadg --hgvs --symbol \
    --custom ${CLINVAR_VCF},ClinVar,vcf,exact,0,${CLINVAR_FIELDS} \
    --dir_plugins ${DIR_PLUGINS} --plugin AlphaMissense,file=${ALPHAMISSENSE_DATA} \
    --plugin CADD,${CADD_FILE},${CADD_INDEL_FILE}


## STEP 4: PREPARING THE SKAT-O FILES .bed .bim .fam

# case + control
# plink --vcf UKBB_MAF_0.01.vcf.gz --max-maf 0.01 --vcf-half-call m --output-chr M --make-bed --out SKAT/UKBB
# plink --bfile SKAT/UKBB  --update-sex sex_UKBB.txt --pheno covar_UKBB.txt --pheno-name Status --make-bed --allow-no-sex --keep covar_UKBB.txt --out SKAT/UKBB --output-chr M 

plink --vcf UKBB_MAF_0.01.vcf.gz \
  # --max-maf 0.01 \
  # Treat heterozygous missing genotypes (e.g., 0/. or 1/.) as missing (. / .)
  --vcf-half-call m \ 
  # Output chromosomes as integers (1, 2, ..., X, Y), without the chr prefix. M refers to "numeric representation."
  --output-chr M \
  # Create PLINK binary format files: .bed, .bim, .fam
  --make-bed \
  # Save the resulting files in the SKAT/ directory with prefix UKBB
  --out SKAT/UKBB

# case + proxy + control
plink --bfile SKAT/UKBB \
  --update-sex covar/sex_UKBB_caseproxy_control_noheader.txt \
  --pheno covar/covar_UKBB_caseproxy_control.txt \
  --pheno-name Status \
  --make-bed \
  --allow-no-sex \
  --keep covar/covar_UKBB_caseproxy_control.txt \
  --out SKAT/UKBB_caseproxy_control \
  --output-chr M 

## STEP 5: PREPARING THE INPUT FOR VARIANT-LEVEL ASSOCIATION ANALYSIS - FISHER TEST, MAF, AND COUNTS PER VARIANT CASE/CONTROL

# for cohort in "UKBB" "AMP_PD"; do

#  plink --bfile SKAT/${cohort} --assoc fisher --counts --freq --out annotation/${cohort} --output-chr M
#  plink --bfile SKAT/${cohort} --assoc fisher --out annotation/${cohort}_MAF --output-chr M

#done

plink --bfile SKAT/UKBB_caseproxy_control --assoc fisher --counts --freq --out annotation/UKBB_caseproxy_control --output-chr M
plink --bfile SKAT/UKBB_caseproxy_control --assoc fisher --out annotation/UKBB_caseproxy_control_MAF --output-chr M
