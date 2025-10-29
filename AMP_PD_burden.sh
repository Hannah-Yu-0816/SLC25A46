#!/bin/bash

mkdir -p METASKAT
mkdir -p SKAT
mkdir -p logs
mkdir -p annotation

### COHORT: AMP-PD
#FILTRATION (REMOVES PPMI AND NON-EUROPEAN ANCESTRY)
plink --bfile ~/runs/senkkon/2022/TIM_TOM/AMP_PD/AMP_PD_FILTERED_ALL_CHR --extract range SLC25A46.bed --keep covar_AMP_PD.txt --max-maf 0.01 --make-bed --out SKAT/AMP_PD
plink --bfile SKAT/AMP_PD --recode vcf bgz --out AMP_PD

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

srun -c 4 --mem=120g -t "1:00:00" vep -i AMP_PD.vcf.gz -o AMP_PD.VEPannotated --force_overwrite --safe \
    --cache --cache_version 111 --buffer_size 100 --dir_cache ${DIR_CACHE} \
    --format vcf --offline  --fasta ${FASTA} --xref_refseq --assembly GRCh38 \
    --pick_allele \
    --show_ref_allele --af_gnomadg --hgvs --symbol \
    --custom ${CLINVAR_VCF},ClinVar,vcf,exact,0,${CLINVAR_FIELDS} \
        --dir_plugins ${DIR_PLUGINS} --plugin AlphaMissense,file=${ALPHAMISSENSE_DATA} \
        --plugin CADD,${CADD_FILE},${CADD_INDEL_FILE}


##FREQ AND SAMPLE SIZES
for cohort in "AMP_PD"; do

	plink --bfile SKAT/${cohort} --assoc fisher --counts --freq --out annotation/${cohort} --output-chr M
	plink --bfile SKAT/${cohort} --assoc fisher --out annotation/${cohort}_MAF --output-chr M

done
