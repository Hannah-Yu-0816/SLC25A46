# Get cleaned bfiles
cp ~/runs/go_lab/GRCh38/RBD_WGS/cleaned_chr5.* .
gunzip cleaned_chr5.*

# Filter out common variants
plink --bfile cleaned_chr5 --maf 0.01 --make-bed --out common
awk '{print $2}' common.bim > common.txt
plink --bfile cleaned_chr5 --exclude common.txt --make-bed --out cleaned_chr5_rare

# Get target Gene position
chr=5
start=110738707   # 110739007 - 300
end=110765457     # 110765157 + 300

plink --bfile cleaned_chr5_rare --chr "$chr" --from-bp "$start" --to-bp "$end" --make-bed --out chr5_temp_SLC25A46

# VEP
plink --bfile chr5_temp_SLC25A46 --recode vcf --out RBD_WGS

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

srun --account=def-grouleau -c 4 --mem=120g -t "1:00:00" vep -i RBD_WGS.vcf -o RBD.VEPannotated --force_overwrite --safe \
    --cache --cache_version 111 --buffer_size 100 --dir_cache ${DIR_CACHE} \
    --format vcf --offline  --fasta ${FASTA} --xref_refseq --assembly GRCh38 \
    --pick_allele \
    --show_ref_allele --af_gnomadg --hgvs --symbol \
    --custom ${CLINVAR_VCF},ClinVar,vcf,exact,0,${CLINVAR_FIELDS} \
        --dir_plugins ${DIR_PLUGINS} --plugin AlphaMissense,file=${ALPHAMISSENSE_DATA} \
        --plugin CADD,${CADD_FILE},${CADD_INDEL_FILE}


#PREPARING THE SKAT-O FILES .bed .bim .fam

#plink --vcf RBD_WGS.vcf --max-maf 0.01 --vcf-half-call m --output-chr M --make-bed --out ../SKAT/RBD
#plink --bfile ../SKAT/RBD  --update-sex ../covar/sex_RBD.txt --make-bed --allow-no-sex --keep ../covar/RBD_AMP-PD_covar.txt --out ../SKAT/RBD --output-chr M 


## PREPARING THE INPUT FOR VARIANT-LEVEL ASSOCIATION ANALYSIS - FISHER TEST, MAF, AND COUNTS PER VARIANT CASE/CONTROL


#plink --bfile ../SKAT/RBD --assoc fisher --counts --freq --out ../annotation/RBD --output-chr M
#plink --bfile ../SKAT/RBD --assoc fisher --out ../annotation/RBD_MAF --output-chr M
