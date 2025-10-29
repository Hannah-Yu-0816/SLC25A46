module load StdEnv/2020
module load gatk/4.2.5.0

mkdir -p logs  # Ensure the logs directory exists

########### step by step ###########

srun -c 1 --mem=18g -t "3:00:00" gatk VariantFiltration \
  -R /lustre03/project/6004655/COMMUN/runs/senkkon/UKBB_DP/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  -V UKBB.vcf.gz \
  -O UKBB_GQ25.vcf.gz \
  --java-options "-Xmx12g" \
  --genotype-filter-expression "GQ < 25" \
  --genotype-filter-name "GQ25" \
  --set-filtered-genotype-to-no-call \
  --missing-values-evaluate-as-failing 2>&1 | tee logs/UKBB_GQ25.log &&

srun -c 1 --mem=18g -t "3:00:00" gatk VariantFiltration \
  -R /lustre03/project/6004655/COMMUN/runs/senkkon/UKBB_DP/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  -V UKBB_GQ25.vcf.gz \
  -O UKBB_GQ25_DP25.vcf.gz \
  --java-options "-Xmx12g" \
  --genotype-filter-expression "DP < 25" \
  --genotype-filter-name "DP25" \
  --set-filtered-genotype-to-no-call \
  --missing-values-evaluate-as-failing 2>&1 | tee logs/UKBB_GQ25_DP25.log &&

#Total AN is 981,098 (sample size of the file = 490,549 x 2 = AN) (For 5% missingness threshold we set AN < 932043)
srun -c 1 --mem=18g -t "3:00:00" gatk VariantFiltration \
  -R /lustre03/project/6004655/COMMUN/runs/senkkon/UKBB_DP/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  -V UKBB_GQ25_DP25.vcf.gz \
  -O UKBB_GQ25_DP25_MISS95.vcf.gz \
  --java-options "-Xmx12g" \
  --filter-expression "AN < 932043" \
  --filter-name "MISS" \
  --missing-values-evaluate-as-failing 2>&1 | tee logs/UKBB_GQ25_DP25_MISS95.log &&

srun -c 1 --mem=18g -t "3:00:00" gatk SelectVariants \
  -R /lustre03/project/6004655/COMMUN/runs/senkkon/UKBB_DP/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  -V UKBB_GQ25_DP25_MISS95.vcf.gz \
  -O UKBB_GQ25_DP25_MISS95_filtered.vcf.gz \
  --java-options "-Xmx12g" \
  --exclude-filtered \
  --exclude-non-variants 2>&1 | tee logs/UKBB_GQ25_DP25_MISS95_filtered.log 

####################################

srun --account=def-grouleau -c 1 --mem=18g -t "8:00:00" gatk VariantFiltration \
  -R /lustre03/project/6004655/COMMUN/runs/senkkon/UKBB_DP/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  -V UKBB.vcf.gz \
  -O UKBB_GQ25_DP25_MISS95.vcf.gz \
  --java-options "-Xmx12g" \
  --genotype-filter-expression "GQ < 25" \
  --genotype-filter-name "GQ25" \
  --genotype-filter-expression "DP < 25" \
  --genotype-filter-name "DP25" \
  --set-filtered-genotype-to-no-call \
  --filter-expression "AN < 932043" \
  --filter-name "MISS" \
  --missing-values-evaluate-as-failing 2>&1 | tee logs/UKBB_GQ25_DP25_MISS95.log &&

srun --account=def-grouleau -c 1 --mem=18g -t "5:00:00" gatk SelectVariants \
  -R /lustre03/project/6004655/COMMUN/runs/senkkon/UKBB_DP/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  -V UKBB_GQ25_DP25_MISS95.vcf.gz \
  -O UKBB_GQ25_DP25_MISS95_filtered.vcf.gz \
  --java-options "-Xmx12g" \
  --exclude-filtered \
  --exclude-non-variants 2>&1 | tee logs/UKBB_GQ25_DP25_MISS95_filtered.log 
