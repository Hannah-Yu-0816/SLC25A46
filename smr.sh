#!/usr/bin/env bash
# ================================================================
# Summary-based Mendelian Randomization (SMR) runner
# 
# This script performs SMR analyses across multiple cis-window 
# definitions using summary-level GWAS and QTL data.
# 
# Author: Han Yu
#
# Input files:
#   - gene_bed: Tab-separated file with columns: chr, start, end, gene
#   - wind_list: File containing window sizes (one per line)
#
# Variables (edit below as needed):
#   gene_bed    : Path to the bed file containing genes info
#   wind_list   : Path to the file containing window sizes (kb)
#   smr_final   : Output for SMR output files
#   gwas        : Path to GWAS summary statistics
#   LD_ref      : Path to LD reference panel prefix
#   qQTL_file   : Path to eQTL summary files
#
# ================================================================

# ---------------------------
# Input arguments
# ---------------------------
gene_bed=$1
wind_list=$2
smr_final=$3
gwas=$4
LD_ref=$5
qQTL_file=$6

# ---------------------------
# Main loop
# ---------------------------
# Read gene list line by line
while IFS=$'\t' read -r chr position_start position_end gene; do

    # Read each window size from the window list
    while IFS= read -r wind; do

        # Run SMR for the current gene and window
        ./smr \
         --bfile "$LD_ref" \
         --gwas-summary "$gwas" \
         --beqtl-summary "$qQTL_file" \
         --chr "$chr" \
         --gene "$gene" \
         --peqtl-smr 5e-8 \
         --smr-multi --cis-wind "$wind" \
         --out "${smr_final}smr_${gene}_chr${chr}_wind${wind}"

    done < "$wind_list"

done < "$gene_bed"

echo "SMR analysis completed."

# End of script
