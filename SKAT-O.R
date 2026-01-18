#!/usr/bin/env Rscript
# ================================================================
# SKAT-O rare variant association analysis runner
#
# This script performs gene- or region-based rare variant
# association testing using SKAT-O and Burden methods based on
# PLINK genotype data and a predefined variant set definition.
#
# Author: Han Yu
#
# Input files:
#   - PLINK files: .bed / .bim / .fam
#   - setID file: Variant-to-gene (or region) mapping
#   - covariate file: Sample-level covariates
#
# Variables (edit below as needed):
#   File.Bed            : PLINK .bed file
#   File.Bim            : PLINK .bim file
#   File.Fam            : PLINK .fam file
#   File.SetID          : Variant set definition file
#   File.SSD            : Output SSD file
#   File.Info           : Output SSD info file
#   File.Mat            : (Reserved) genotype matrix file
#   File.SetInfo        : (Reserved) set annotation file
#   File.cov            : Covariate file
#   File.Results.SKATO  : Output file for SKAT-O results
#   File.Results.BURDEN : Output file for Burden test results
#   File.Results.FDR    : Output file for FDR correction results
#
# ================================================================

# Load library
library(SKAT)
library(Matrix)
library(data.table)
library(readr)
library(dplyr)

# Input arguments
File.Bed <- $1
File.Bim <- $2
File.Fam <- $3
File.SetID <- $4
File.SSD <- $5
File.Info <- $6
File.Mat <- $7
File.SetInfo <- $8
File.cov <- $9
File.Results.SKATO <- $10
File.Results.BURDEN <- $11
File.Results.FDR <- $12

# Create SSD file
Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)
SSD.INFO <- Open_SSD(File.SSD, File.Info)
 
# Load covariates and phenotype
FAM <- Read_Plink_FAM_Cov(File.Fam,
                          File.cov,
                          Is.binary = TRUE,
                          cov_header = TRUE)
y <- FAM$Phenotype
Age <- FAM$Age
Sex <- FAM$Sex.y
PC1<-FAM$PC1
PC2<-FAM$PC2
PC3<-FAM$PC3
PC4<-FAM$PC4
PC5<-FAM$PC5

# Create null model
obj <- SKAT_Null_Model(y ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5, out_type = "D")

# Run SKAT-O 
out.skato <- SKATBinary_Robust.SSD.All(SSD.INFO, obj, method = "SKATO")  
write.table(
    out.skato$results,
    file = File.Results.SKATO,
    col.names = TRUE,
    row.names = FALSE
  )
print(paste0("SKATO of ", cohort, " is done successfully."))

# Run Burden test
out.skato.burden <- SKATBinary_Robust.SSD.All(SSD.INFO, obj, method = "Burden")
write.table(
    out.skato.burden$results,
    file = File.Results.BURDEN,
    col.names = TRUE,
    row.names = FALSE
  )
print(paste0("Burden of ", cohort, " is done successfully."))

Close_SSD(File.SSD)

# FDR correction
result_FDR <- out.skato$results %>% mutate(p.fdr = p.adjust(P.value, method = "fdr"))
write.table(result_FDR, file = File.Results.FDR, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


# End of script