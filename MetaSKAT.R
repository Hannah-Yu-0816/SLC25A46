#!/usr/bin/env Rscript
# ================================================================
# MetaSKAT rare variant association meta-analysis runner
#
# This script performs cohort-level SKAT-O and Burden tests
# followed by meta-analysis across cohorts using MetaSKAT,
# based on summary statistics generated from SKAT.
#
# Author: Han Yu
#
# Input files:
#   - PLINK files: .bed / .bim / .fam (per cohort)
#   - setID file: Variant-to-gene (or region) mapping
#   - covariate files: Sample-level covariates
#
# Variables (edit below as needed):
#   Base.SKAT.Dir     : Base directory containing SKAT-related input and output files
#   Base.COV.Dir      : Directory containing covariate files
#   Cohort.File       : File listing cohort identifiers (one per line)
#   File.Meta.SKATO   : Output file for MetaSKAT-O results
#   File.Meta.BURDEN  : Output file for MetaSKAT Burden results
#
# ================================================================

# Load libraries
library(SKAT)
library(MetaSKAT)

# ---------------------------
# Input arguments
# ---------------------------
Base.SKAT.Dir    <- $1
Base.COV.Dir     <- $2
Cohort.File      <- $3
File.Meta.SKATO  <- $4
File.Meta.BURDEN <- $5

# ---------------------------
# Cohorts
# ---------------------------
all_cohorts <- scan(Cohort.File, what = character())
filtered_cohorts <- c()

# Identify available cohorts
for (cohort in all_cohorts) {

  File.SetID <- paste0(Base.SKAT.Dir, "/", cohort, ".SETID")

  if (file.exists(File.SetID)) {
    filtered_cohorts <- c(filtered_cohorts, cohort)
  } else {
    print(paste0(cohort, " is excluded as it does not have the SETID file."))
  }
}

File.Mat.vec <- rep("", length(filtered_cohorts))
File.SetInfo.vec <- rep("", length(filtered_cohorts))

# ---------------------------
# Cohort-level analysis
# ---------------------------
for (i in seq_along(filtered_cohorts)) {

  cohort <- filtered_cohorts[i]

  File.Bed <- paste0(Base.SKAT.Dir, "/", cohort, ".bed")
  File.Bim <- paste0(Base.SKAT.Dir, "/", cohort, ".bim")
  File.Fam <- paste0(Base.SKAT.Dir, "/", cohort, ".fam")
  File.SetID <- paste0(Base.SKAT.Dir, "/", cohort, ".SETID")
  File.SSD <- paste0(Base.SKAT.Dir, "/", cohort, ".SSD")
  File.Info <- paste0(Base.SKAT.Dir, "/", cohort, ".info")
  File.Mat <- paste0(Base.SKAT.Dir, "/", cohort, ".mat")
  File.SetInfo <- paste0(Base.SKAT.Dir, "/", cohort, ".MInfo")
  File.cov <- paste0(Base.COV.Dir, "/", cohort, ".txt")

  File.Results.SKATO <- paste0(Base.SKAT.Dir, "/", cohort, ".results.skato")
  File.Results.BURDEN <- paste0(Base.SKAT.Dir, "/", cohort, ".results.burden")

  # Create SSD file
  Generate_SSD_SetID(
    File.Bed,
    File.Bim,
    File.Fam,
    File.SetID,
    File.SSD,
    File.Info
  )

  SSD.INFO <- Open_SSD(File.SSD, File.Info)

  # Load covariates and phenotype
  FAM <- Read_Plink_FAM_Cov(
    File.Fam,
    File.cov,
    Is.binary = TRUE,
    cov_header = TRUE
  )

  y <- FAM$Phenotype
  Age <- FAM$Age
  Sex <- FAM$Sex.y
  PC1<-FAM$PC1
  PC2<-FAM$PC2
  PC3<-FAM$PC3
  PC4<-FAM$PC4
  PC5<-FAM$PC5
  N.Sample <- length(y)

  # Create null model
  obj <- SKAT_Null_Model(y ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5, out_type = "D")

  # Run SKAT-O
  out.skato <- SKATBinary_Robust.SSD.All(
    SSD.INFO,
    obj,
    method = "SKATO"
  )

  write.table(
    out.skato$results,
    file = File.Results.SKATO,
    col.names = TRUE,
    row.names = FALSE
  )

  print(paste0("SKATO of ", cohort, " is done successfully."))

  # Run Burden test
  out.skato.burden <- SKATBinary_Robust.SSD.All(
    SSD.INFO,
    obj,
    method = "Burden"
  )

  write.table(
    out.skato.burden$results,
    file = File.Results.BURDEN,
    col.names = TRUE,
    row.names = FALSE
  )

  print(paste0("Burden of ", cohort, " is done successfully."))

  # Generate MetaSKAT input files
  Generate_Meta_Files(
    obj,
    File.Bed,
    File.Bim,
    File.SetID,
    File.Mat,
    File.SetInfo,
    N.Sample
  )

  File.Mat.vec[i] <- File.Mat
  File.SetInfo.vec[i] <- File.SetInfo
}

# ---------------------------
# Meta-analysis
# ---------------------------
Cohort.Info <- Open_MSSD_File_2Read(
  File.Mat.vec,
  File.SetInfo.vec
)

print("Cohort.Info generated successfully.")

out.skato.burden <- MetaSKAT_MSSD_ALL(
  Cohort.Info,
  method = "Burden"
)

write.table(
  out.skato.burden,
  file = File.Meta.BURDEN,
  col.names = TRUE,
  row.names = FALSE
)

print("Meta Burden is done successfully.")

out.skato <- MetaSKAT_MSSD_ALL(
  Cohort.Info,
  method = "optimal"
)

write.table(
  out.skato,
  file = File.Meta.SKATO,
  col.names = TRUE,
  row.names = FALSE
)

print("Meta SKATO is done successfully.")


# End of scripts