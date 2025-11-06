wd <- getwd()

packrat::init("~/runs/eyu8/library/SKAT")
library(SKAT)
library(MetaSKAT)

setwd(wd)

require(methods)
require(Xmisc)
parser <- ArgumentParser$new()

cohort = "UKBB_caseproxy_subcontrol"

File.Bed <- paste0("SKAT/", cohort, ".bed")
File.Bim <- paste0("SKAT/", cohort, ".bim")
File.Fam <- paste0("SKAT/", cohort, ".fam")
File.SetID <- paste0("SKAT/UKBB.SETID")
File.SSD <- paste0("SKAT/", cohort, ".SSD")
File.Info <- paste0("SKAT/", cohort, ".info")
File.Mat <- paste0("SKAT/", cohort, ".mat")
File.SetInfo <- paste0("SKAT/", cohort, ".MInfo")
File.cov <- paste0("covar/covar_", cohort, ".txt")
File.Results.SKATO <- paste0("SKAT/", cohort, ".results.skato")
File.Results.BURDEN <- paste0("SKAT/", cohort, ".results.burden")

Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)
SSD.INFO <- Open_SSD(File.SSD, File.Info)
  
FAM <- Read_Plink_FAM_Cov(File.Fam,
                          File.cov,
                          Is.binary = TRUE,
                          cov_header = TRUE)
y <- FAM$Phenotype
Age <- FAM$Age
Sex <- FAM$Sex.y
N.Sample <- length(y)
obj <- SKAT_Null_Model(y ~ Sex + Age, out_type = "D")
  
out.skato <- SKATBinary_Robust.SSD.All(SSD.INFO, obj, method = "SKATO")  
write.table(
    out.skato$results,
    file = File.Results.SKATO,
    col.names = TRUE,
    row.names = FALSE
  )
print(paste0("SKATO of ", cohort, " is done successfully."))

out.skato.burden <- SKATBinary_Robust.SSD.All(SSD.INFO, obj, method = "Burden")
write.table(
    out.skato.burden$results,
    file = File.Results.BURDEN,
    col.names = TRUE,
    row.names = FALSE
  )
print(paste0("Burden of ", cohort, " is done successfully."))

re1 <- Generate_Meta_Files(obj,
                            File.Bed,
                            File.Bim,
                            File.SetID,
                            File.Mat,
                            File.SetInfo,
                            N.Sample)