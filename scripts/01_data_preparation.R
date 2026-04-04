# ============================================================
# PROJECT   : Multiomics integration - TCGA-COAD
# SCRIPT    : 01 - Data preparation
# GOAL      : Download, clean and align RNA-seq + methylation
#             data from TCGA for MOFA+ integration
# AUTHOR    : Fadwa EL KHADDAR
# ============================================================

set.seed(42)

# 0. Libraries 
library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)
library(DESeq2)
library(minfi)
library(matrixStats)
library(tidyverse)


# 1. Data download 
# Browse available COAD assays without downloading
curatedTCGAData("COAD", assays = "*", dry.run = TRUE, version = "2.0.1")

# Download RNA-seq and methylation only (~500 Mb total)
# RNASeq2GeneNorm      : already RSEM-normalized by TCGA
# Methylation_methyl450 : beta values (0-1), 485k CpG probes
coad <- curatedTCGAData(
  diseaseCode = "COAD",
  assays      = c("RNASeq2GeneNorm_illuminahiseq", "Methylation_methyl450"),
  dry.run     = FALSE,
  version     = "2.0.1"
)

coad


# 2. Clinical data 
# Extract all clinical variables (1615 columns total)
clinical_raw <- as.data.frame(colData(coad))

# WHY OVERALL SURVIVAL ?
# MSI status was available for only 90/313 patients -> too few
# Overall survival uses all 308 patients and is the key
# clinical endpoint in oncology

# Build a clean survival table with only 4 useful columns
clinical_surv <- clinical_raw %>%
  dplyr::mutate(
    
    # Shorten TCGA barcode to 12 chars (patient-level ID)
    # e.g. "TCGA-A6-2671-01A-..." -> "TCGA-A6-2671"
    patient_id = substr(rownames(clinical_raw), 1, 12),
    
    # Vital status : 1 = deceased, 0 = censored (still alive)
    status = dplyr::case_when(
      vital_status.x == "1" ~ 1,
      vital_status.x == "0" ~ 0,
      TRUE                  ~ NA_real_
    ),
    
    # Follow-up time in days
    # If deceased : days from diagnosis to death
    # If censored : days from diagnosis to last contact
    time_days = dplyr::case_when(
      status == 1 ~ as.numeric(days_to_death.x),
      status == 0 ~ as.numeric(days_to_last_followup.x),
      TRUE        ~ NA_real_
    ),
    
    # Tumor stage (binary) - used later as secondary validation
    # early    = stage I / II
    # advanced = stage III / IV
    stage_bin = dplyr::case_when(
      pathologic_stage %in%
        c("stage i", "stage ia", "stage ii",
          "stage iia", "stage iib", "stage iic")      ~ "early",
      pathologic_stage %in%
        c("stage iii", "stage iiia", "stage iiib",
          "stage iiic", "stage iv", "stage iva",
          "stage ivb")                                 ~ "advanced",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(patient_id, status, time_days, stage_bin) %>%
  dplyr::filter(!is.na(status), !is.na(time_days), time_days > 0)

# Summary
cat("Patients with complete survival data :", nrow(clinical_surv), "\n")
cat("Vital status - 0=alive / 1=deceased :\n")
print(table(clinical_surv$status))
cat("Tumor stage :\n")
print(table(clinical_surv$stage_bin, useNA = "always"))
cat("Median follow-up (days) :", median(clinical_surv$time_days), "\n")
cat("Median follow-up (years):",
    round(median(clinical_surv$time_days) / 365, 1), "\n")


# 3. Extract omics matrices 
rna_mat  <- assay(coad[["COAD_RNASeq2GeneNorm_illuminahiseq-20160128"]])
meth_mat <- assay(coad[["COAD_Methylation_methyl450-20160128"]])

cat("RNA-seq raw :", nrow(rna_mat),  "genes x", ncol(rna_mat),  "samples\n")
cat("Methylation :", nrow(meth_mat), "CpG x",   ncol(meth_mat), "samples\n")


# 4. Keep primary tumors only (tissue code = 01) 
# TCGA barcodes encode the tissue type at positions 14-15 :
#   01 = primary solid tumor   -> keep
#   11 = solid tissue normal   -> remove
#   02 = recurrent tumor       -> remove
#   06 = metastasis            -> remove
# Keeping only code 01 ensures a biologically homogeneous cohort

keep_primary_tumor <- function(mat) {
  tissue_code      <- substr(colnames(mat), 14, 15)
  cat("Tissue codes found :", unique(tissue_code), "\n")
  mat_01           <- mat[, tissue_code == "01", drop = FALSE]
  colnames(mat_01) <- substr(colnames(mat_01), 1, 12)
  mat_01           <- mat_01[, !duplicated(colnames(mat_01)), drop = FALSE]
  return(mat_01)
}

rna_tumor  <- keep_primary_tumor(rna_mat)
meth_tumor <- keep_primary_tumor(meth_mat)

cat("RNA  primary tumors :", ncol(rna_tumor),  "patients\n")
cat("Meth primary tumors :", ncol(meth_tumor), "patients\n")


# 5. Align the 3 sources on common patients
# Only patients present in RNA + methylation + clinical are kept
# Critical : MOFA+ requires identical patient order in all modalities

overlap_final <- Reduce(intersect, list(
  colnames(rna_tumor),
  colnames(meth_tumor),
  clinical_surv$patient_id
))
cat("Final overlap :", length(overlap_final), "patients\n")

rna_final  <- rna_tumor[,  overlap_final]
meth_final <- meth_tumor[, overlap_final]
clin_final <- clinical_surv[match(overlap_final, clinical_surv$patient_id), ]

# Sanity checks - both must be TRUE before proceeding
cat("Same order RNA / Methylation :",
    all(colnames(rna_final) == colnames(meth_final)), "\n")
cat("Same order RNA / Clinical    :",
    all(colnames(rna_final) == clin_final$patient_id), "\n")


# 6. Dimension reduction 
# Keep only the most variable features per modality
# Why : MOFA+ needs manageable input size for local RAM
#       Variable features carry the most biological signal
#       Stable features (low variance) are mostly noise

# RNA-seq : top 3000 most variable genes
vars_rna <- rowVars(rna_final)
top_rna  <- order(vars_rna, decreasing = TRUE)[1:3000]
rna_top  <- rna_final[top_rna, ]
cat("RNA reduced :", nrow(rna_top), "genes x", ncol(rna_top), "patients\n")

# Methylation : remove CpG with > 20% missing values first
na_frac    <- rowMeans(is.na(meth_final))
meth_clean <- meth_final[na_frac < 0.2, ]
cat("CpG after NA filter :", nrow(meth_clean), "\n")

# Methylation : top 3000 most variable CpG probes
vars_meth <- rowVars(meth_clean, na.rm = TRUE)
top_meth  <- order(vars_meth, decreasing = TRUE)[1:3000]
meth_top  <- meth_clean[top_meth, ]
cat("Methylation reduced :", nrow(meth_top), "CpG x", ncol(meth_top), "patients\n")


# 7. Save 
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

saveRDS(rna_top,    "data/processed/rna_top3000.rds")
saveRDS(meth_top,   "data/processed/meth_top3000.rds")
saveRDS(clin_final, "data/processed/clinical_surv.rds")

cat("\nAll files saved in data/processed/\n")
cat("\n=== FINAL DATASET SUMMARY ===\n")
cat("RNA-seq     :", nrow(rna_top),    "genes x", ncol(rna_top),    "patients\n")
cat("Methylation :", nrow(meth_top),   "CpG x",   ncol(meth_top),   "patients\n")
cat("Clinical    :", nrow(clin_final), "patients -",
    sum(clin_final$status == 1), "deceased /",
    sum(clin_final$status == 0), "alive\n")
cat("\nNext step : run 02_mofa_analysis.R\n")

# Reproducibility
sessionInfo() 