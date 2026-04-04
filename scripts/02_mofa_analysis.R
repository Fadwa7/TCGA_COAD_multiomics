# ============================================================
# PROJECT   : Multiomics integration - TCGA-COAD
# SCRIPT    : 02 - MOFA+ analysis and biological interpretation
# GOAL      : Identify multiomics factors associated with
#             colorectal cancer tumor stage progression
# AUTHOR    : Fadwa EL KHADDAR
# ============================================================


# 0. Libraries 
library(MOFA2)
library(survival)
library(survminer)
library(ggplot2)
library(gridExtra)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(gprofiler2)
library(tidyverse)


# 1. Load preprocessed data 
rna_top    <- readRDS("data/processed/rna_top3000.rds")
meth_top   <- readRDS("data/processed/meth_top3000.rds")
clin_final <- readRDS("data/processed/clinical_surv.rds")

cat("RNA     :", nrow(rna_top),    "genes x",    ncol(rna_top),    "patients\n")
cat("Meth    :", nrow(meth_top),   "CpG x",      ncol(meth_top),   "patients\n")
cat("Clinical:", nrow(clin_final), "patients\n")


# 2. Build and train MOFA+ model 
# Convert to standard R matrices (required by MOFA2)
rna_mat_mofa  <- as.matrix(rna_top)
meth_mat_mofa <- as.matrix(meth_top)

# Create MOFA object
mofa_data <- list(RNA = rna_mat_mofa, Methylation = meth_mat_mofa)
mofa      <- create_mofa(mofa_data)

# Model options : 15 latent factors
model_opts             <- get_default_model_options(mofa)
model_opts$num_factors <- 15

# Training options
train_opts                  <- get_default_training_options(mofa)
train_opts$convergence_mode <- "slow"
train_opts$seed             <- 42

# Data options
data_opts <- get_default_data_options(mofa)

# Train
mofa <- prepare_mofa(
  mofa,
  data_options     = data_opts,
  model_options    = model_opts,
  training_options = train_opts
)
mofa <- run_mofa(mofa, use_basilisk = TRUE)

saveRDS(mofa, "data/processed/mofa_model.rds")
cat("MOFA+ model trained and saved\n")


# 3. Variance explained 
# Per factor and per modality
plot_variance_explained(mofa, max_r2 = 15)
# Factor1  : dominant in methylation (~15%) - strong epigenetic
#            signal, not yet transcriptomic
# Factor5  : most interesting biologically - significant in both
#            modalities simultaneously, shared signal between
#            methylation and expression
# Factor7-15 : low variance explained - mostly noise

# Total variance per modality
plot_variance_explained(mofa, plot_total = TRUE)
# RNA-seq     : ~13% - correct for 3000 genes, transcriptomics
#               is noisier
# Methylation : ~47% - excellent, methylation is highly
#               structured in colorectal cancer


# 4. Extract factor scores and merge with clinical data 
factors_raw           <- get_factors(mofa, factors = "all")[[1]]
factors_df            <- as.data.frame(factors_raw)
factors_df$patient_id <- rownames(factors_df)

clin_mofa <- clin_final %>%
  dplyr::inner_join(factors_df, by = "patient_id")

cat("Patients in analysis :", nrow(clin_mofa), "\n")


# 5. Univariate Cox regression - overall survival 
# HR (Hazard Ratio) interpretation :
#   HR = 1.0 : no difference between groups
#   HR > 1.0 : high score group has higher death risk
#   HR < 1.0 : high score group has lower death risk (protective)
# Both HR and p-value are needed to draw conclusions

results_cox <- data.frame()

for (f in paste0("Factor", 1:15)) {
  formula  <- as.formula(paste("Surv(time_days, status) ~", f))
  cox_fit  <- coxph(formula, data = clin_mofa)
  cox_sum  <- summary(cox_fit)
  results_cox <- rbind(results_cox, data.frame(
    Factor  = f,
    HR      = round(cox_sum$coefficients[, "exp(coef)"], 3),
    p_value = round(cox_sum$coefficients[, "Pr(>|z|)"],  4)
  ))
}

results_cox <- results_cox[order(results_cox$p_value), ]
print(results_cox)

# RESULT : No factor reaches p < 0.05 with n=271 (74 events only)
# Factor4 (p=0.065, HR=0.924) shows the strongest trend :
# patients with high Factor4 score have 7.6% lower death risk
# This is common in multiomics - individual factors are rarely
# perfect predictors; biological interpretation matters more
# -> Switch to tumor stage validation (more statistical power)


# 6. Wilcoxon test - tumor stage association
# Biological question : which MOFA+ factors distinguish
# early (stage I-II) from advanced (stage III-IV) tumors ?
# Wilcoxon test : non-parametric, robust for MOFA+ scores

clin_stage <- clin_mofa %>%
  dplyr::filter(!is.na(stage_bin))

cat("Patients with known stage :", nrow(clin_stage), "\n")
print(table(clin_stage$stage_bin))

results_stage <- data.frame()

for (f in paste0("Factor", 1:15)) {
  # wilcox.test(values ~ groups) : test if values differ between groups
  test  <- wilcox.test(clin_stage[[f]] ~ clin_stage$stage_bin)
  means <- tapply(clin_stage[[f]], clin_stage$stage_bin, mean)
  results_stage <- rbind(results_stage, data.frame(
    Factor        = f,
    mean_early    = round(means["early"],    3),
    mean_advanced = round(means["advanced"], 3),
    p_value       = round(test$p.value,      4)
  ))
}

results_stage <- results_stage[order(results_stage$p_value), ]
print(results_stage)

# RESULT : 7 significant factors (p < 0.05)
# Higher in advanced tumors (increase with progression) :
#   Factor14, Factor3, Factor9, Factor6, Factor10
# Higher in early tumors (decrease with progression) :
#   Factor5 (p=0.003, most significant), Factor12
# Factor5 is the focus of biological interpretation :
#   high score in early tumors suggests an active mechanism
#   in well-differentiated cells that is lost during progression


# 7. Visualization - Factor5 and Factor6 boxplots 
p1 <- ggplot(clin_stage, aes(x = stage_bin, y = Factor5,
                             fill = stage_bin)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("early" = "#2ECC71",
                               "advanced" = "#E74C3C")) +
  labs(title    = "MOFA+ Factor5 by tumor stage",
       subtitle = "Wilcoxon p = 0.003",
       x = "Tumor stage", y = "Factor5 score") +
  theme_classic() +
  theme(legend.position = "none")

p2 <- ggplot(clin_stage, aes(x = stage_bin, y = Factor6,
                             fill = stage_bin)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("early" = "#2ECC71",
                               "advanced" = "#E74C3C")) +
  labs(title    = "MOFA+ Factor6 by tumor stage",
       subtitle = "Wilcoxon p = 0.032",
       x = "Tumor stage", y = "Factor6 score") +
  theme_classic() +
  theme(legend.position = "none")

grid.arrange(p1, p2, ncol = 2)

# Factor5 : clear separation, early tumors have positive median
#           score, advanced tumors have negative median score
#           -> captures a mechanism lost during tumor progression
# Factor6 : separation exists but less clear, high outliers in
#           advanced group -> heterogeneous signal, suggests a
#           molecular subgroup with specific epigenetic profile


# 8. Biological interpretation - Factor5 weights 
# Weights = contribution of each gene/CpG to the factor
# Positive weight : feature increases with Factor5 score (early)
# Negative weight : feature decreases with Factor5 score (advanced)

weights_rna  <- get_weights(mofa, views = "RNA",
                            factors = "Factor5")[[1]]
weights_meth <- get_weights(mofa, views = "Methylation",
                            factors = "Factor5")[[1]]

# Top 20 RNA genes per direction
top_rna_pos <- head(sort(weights_rna[, 1], decreasing = TRUE),  20)
top_rna_neg <- head(sort(weights_rna[, 1], decreasing = FALSE), 20)

cat("Top 20 genes HIGH in early tumors (positive weights):\n")
print(top_rna_pos)
# Result : RPL/RPS ribosomal genes -> active translation machinery
#          signature of well-differentiated functional cells

cat("\nTop 20 genes HIGH in advanced tumors (negative weights):\n")
print(top_rna_neg)
# Result : MYH9/MYH14 (motility), PRKDC (DNA repair),
#          RNF43 (Wnt pathway), HNF4A (differentiation)
#          -> dedifferentiation and invasion signature

# Top 20 CpG probes per direction
top_meth_pos <- head(sort(weights_meth[, 1], decreasing = TRUE),  20)
top_meth_neg <- head(sort(weights_meth[, 1], decreasing = FALSE), 20)

cat("\nTop 20 CpG hypermethylated in early tumors :\n")
print(top_meth_pos)

cat("\nTop 20 CpG hypermethylated in advanced tumors :\n")
print(top_meth_neg)


# 9. CpG annotation 
# Relation_to_Island : position relative to CpG island
#   Island  -> inside a CpG island
#   S_Shore -> south shore (up to 2kb)
#   N_Shore -> north shore (up to 2kb)
#   S_Shelf -> south shelf (2-4kb)
#   OpenSea -> far from any island
# Key rule : methylation in Island/Shore at promoter = gene silencing

# UCSC_RefGene_Group : position within the gene
#   TSS200  -> 200bp upstream of transcription start
#   TSS1500 -> 200-1500bp upstream
#   1stExon -> first exon
#   5'UTR   -> 5' untranslated region
#   Body    -> gene body
#   3'UTR   -> 3' untranslated region
# Key rule : methylation at TSS200/TSS1500/1stExon = transcription repression

anno_450k  <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
cols_keep  <- c("chr", "pos", "UCSC_RefGene_Name",
                "Relation_to_Island", "UCSC_RefGene_Group")

cpg_pos <- names(top_meth_pos)
cpg_neg <- names(top_meth_neg)

anno_pos <- as.data.frame(anno_450k[cpg_pos, cols_keep])
anno_neg <- as.data.frame(anno_450k[cpg_neg, cols_keep])

cat("CpG hypermethylated in EARLY tumors :\n")
print(anno_pos[, c("UCSC_RefGene_Name", "Relation_to_Island",
                   "UCSC_RefGene_Group")])
# Key genes : XKR6 (apoptosis), FZD10 (Wnt receptor),
#             ZSCAN18 (transcription factor)

cat("\nCpG hypermethylated in ADVANCED tumors :\n")
print(anno_neg[, c("UCSC_RefGene_Name", "Relation_to_Island",
                   "UCSC_RefGene_Group")])
# Key genes : NF1 (tumor suppressor), HOXD10 (differentiation),
#             TJP2 (tight junctions), BAG3 (apoptosis regulation)


# 0. Functional enrichment - gprofiler2 
# Goal : identify biological pathways enriched in Factor5 genes
# Threshold : 2 standard deviations (statistically justified)
# Sources : GO:BP (biological process), KEGG, Reactome
# Correction : FDR (False Discovery Rate) - standard in literature

seuil          <- 2 * sd(weights_rna[, 1])
genes_early    <- names(weights_rna[weights_rna[, 1] >  seuil, 1])
genes_advanced <- names(weights_rna[weights_rna[, 1] < -seuil, 1])

cat("Genes early (2SD)    :", length(genes_early), "\n")
cat("Genes advanced (2SD) :", length(genes_advanced), "\n")

gost_early <- gost(
  query             = genes_early,
  organism          = "hsapiens",
  sources           = c("GO:BP", "KEGG", "REAC"),
  significant       = TRUE,
  correction_method = "fdr"
)

gost_advanced <- gost(
  query             = genes_advanced,
  organism          = "hsapiens",
  sources           = c("GO:BP", "KEGG", "REAC"),
  significant       = TRUE,
  correction_method = "fdr"
)

gostplot(gost_early,    interactive = FALSE, capped = TRUE)
gostplot(gost_advanced, interactive = FALSE, capped = TRUE)

cat("\nTop pathways - early tumor genes :\n")
print(head(gost_early$result[,
                             c("term_name", "p_value", "term_size", "intersection_size")], 10))
# Result : cytoplasmic translation (p=1e-100), ribosome biogenesis
#          -> highly active and coordinated translational program
#          -> signature of well-differentiated early-stage cells

cat("\nTop pathways - advanced tumor genes :\n")
print(head(gost_advanced$result[,
                                c("term_name", "p_value", "term_size", "intersection_size")], 10))
# Result : anatomical structure morphogenesis, macromolecule
#          localization, cellular response to nitrogen compound
#          -> cell remodeling and migration signature
#          -> consistent with invasion in advanced tumors


# 11. Save all results 
saveRDS(results_stage, "data/processed/results_stage.rds")

write.csv(anno_pos,
          "data/processed/cpg_early_annotated.csv",
          row.names = FALSE)
write.csv(anno_neg,
          "data/processed/cpg_advanced_annotated.csv",
          row.names = FALSE)

weights_df <- data.frame(
  gene           = rownames(weights_rna),
  weight_Factor5 = weights_rna[, 1]
)
weights_df <- weights_df[order(weights_df$weight_Factor5,
                               decreasing = TRUE), ]
write.csv(weights_df,
          "data/processed/rna_weights_Factor5.csv",
          row.names = FALSE)

write.csv(
  gost_early$result[, c("term_name", "p_value",
                        "term_size", "intersection_size")],
  "data/processed/enrichment_early.csv",
  row.names = FALSE
)
write.csv(
  gost_advanced$result[, c("term_name", "p_value",
                           "term_size", "intersection_size")],
  "data/processed/enrichment_advanced.csv",
  row.names = FALSE
)

cat("\nAll results saved in data/processed/\n")

# BIOLOGICAL CONCLUSION
# Factor5 identified by MOFA+ captures an opposition between :
#   Early tumors  : active ribosomal translation program (RPL/RPS)
#                   methylation of FZD10 and XKR6
#                   -> functional well-differentiated cells
#   Advanced tumors : cell motility and remodeling (MYH9, MYH14)
#                     methylation of NF1, HOXD10, TJP2
#                     -> loss of tumor suppressors, dedifferentiation
# Validated by functional enrichment on GO:BP, KEGG and Reactome


