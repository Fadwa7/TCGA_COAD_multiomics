# Multiomics integration of colorectal cancer - TCGA-COAD

**Author** : Fadwa EL KHADDAR  
**Tools** : R, MOFA+, curatedTCGAData, gprofiler2  
**Data** : TCGA-COAD - 271 patients, RNA-seq + DNA methylation (450k)

---
**[View the full interactive report](https://fadwa7.github.io/TCGA_COAD_multiomics/)**

## What is multiomics ?

Biology does not happen at a single level. A cell is simultaneously shaped by its genome, its epigenome, its transcriptome, its proteome, and its metabolome. Each of these layers captures a different aspect of cellular state - and none of them alone tells the full story.

**Multiomics** is the simultaneous integration of multiple molecular measurement types from the same biological system. Instead of asking "which genes are differentially expressed ?", multiomics asks "which molecular patterns, shared across DNA methylation, gene expression, and other layers, distinguish biological states ?"

The core challenge is not data volume - it is heterogeneity. Each omics layer has its own scale, its own noise structure, its own sparsity. Naive concatenation of matrices destroys the biological signal. Proper integration requires methods that can find shared latent structure across fundamentally different data types.

---

## What is MOFA+ ?

**MOFA+** (Multi-Omics Factor Analysis, version 2) is a statistical framework for unsupervised integration of multiple omics datasets. It was published in *Genome Biology* (Argelaguet et al., 2020) and has become a reference method in the field.

### The core idea

MOFA+ learns a set of **latent factors** - compressed representations that capture the main sources of variation across all omics layers simultaneously. Each factor represents a biological signal that is expressed, to varying degrees, across modalities.

```
RNA-seq (3000 genes)       \
                            --> MOFA+ --> 15 latent factors
Methylation (3000 CpG)     /
```

Each factor has two components :

- **Scores** - one value per patient, indicating how strongly each patient expresses that factor
- **Weights** - one value per feature (gene or CpG for example ), indicating how much each feature contributes to that factor

### Why MOFA+ and not other methods ?

| Method | Limitation |
|--------|-----------|
| PCA | Single modality only |
| Concatenation + PCA | Does not account for modality-specific noise |
| DIABLO | Requires supervised labels |
| iCluster | Does not scale well, complex to interpret |
| **MOFA+** | Unsupervised, multi-group, interpretable weights, handles missing data |

MOFA+ is particularly suited for exploratory analysis where you do not know in advance which biological signal drives the heterogeneity - which is exactly the case in cancer multiomics.

---

## Biological question

> In colorectal cancer (TCGA-COAD cohort), which shared molecular patterns between DNA methylation and gene expression distinguish early-stage tumors (stage I-II) from advanced tumors (stage III-IV) ?

### Why this question ?

Colorectal cancer is molecularly heterogeneous , two patients with the same clinical stage can have very different molecular profiles and outcomes. Understanding what drives this heterogeneity at the multiomics level is essential for identifying robust biomarkers and understanding tumor progression mechanisms.

We initially aimed to predict overall survival, but with only 74 death events in 271 patients, statistical power was insufficient for individual factor detection. Tumor stage (early vs advanced) provided a more robust and biologically grounded validation endpoint.

---

## Data

| Source | Description |
|--------|-------------|
| TCGA-COAD | Colon adenocarcinoma cohort |
| RNA-seq | HiSeq IlluminaGA, RSEM-normalized, 326 samples |
| DNA methylation | Illumina 450k array, beta values, 333 samples |
| Clinical | Vital status, follow-up time, tumor stage |
| Final cohort | **271 patients** with complete data across all modalities |

Data accessed via the `curatedTCGAData` Bioconductor package (version 2.0.1). Only primary solid tumors (TCGA tissue code `01`) were retained to ensure biological homogeneity.

---

## Methods summary

```
1. Data download      curatedTCGAData - RNA-seq + methylation 450k
        |
2. Clinical table     Overall survival + tumor stage (early / advanced)
        |
3. QC and filtering   Primary tumors only (code 01)
                      Top 3000 most variable features per modality
                      CpG with > 20% missing values removed
        |
4. Patient alignment  271 patients common to all 3 sources
        |
5. MOFA+ training     15 latent factors, seed = 42
        |
6. Survival analysis  Univariate Cox regression - 15 factors
        |
7. Stage validation   Wilcoxon test - early vs advanced tumors
        |
8. Interpretation     RNA weights, CpG annotation (Illumina 450k)
        |
9. Enrichment         gprofiler2 - GO:BP, KEGG, Reactome (FDR correction)
```

---

## Key results

### Variance explained by MOFA+

- **Methylation** : ~47% of total variance captured -highly structured signal consistent with the known epigenetic heterogeneity of colorectal cancer
- **RNA-seq** : ~13% of total variance captured — expected for transcriptomics data with 3000 genes

### Factor5 - the most biologically informative factor

Factor5 significantly distinguishes early from advanced tumors (Wilcoxon p = 0.003) and captures a clear biological opposition :

**Early tumors (high Factor5 score)**
- RNA : overexpression of ribosomal genes (RPL/RPS family) : active and coordinated translational machinery, signature of well-differentiated cells
- Methylation : hypermethylation of FZD10 (Wnt receptor), XKR6 (apoptosis regulator), ZSCAN18 (transcription factor)
- Enriched pathways : cytoplasmic translation (p = 10⁻¹⁰⁰), ribosome biogenesis, gene expression

**Advanced tumors (low Factor5 score)**
- RNA : overexpression of MYH9/MYH14 (cell motility), PRKDC (DNA repair), RNF43 (Wnt pathway), HNF4A (differentiation)
- Methylation : hypermethylation of NF1 (tumor suppressor), HOXD10 (differentiation), TJP2 (tight junctions), BAG3 (apoptosis)
- Enriched pathways : anatomical structure morphogenesis, macromolecule localization, cellular remodeling

> **Conclusion** : Factor5 captures a transition from an active ribosomal translation program in early tumors to a cell motility and dedifferentiation program in advanced tumors  consistent with known mechanisms of colorectal cancer progression.

---

## Repository structure

```
TCGA_COAD_multiomics/
├── README.md
├── scripts/
│   ├── 01_data_preparation.R     # Download, QC, alignment, dimension reduction
│   └── 02_mofa_analysis.R        # MOFA+ training, survival, stage, interpretation
├── data/
│   └── processed/                # Generated by scripts (not tracked by git)
│       ├── rna_top3000.rds
│       ├── meth_top3000.rds
│       ├── clinical_surv.rds
│       ├── mofa_model.rds
│       ├── results_stage.rds
│       ├── rna_weights_Factor5.csv
│       ├── cpg_early_annotated.csv
│       ├── cpg_advanced_annotated.csv
│       ├── enrichment_early.csv
│       └── enrichment_advanced.csv
└── report/
    └── report.Rmd                # Full reproducible report with figures
```

---

## How to reproduce

```r
# 1. Install required packages
BiocManager::install(c(
  "curatedTCGAData", "MultiAssayExperiment", "TCGAutils",
  "DESeq2", "minfi", "matrixStats", "MOFA2",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19"
))
install.packages(c("tidyverse", "survival", "survminer",
                   "ggplot2", "gridExtra", "gprofiler2"))

# 2. Run scripts in order
source("scripts/01_data_preparation.R")
source("scripts/02_mofa_analysis.R")
```

> Note : MOFA+ requires Python via `basilisk`. The first run will automatically install a conda environment - this may take 5-10 minutes.

---

## References

- Argelaguet R. et al. (2020). MOFA+ : a statistical framework for comprehensive integration of multi-modal single-cell data. *Genome Biology*, 21, 111.
- Colaprico A. et al. (2016). TCGAbiolinks : an R/Bioconductor package for integrative analysis with GDC data. *Nucleic Acids Research*, 44(8).
- Rainer J. et al. (2022). curatedTCGAData : Curated Data From The Cancer Genome Atlas as MultiAssayExperiment Objects. *Bioconductor*.
- Kolberg L. et al. (2023). g:Profiler — interoperable web service for functional enrichment analysis and gene identifier mapping. *Nucleic Acids Research*, 51(W1).
