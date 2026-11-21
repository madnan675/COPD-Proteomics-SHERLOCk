# COPD-Proteomics-SHERLOCk

Analysis of clinical and Olink proteomics data from COPD and Control
cohorts, relating blood eosinophil status (Eosinophilic vs
Non-Eosinophilic, threshold 0.3 x10⁹/L) to baseline clinical
characteristics and proteomics from the bronchial epithelial lining fluid. The pipeline builds
baseline characteristics tables, runs differential protein expression,
visualizes and enriches the results, and identifies hub proteins from a
protein-protein interaction network.

## Overview

The analysis proceeds through six stages:

1. **Baseline characteristics** — clinical Table 1 for the Control and
   COPD cohorts, split by Eosinophilic/Non-Eosinophilic, with
   t-test/Mann-Whitney/Fisher's Exact as appropriate.
2. **Data merging** — clinical data merged with Olink proteomics (NPX)
   values by Subject_ID, patients classified as
   Eosinophilic/Non-Eosinophilic, and split into COPD-only/Control-only
   subsets.
3. **Differential protein expression** — `limma` moderated t-test
   comparing Eosinophilic vs Non-Eosinophilic, adjusting for age, sex,
   and smoking status. Primary comparison in COPD, secondary in
   Control.
4. **Visualisation & GSEA** — boxplots of top hits, a volcano plot, and
   Gene Set Enrichment Analysis (GO Biological Process, KEGG) on the
   full ranked protein list.
5. **Pathway enrichment** — over-representation analysis (GO
   Biological Process, KEGG) on the significant and disease-specific
   protein sets, with network visualisations.
6. **PPI hub protein analysis** — protein-protein interaction network
   built from STRING data, proteins ranked by degree
   and annotated with differential expression results.

## Repository structure

```
.
├── scripts/
│   ├── 01_Baseline_Control.R
│   ├── 02_Baseline_COPD.R
│   ├── 03_Merged_Clinical_Proteomics.R
│   ├── 04_Differential_Protein_Expression.R
│   ├── 05_Updated_DEP_GSEA.R
│   ├── 06_Pathway_Enrichment_COPD.R
│   └── 07_PPI_Hub_Proteins.R
├── data/
│   ├── input/        # source Excel/CSV/TSV files + merged dataset
│   └── results/       # all generated tables, figures, and Excel outputs
└── README.md
```

## Requirements

This is an R analysis (R ≥ 4.4 recommended). Install the required
packages with:

```r
# CRAN packages
install.packages(c(
  "tidyverse", "dplyr", "readxl", "writexl", "openxlsx", "huxtable",
  "ggplot2", "ggrepel"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "clusterProfiler", "org.Hs.eg.db", "fgsea", "enrichplot"))
```

## Data

The raw data are patient-level clinical and proteomics data and are
not included in this repository (see `.gitignore`). To run the
pipeline you will need access to:

- Control cohort clinical/CRF data (`control_Group.xlsx`)
- COPD cohort clinical/CRF data (`COPD_Group.xlsx`)
- Clinical data with COPD/Control classification (`clinical_data_with_classification.xlsx`)
- Olink proteomics NPX values (`proteomics_patients_x_2395proteins.csv`)
- STRING database interaction export (`string_interactions.tsv`)

## Setup

Clone the repository, then create the expected data folders at the
repo root:

```r
dir.create("data/input", recursive = TRUE)
dir.create("data/results", recursive = TRUE)
```

Place the raw input files listed above into `data/input/`. Open the
project in RStudio (or set your working directory to the repo root)
and run the scripts as described below — each script uses relative
paths (`./data/input/...`, `./data/results/...`), so they resolve
correctly as long as R is running from the repo root.

## Usage

Run the scripts in this order:

| Step | Script | Description |
|---|---|---|
| 1 | `scripts/01_Baseline_Control.R` | Baseline characteristics table for the Control cohort |
| 2 | `scripts/02_Baseline_COPD.R` | Baseline characteristics table for the COPD cohort |
| 3 | `scripts/03_Merged_Clinical_Proteomics.R` | Merges clinical and proteomics data, classifies Eosinophilic status |
| 4 | `scripts/04_Differential_Protein_Expression.R` | `limma` differential protein expression (COPD primary, Control secondary) |
| 5 | `scripts/05_Updated_DEP_GSEA.R` | Boxplots, volcano plot, and GSEA (GO Biological Process, KEGG) |
| 6 | `scripts/06_Pathway_Enrichment_COPD.R` | GO/KEGG over-representation analysis and network plots |
| 7 | `scripts/07_PPI_Hub_Proteins.R` | STRING PPI network hub protein analysis |

Step 3 depends on the raw clinical and proteomics files only (not on
Steps 1–2). Steps 4–7 each depend on the `.xlsx` outputs of the step(s)
before them in `data/`.

## Output

Key outputs (all saved to `data/results/`, except the merged dataset
which is saved to `data/input/`) include:

- Baseline characteristics tables for Control and COPD (`.xlsx`)
- Merged clinical + proteomics dataset (`.xlsx`)
- Full and significant differential protein expression results
  (`.xlsx`) — COPD significance: FDR < 0.20 & |log2FC| > 0.5; Control
  did not reach FDR < 0.20, so its top 200 proteins by raw p-value are
  reported instead as an exploratory reference set, not as
  FDR-significant hits
- Boxplots and volcano plot of top differentially expressed proteins
  (`.pdf`/`.png`)
- GSEA and pathway enrichment results and dotplots (`.xlsx`, `.pdf`/`.png`)
- Protein-pathway network plots (`.pdf`/`.png`)
- Hub protein table from the PPI network, annotated with differential
  expression results (`.xlsx`)

## Authors

Muhammad Adnan, GRIAC Bioinformatics

## Citation

(citation to be added on publication).

