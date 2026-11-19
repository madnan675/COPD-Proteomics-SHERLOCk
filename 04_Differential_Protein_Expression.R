
#DIFFERENTIAL PROTEIN EXPRESSION ANALYSIS

# PURPOSE:
# Tests which proteins differ between Eosinophilic and Non-Eosinophilic
# subgroups, adjusted for age, sex, and smoking status, using limma's
# moderated t-test (empirical Bayes).
#   Primary:   COPD Eosinophilic vs COPD Non-Eosinophilic
#   Secondary: Control Eosinophilic vs Control Non-Eosinophilic
#
# WHAT THE SCRIPT DOES:
#   1. Loads merged clinical + proteomics data (full cohort, COPD-only,
#      and Control-only sheets) and confirms group sizes.
#   2. Identifies the protein (NPX) columns and converts them to numeric.
#   3. Builds protein-by-patient expression matrices for COPD and Control
#      (limma requires proteins as rows, patients as columns).
#   4. Sets Non_Eosinophilic as the reference level, so a positive logFC
#      means a protein is higher in the Eosinophilic group.
#   5. Builds design matrices adjusting for age, sex, and smoking status.
#   6. Runs limma (lmFit -> eBayes -> topTable) separately for COPD and
#      Control, producing full ranked results for every protein.
#   7. Flags significant proteins:
#        - COPD:    FDR (adj.P.Val) < 0.20 (no protein reached FDR < 0.05)
#        - Control: no protein reached FDR < 0.20, so instead the top 200
#                    proteins by raw p-value are reported as exploratory
#                    hits (NOT FDR-significant)
#   8. Cross-compares COPD and Control results to find:
#        - Shared_Both:   proteins significant in COPD that also appear
#                          in Control's top 200 (with full stats from
#                          both groups side by side)
#        - COPD_Specific: proteins significant in COPD only — the most
#                          disease-relevant hits (with full stats from
#                          both groups side by side)
#   9. Exports everything to one Excel file with a sheet per result set.
#
# OUTPUT:
#   ./data/results/DEP_Results_Stage3.xlsx



#STEP 1 — INSTALL AND LOAD REQUIRED PACKAGES


install.packages("tidyverse")    
install.packages("readxl")       
install.packages("writexl")      

#limma comes from Bioconductor needs its own installer
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")

# Load libraries into current session
library(tidyverse)   
library(readxl)      
library(writexl)    
library(limma)      



#STEP 2 — LOAD MERGED Clinical and proteomics EXCEL FILE

merged_data  <- read_excel("data/input/Merged_Clinical_Proteomics.xlsx",
                            sheet = "Full_Merged")

copd_data    <- read_excel("data/input/Merged_Clinical_Proteomics.xlsx",
                            sheet = "COPD_only")

control_data <- read_excel("data/input/Merged_Clinical_Proteomics.xlsx",
                            sheet = "Control_only")

#Verify all three sheets loaded correctly ---
dim(merged_data)     
dim(copd_data)        
dim(control_data)     

#Confirm group column and eosinophilic subgroup column exist
colnames(merged_data)[1:15]   #shows first 15 column names

#Check subgroup sizes
table(copd_data$Eos_Group)
table(control_data$Eos_Group)


# STEP 3 — IDENTIFY PROTEIN COLUMNS
colnames(merged_data)

#9 clinical columns are:
# 1=Subject_ID, 2=classification, 3=eosinophils, 4=age, 5=sex,
# 6=smoking_status, 7=Eos_Group
#proteins start at column 8 

protein_names <- colnames(merged_data)[8:ncol(merged_data)]

#CHECK ---
length(protein_names)       
head(protein_names, 10)     #preview first 10 protein names
tail(protein_names, 10)     #preview last 10 protein names



#STEP 4 — CONVERT PROTEIN COLUMNS TO NUMERIC

# After reading from Excel, protein NPX values can import as text.
# mutate(across()) applies a function across multiple columns at once.
# as.numeric() converts each value to a number.
# -Subject_ID means: do this to all columns EXCEPT Subject_ID.

merged_data <- merged_data %>%
  mutate(across(all_of(protein_names), as.numeric))

copd_data <- copd_data %>%
  mutate(across(all_of(protein_names), as.numeric))

control_data <- control_data %>%
  mutate(across(all_of(protein_names), as.numeric))



# STEP 5 — BUILD PROTEIN EXPRESSION MATRICES FOR LIMMA

#limma requires a specific matrix format:
#   Rows    = proteins  (one row per protein)
#   Columns = patients  (one column per patient)
#my dataframe is the opposite (rows=patients, cols=proteins)
#so extract the protein columns then TRANSPOSE (flip rows/cols).

# ── COPD expression matrix ───────────────────────────────────────────
# select() pulls only the protein columns from the COPD dataframe
copd_matrix <- copd_data %>%
  select(all_of(protein_names)) %>%
  as.matrix()

# Assign patient IDs as row names  needed for t() transposition
rownames(copd_matrix) <- copd_data$Subject_ID

# t() transposes: now rows=proteins, columns=patients
copd_expr <- t(copd_matrix)

# ── Control expression matrix ────────────────────────────────────────
control_matrix <- control_data %>%
  select(all_of(protein_names)) %>%
  as.matrix()

rownames(control_matrix) <- control_data$Subject_ID
control_expr <- t(control_matrix)

# CHECK ---
dim(copd_expr)      #proteins x COPD patients
dim(control_expr)   #proteins x Control patients


#STEP 6 — DEFINE EOSINOPHILIC GROUP FACTORS

# limma needs the comparison group as a factor (categorical variable).
# levels = sets the ORDER of groups.
# The FIRST level is the REFERENCE group (baseline for comparison).
# Setting Non_Eosinophilic first means:
#   positive logFC = protein is HIGHER in Eosinophilic group
#   negative logFC = protein is LOWER in Eosinophilic group

# ── COPD group factor ────────────────────────────────────────────────
copd_group <- factor(copd_data$Eos_Group,
                     levels = c("Non_Eosinophilic", "Eosinophilic"))

# ── Control group factor ─────────────────────────────────────────────
control_group <- factor(control_data$Eos_Group,
                        levels = c("Non_Eosinophilic", "Eosinophilic"))

#CHECK COPD & control groups Eosinophilic vs Non-Eosinophilic
levels(copd_group)     
levels(control_group)   
table(copd_group)     
table(control_group)



#STEP 7 — BUILD DESIGN MATRICES WITH COVARIATES

# The design matrix is the mathematical description of model.
# ~ means "modelled by"
# The first term (copd_group) is main effect — what you are testing
# Each + term is a covariate — what you are controlling for
# model.matrix() converts this formula into a numeric matrix for limma
#
# Covariates included:
#   age           — biological confounder, always adjust for
#   sex           — biological confounder, always adjust for
#   smoking_status — major COPD confounder


#PRIMARY: COPD design matrix ──────────────────────────────────────
design_copd <- model.matrix(~ copd_group +
                               copd_data$age +
                               copd_data$sex +
                               copd_data$smoking_status)
                               

#SECONDARY: Control design matrix ─────────────────────────────────
design_control <- model.matrix(~ control_group +
                                  control_data$age +
                                  control_data$sex +
                                  control_data$smoking_status)
                              

#CHECK colunm names of design matrices to ensure covariates are included
colnames(design_copd)
colnames(design_control)


#STEP 8 — RUN LIMMA: PRIMARY ANALYSIS (COPD)

#limma works in three sequential steps:
#
# Step A: lmFit()
#   Fits a linear model for EVERY protein simultaneously.
#   This is equivalent to running a separate linear regression
#   for each protein, adjusted for covariates.
#
# Step B: eBayes()
#   Applies Empirical Bayes moderation.
#   This borrows variance information ACROSS all proteins to
#   stabilise estimates for individual proteins.
#   It shrinks extreme variance values toward the overall mean.
#   This is what makes limma far more reliable than running
#   individual t-tests, especially with small sample sizes.
#
# Step C: topTable()
#   Extracts the ranked results table.
#   coef=2 refers to the second column of the design matrix
#   which is "copd_groupEosinophilic" — main comparison.
#   number=Inf returns ALL proteins (not just top 10).
#   adjust.method="BH" applies Benjamini-Hochberg FDR correction.

# Step A: Fit linear model
fit_copd <- lmFit(copd_expr, design_copd)

# Step B: Apply empirical Bayes moderation
fit_copd <- eBayes(fit_copd)

# Step C: Extract full results table
results_copd <- topTable(fit_copd,
                         coef          = 2,
                         number        = Inf,
                         adjust.method = "BH",
                         sort.by       = "P")

#Add protein name as a proper column (currently stored as row names)
results_copd$Protein <- rownames(results_copd)

#Preview top 20 results ---
head(results_copd, 20)

#Understanding each column in the output:
#
# Protein    → name of the protein
# logFC      → log2 fold change
#              positive = higher in Eosinophilic group
#              negative = lower in Eosinophilic group
#              logFC of 1.0 means 2-fold higher in Eosinophilic
#              logFC of -1.0 means 2-fold lower in Eosinophilic
# AveExpr    → average NPX expression across all patients in this group
# t          → moderated t-statistic
#              larger absolute value = stronger effect relative to noise
# P.Value    → unadjusted (raw) p-value — do NOT use for significance
# adj.P.Val  → FDR-corrected p-value (Benjamini-Hochberg)
#              USE THIS for deciding significance
#              threshold: adj.P.Val < 0.05
# B          → log-odds that the protein is truly differentially expressed
#              higher B = more confident the protein is a real hit



#STEP 9 — RUN LIMMA: SECONDARY ANALYSIS (CONTROL)

#Exact same three steps — using control data

fit_control <- lmFit(control_expr, design_control)
fit_control <- eBayes(fit_control)

results_control <- topTable(fit_control,
                            coef          = 2,
                            number        = Inf,
                            adjust.method = "BH",
                            sort.by       = "P")

results_control$Protein <- rownames(results_control)

#preview top 20 results
head(results_control, 20)



#STEP 10 — IDENTIFY SIGNIFICANT PROTEINS

# Standard significance thresholds:
#   FDR (adj.P.Val) < 0.05   — statistical significance after correction

#
# If NO proteins reach FDR < 0.05 (common with small samples):
# Relax to FDR < 0.20 — acceptable in exploratory studies


library(dplyr) 

#PRIMARY: COPD significant proteins ──────────────────────────────
sig_copd <- results_copd %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 0.5) %>%
  arrange(adj.P.Val)

#use Fdr< 0.20 because no proteins reached at fdr 0.05
sig_copd <- results_copd %>%
  filter(adj.P.Val < 0.20) %>%
  arrange(adj.P.Val)



# Split by direction of change
sig_copd_up   <- sig_copd %>% filter(logFC > 0)   # higher in Eosinophilic
sig_copd_down <- sig_copd %>% filter(logFC < 0)   # lower in Eosinophilic


cat("=== PRIMARY ANALYSIS: COPD ===\n")
cat("Total significant proteins:         ", nrow(sig_copd), "\n")
cat("Upregulated in COPD Eosinophilic:   ", nrow(sig_copd_up), "\n")
cat("Downregulated in COPD Eosinophilic: ", nrow(sig_copd_down), "\n\n")

#SECONDARY: Control significant proteins ──────────────────────────
sig_control <- results_control %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 0.5) %>%
  arrange(adj.P.Val)

#SELECT TOP 200 CONTROL PROTEINS BY P-VALUE REGARDLESS OF FDR SIGNIFICANCE
top200_control <- results_control %>%
  arrange(P.Value) %>%         # rank by raw p-value smallest first
  head(200) %>%                # take top 200
  dplyr::select(Protein, logFC, P.Value, adj.P.Val)


#Proteins in BOTH analyses ─────────────────────────────────────────
# sig_control is now the top 200 by p-value (not FDR-significant), so
# this identifies COPD-significant proteins that are ALSO among the
# top 200 in Control by nominal p-value — exploratory overlap


shared_proteins <- intersect(sig_copd$Protein, top200_control$Protein)
cat("COPD-significant proteins also in Control's top 200:", length(shared_proteins), "\n")
if(length(shared_proteins) > 0) print(shared_proteins)

shared_stats <- results_copd %>%
  filter(Protein %in% shared_proteins) %>%
  select(Protein, logFC, AveExpr, t, P.Value, adj.P.Val, B) %>%
  rename_with(~ paste0(., "_COPD"), -Protein) %>%
  left_join(
    results_control %>%
      filter(Protein %in% shared_proteins) %>%
      select(Protein, logFC, AveExpr, t, P.Value, adj.P.Val, B) %>%
      rename_with(~ paste0(., "_Control"), -Protein),
    by = "Protein"
  ) %>%
  arrange(adj.P.Val_COPD)



#COPD-specific proteins ───────────────────────────────────────────
# Significant in COPD but NOT in Control
# These are the most disease-relevant hits for your primary question

copd_specific <- setdiff(sig_copd$Protein, top200_control$Protein)
cat("Proteins significant in COPD only:", length(copd_specific), "\n")
if(length(copd_specific) > 0) print(copd_specific)



#STEP 11 — SAVE STAGE 3 RESULTS TO EXCEL

library(writexl)
write_xlsx(
  list(
    #PRIMARY results
    "COPD_All_Proteins"    = results_copd,
    "COPD_Significant"     = sig_copd,
    "COPD_Up"              = sig_copd_up,
    "COPD_Down"            = sig_copd_down,

    #SECONDARY results
    "top200 Control" = top200_control,
   

    #Cross-comparison
    "Shared_Both"          = shared_stats,
    "COPD_Specific"        = data.frame(Protein = copd_specific)
  ),
  path = "./data/results/DEP_Results_Stage3.xlsx"
)


