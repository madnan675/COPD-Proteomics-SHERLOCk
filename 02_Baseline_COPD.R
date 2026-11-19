
#BASELINE CHARACTERISTICS TABLE — COPD GROUP (Eosinophilic vs Non-Eosinophilic)

#PURPOSE:
  # Builds a "Table 1"-style baseline characteristics table comparing the COPD
  # group split into two subgroups by blood eosinophil count:
  #   - Eosinophilic:      blood eosinophils >= 0.3 x10^9/L
  #   - Non-Eosinophilic:  blood eosinophils <  0.3 x10^9/L
  #
  # WHAT THE SCRIPT DOES:
  #1. Loads the COPD group data, removes a known duplicate column, and
  #      calculates BMI from weight/height.
  #2. Renames raw dataset columns to short, consistent internal names.
  #3. Removes subjects with missing blood eosinophil values (required for
  #      group classification).
  #4. Classifies each subject as Eosinophilic / Non-Eosinophilic and flags
  #      subjects with >0 exacerbations vs 0.
  #5. Runs Shapiro-Wilk tests to check normality of continuous variables.
  #6. Summarizes each variable by group:
  #        - Normally distributed continuous variables -> Mean +/- SD, t-test
  #        - Non-normal continuous variables            -> Median [IQR], Mann-Whitney U
  #        - Categorical variables                       -> n (%), Fisher's Exact test
  #      Categorical variables (sex, ICS use, oral steroids, smoking status,
  #      exacerbations) get a header row (n/%, p-value) plus indented sub-rows
  #      for each category level.
  #7. Runs diagnostic checks (row presence, observation counts) to confirm
  #      every variable was pulled correctly before assembling the table.
  #8. Assembles all rows into one final table, organized into sections:
  #      demographics, post-bronchodilator spirometry, diffusion capacity,
  #      blood cell counts, HRCT, and impulse oscillometry — with a subject
  #      count (n) row at the top.
  #9. Exports the final table to a formatted Excel file and
  #      significant p-values (<0.05) highlighted in red.
  #
  # OUTPUT:
  #   ./data/results/Table1_COPDFinal_baseline.xlsx

#INSTALL AND LOAD PACKAGES


packages <- c("readxl", "dplyr", "openxlsx")

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
invisible(lapply(packages, install_if_missing))

library(readxl)
library(dplyr)
library(openxlsx)


# STEP 1 — LOAD COPD Clinical DATA AND CHECK COLUMN NAMES
copd_data <- read_excel(
"data/input/COPD_Group.xlsx"
)
cat("Rows:", nrow(copd_data), "| Columns:", ncol(copd_data), "\n")
cat("\nColumn names in your file:\n")
print(names(copd_data))

# STEP 2 — RENAME COLUMNS TO INTERNAL SHORT NAMES

#remove duplicate column 
copd_data <- copd_data[,-617]

#calculate the BMI
copd_data$BMI <- copd_data$crf_kg / (copd_data$height^2)

copd_data <- copd_data %>%
  rename(
    # Demographics
    age             = "crf_age",
    sex             = "crf_gender",
    bmi             = "BMI",

    # Clinical variables
    pack_years      = "crf_packyears",
    exacerbations   = "crf_number_exa",
    oral_steroids   = "oral_steroids",
    ics_use         = "ics_use",
    crf_smoking     = "crf_smoking",
    cat_score       = "CATtotalscore_def",

    # Post-bronchodilator spirometry
    fev1_pct        = "post_FEV1_pctpred",
    fev1_fvc        = "post_FEV1FVC_pctpred",
    rv_tlc          = "post_RVTLC_pctpred",

    # Diffusion capacity
    kco             = "post_KCO_pctpred",
    tlco            = "post_TLCO_pctpred",

    # Blood cell counts
    eosinophils     = "blood_eusinophils",
    neutrophils     = "blood_neutrophils",
    leukocytes      = "blood_leukocyten",
    lymphocytes     = "blood_lymphocytes",
    monocytes       = "blood_monocytes",
    crp             = "blood_crp",

    # HRCT
    air_trapping    = "Lungs_airtrapping_%",
    emphysema       = "Lungs_emphysema_%",
    mucus_plugs     = "Lungs_mucus_plugs_count",
    mucus_volume    = "Lungs_mucus_vol_ml",
    pi10            = "Lungs_Pi10",

    # Impulse oscillometry
    ax              = "average_ios_ax",
    x5              = "average_ios_x5",
  )

cat("\nRenamed columns confirmed:\n")
print(names(copd_data))



#STEP3 CLEAN MISSING VALUES in eosinophils variables

# Remove rows where eosinophils is missing
copd_data <- copd_data %>%
  filter(!is.na(eosinophils))

cat("\nRows remaining after removing missing eosinophils:",
    nrow(copd_data), "\n")

# STEP 4 — CLASSIFY EOSINOPHILIC / NON-EOSINOPHILIC
#threshold. table() at the end shows how many patients fall into each group.

copd_data <- copd_data %>%
  mutate(eos_class = if_else(eosinophils >= 0.3,
                              "COPD-Eosinophilic",
                              "COPD-Non_Eosinophilic"))

cat("\nGroup sizes:\n")
print(table(copd_data$eos_class))


# STEP 5 — SET GROUP ORDER AS FACTOR


GROUP_ORDER <- c("COPD-Eosinophilic", "COPD-Non_Eosinophilic")

copd_data <- copd_data %>%
  mutate(eos_class = factor(eos_class, levels = GROUP_ORDER))


# STEP 6 — CREATE EXACERBATIONS FLAG (>0 vs 0)

copd_data <- copd_data %>%
  mutate(exacerbations_flag = if_else(exacerbations > 0, ">0", "0"))

cat("\nExacerbations flag distribution:\n")
print(table(copd_data$exacerbations_flag))


# STEP 7 — DEFINE VARIABLE LISTS AND DISPLAY LABELS

# Normally distributed — reported as Mean +/- SD
CONTINUOUS_MEAN_SD <- c("age", "bmi", "fev1_pct", "fev1_fvc",
                         "rv_tlc", "tlco", "kco", "pi10", "cat_score")

# Non-normally distributed — reported as Median [IQR]
CONTINUOUS_MEDIAN_IQR <- c("pack_years", "eosinophils", "neutrophils",
                            "leukocytes", "lymphocytes", "monocytes", "crp",
                            "air_trapping", "emphysema", "mucus_plugs",
                            "mucus_volume", "ax", "x5")

# Categorical variables
CATEGORICAL_VARS <- c("sex", "ics_use", "oral_steroids",
                      "crf_smoking", "exacerbations_flag")

# Display labels — exactly as they appear in the final table
LABELS <- c(
  age                = "Age, years",
  bmi                = "BMI, kg/m2",
  pack_years         = "Pack-years",
  cat_score          = "CAT total score",
  fev1_pct           = "FEV1% pred",
  fev1_fvc           = "FEV1/FVC%",
  rv_tlc             = "RV/TLC %pred",
  tlco               = "TLCO% pred",
  kco                = "KCO% pred",
  pi10               = "Pi10, mm",
  eosinophils        = "Eosinophils x10^9/L",
  neutrophils        = "Neutrophils x10^9/L",
  leukocytes         = "Leukocytes x10^9/L",
  lymphocytes        = "Lymphocytes x10^9/L",
  monocytes          = "Monocytes x10^9/L",
  crp                = "CRP, mg/L",
  air_trapping       = "Air trapping, %",
  emphysema          = "Emphysema, % (LAA < -950 HU)",
  mucus_plugs        = "Mucus plugs count",
  mucus_volume       = "Mucus volume, mL",
  ax                 = "AX, kPa.s/L",
  x5                 = "X5, kPa.s/L",
  sex                = "Sex, n (%)",
  ics_use            = "ICS use, n (%)",
  oral_steroids      = "Oral steroids, n (%)",
  crf_smoking        = "Smoking status, n (%)",
  exacerbations_flag = "Number of exacerbations, n (%)"
)

# STEP 8 — NORMALITY CHECK (SHAPIRO-WILK TEST)
# What this does: tests each continuous variable for normality.
# NORMAL (p > 0.05)     → keep in CONTINUOUS_MEAN_SD (Mean +/- SD)
# NON-NORMAL (p < 0.05) → move to CONTINUOUS_MEDIAN_IQR (Median [IQR])
# After updating lists in Step 7, re-run from Step 7 downward.


all_continuous <- c(CONTINUOUS_MEAN_SD, CONTINUOUS_MEDIAN_IQR)

cat("\nShapiro-Wilk Normality Test Results:\n")
cat(strrep("-", 65), "\n")

for (var in all_continuous) {
  if (!var %in% names(copd_data)) next
  x <- na.omit(copd_data[[var]])
  if (length(x) < 3) {
    cat(sprintf("%-25s : too few values\n", var)); next
  }
  result <- tryCatch(shapiro.test(x), error = function(e) NULL)
  if (!is.null(result)) {
    decision <- if (result$p.value > 0.05) {
      "NORMAL     --> Mean +/- SD ok"
    } else {
      "NON-NORMAL --> use Median [IQR]"
    }
    cat(sprintf("%-25s : p = %.4f  -->  %s\n",
                var, result$p.value, decision))
  }
}


# STEP 9 — DEFINE HELPER FUNCTIONS

# 9a — Format p-value for display
# Converts raw decimal p-values to clean text: <0.001 or 3 decimal places
fmt_p <- function(p) {
  if (is.na(p)) return("—")
  if (p < 0.001) return("<0.001")
  formatC(p, digits = 3, format = "f")
}

# 9b — Mean +/- SD with independent t-test
# Used for normally distributed continuous variables
summarise_mean_sd <- function(data, var) {
  g1 <- na.omit(data[[var]][data$eos_class == GROUP_ORDER[1]])
  g2 <- na.omit(data[[var]][data$eos_class == GROUP_ORDER[2]])

  cell1 <- sprintf("%.2f \u00B1 %.2f", mean(g1), sd(g1))
  cell2 <- sprintf("%.2f \u00B1 %.2f", mean(g2), sd(g2))

  p <- if (length(g1) >= 2 && length(g2) >= 2) {
    tryCatch(t.test(g1, g2)$p.value, error = function(e) NA)
  } else NA

  c(Group1    = cell1,
    Group2    = cell2,
    `p-value` = fmt_p(p),
    Test      = "t-test",
    Stat      = "Mean \u00B1 SD")
}

# 9c — Median [IQR] with Mann-Whitney U test
# Used for non-normally distributed continuous variables
summarise_median_iqr <- function(data, var) {
  g1 <- na.omit(data[[var]][data$eos_class == GROUP_ORDER[1]])
  g2 <- na.omit(data[[var]][data$eos_class == GROUP_ORDER[2]])

  cell1 <- sprintf("%.2f [%.2f\u2013%.2f]",
                   median(g1), quantile(g1, 0.25), quantile(g1, 0.75))
  cell2 <- sprintf("%.2f [%.2f\u2013%.2f]",
                   median(g2), quantile(g2, 0.25), quantile(g2, 0.75))

  p <- if (length(g1) >= 2 && length(g2) >= 2) {
    tryCatch(wilcox.test(g1, g2)$p.value, error = function(e) NA)
  } else NA

  c(Group1    = cell1,
    Group2    = cell2,
    `p-value` = fmt_p(p),
    Test      = "Mann-Whitney U",
    Stat      = "Median [IQR]")
}

# 9d — Detect positive category (handles Male/Female, Yes/No, 1/0)
positive_level <- function(x) {
  vals <- unique(na.omit(as.character(x)))
  pos  <- c("1", "yes", "y", "true", "male", "m")
  hit  <- vals[tolower(vals) %in% pos]
  if (length(hit) > 0) return(hit[1])
  sort(vals)[length(vals)]
}

# 9e — Categorical header row with overall count/% and Fisher p-value
# Produces the main variable row showing the primary category count
make_header_row <- function(data, var, label, primary_level) {
  g1 <- na.omit(data[[var]][data$eos_class == GROUP_ORDER[1]])
  g2 <- na.omit(data[[var]][data$eos_class == GROUP_ORDER[2]])

  n1 <- sum(as.character(g1) == primary_level); t1 <- length(g1)
  n2 <- sum(as.character(g2) == primary_level); t2 <- length(g2)

  cell1 <- sprintf("%d (%.1f%%)", n1, if (t1 > 0) 100*n1/t1 else 0)
  cell2 <- sprintf("%d (%.1f%%)", n2, if (t2 > 0) 100*n2/t2 else 0)

  ct <- matrix(c(n1, t1 - n1, n2, t2 - n2), nrow = 2, byrow = TRUE)
  p  <- tryCatch(fisher.test(ct)$p.value, error = function(e) NA)

  data.frame(
    Variable         = label,
    Group1           = cell1,
    Group2           = cell2,
    `p-value`        = fmt_p(p),
    Test             = "Fisher's Exact",
    Stat             = "n (%)",
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )
}

# 9f — Categorical sub-rows (indented with two spaces, no p-value)
# Produces one row per category level below the header row
summarise_categorical_both <- function(data, var, level_labels) {
  rows <- list()
  for (lev in level_labels) {
    g1 <- na.omit(data[[var]][data$eos_class == GROUP_ORDER[1]])
    g2 <- na.omit(data[[var]][data$eos_class == GROUP_ORDER[2]])

    n1 <- sum(as.character(g1) == lev); t1 <- length(g1)
    n2 <- sum(as.character(g2) == lev); t2 <- length(g2)

    cell1 <- sprintf("%d (%.1f%%)", n1, if (t1 > 0) 100*n1/t1 else 0)
    cell2 <- sprintf("%d (%.1f%%)", n2, if (t2 > 0) 100*n2/t2 else 0)

    rows[[lev]] <- data.frame(
      Variable         = paste0("  ", lev),
      Group1           = cell1,
      Group2           = cell2,
      `p-value`        = "",
      Test             = "",
      Stat             = "n (%)",
      check.names      = FALSE,
      stringsAsFactors = FALSE
    )
  }
  bind_rows(rows)
}

# 9g — Fix Group1/Group2 column names to match results_df
# Prevents the duplicate column / NA issue when binding categorical blocks
fix_group_names <- function(block, real_names) {
  names(block)[names(block) == "Group1"] <- real_names[1]
  names(block)[names(block) == "Group2"] <- real_names[2]
  block
}


# STEP 10 — BUILD CONTINUOUS VARIABLE ROWS
# What this does: loops through both continuous variable lists, calls the
# correct summary function for each variable, and stacks all results into
# results_df. The if() check skips any column not found in copd_data
# instead of crashing. Uses LABELS[] key lookup to avoid string mismatches.


results_list <- list()

# Mean +/- SD variables
for (var in CONTINUOUS_MEAN_SD) {
  if (!var %in% names(copd_data)) {
    cat("Skipping (not found in data):", var, "\n"); next
  }
  row_df              <- as.data.frame(
                           t(summarise_mean_sd(copd_data, var)),
                           stringsAsFactors = FALSE)
  row_df$Variable     <- LABELS[var]
  results_list[[var]] <- row_df
}

# Median [IQR] variables
for (var in CONTINUOUS_MEDIAN_IQR) {
  if (!var %in% names(copd_data)) {
    cat("Skipping (not found in data):", var, "\n"); next
  }
  row_df              <- as.data.frame(
                           t(summarise_median_iqr(copd_data, var)),
                           stringsAsFactors = FALSE)
  row_df$Variable     <- LABELS[var]
  results_list[[var]] <- row_df
}

# Stack all rows — bind_rows first, then select to avoid dropping variables
results_df <- bind_rows(results_list)

cat("\nColumns after bind_rows:\n")
print(names(results_df))

# Reorder columns
results_df <- results_df %>%
  select(Variable, Stat, Group1, Group2, `p-value`, Test)

# Rename Group1/Group2 to actual group names
names(results_df)[names(results_df) == "Group1"] <- GROUP_ORDER[1]
names(results_df)[names(results_df) == "Group2"] <- GROUP_ORDER[2]

# Store real group column names for categorical block fixing
real_group_names <- setdiff(names(results_df),
                             c("Variable", "Stat", "p-value", "Test"))

# Verify all continuous variables made it into results_df
cat("\nVariables confirmed in results_df:\n")
print(unique(results_df$Variable))



# STEP 11 — BUILD CATEGORICAL BLOCKS (HEADER + SUB-ROWS)
# What this does: for each categorical variable, builds a header row showing
# the overall count/% with a p-value, plus indented sub-rows for each
# category level. fix_group_names() ensures column names match results_df
# to prevent the NA/duplicate column issue.


cat("\nChecking category labels in your data:\n")
cat("sex:           "); print(unique(copd_data$sex))
cat("ics_use:       "); print(unique(copd_data$ics_use))
cat("oral_steroids: "); print(unique(copd_data$oral_steroids))
cat("crf_smoking:   "); print(unique(copd_data$crf_smoking))
cat("exacerbations_flag: "); print(unique(copd_data$exacerbations_flag))

# --- Sex ---
sex_header  <- make_header_row(copd_data, "sex",
                                "Sex, n (%)",
                                primary_level = "Male")
sex_subrows <- summarise_categorical_both(copd_data, "sex",
                                           level_labels = c("Male", "Female"))
sex_block   <- fix_group_names(bind_rows(sex_header, sex_subrows),
                                real_group_names)

# --- ICS use ---
ics_header  <- make_header_row(copd_data, "ics_use",
                                "ICS use, n (%)",
                                primary_level = "Yes")
ics_subrows <- summarise_categorical_both(copd_data, "ics_use",
                                           level_labels = c("Yes", "No"))
ics_block   <- fix_group_names(bind_rows(ics_header, ics_subrows),
                                real_group_names)

# --- Oral steroids ---
oral_header  <- make_header_row(copd_data, "oral_steroids",
                                 "Oral steroids, n (%)",
                                 primary_level = "Yes")
oral_subrows <- summarise_categorical_both(copd_data, "oral_steroids",
                                            level_labels = c("Yes", "No"))
oral_block   <- fix_group_names(bind_rows(oral_header, oral_subrows),
                                 real_group_names)

# --- Smoking status ---
# *** Run unique(copd_data$crf_smoking) above and adjust level_labels ***
# *** to match exactly what your data contains                        ***
smoking_header  <- make_header_row(copd_data, "crf_smoking",
                                    "Smoking status, n (%)",
                                    primary_level = "Current smoker")
smoking_subrows <- summarise_categorical_both(copd_data, "crf_smoking",
                                               level_labels = c("Current smoker",
                                                                 "Ex-smoker"))
smoking_block   <- fix_group_names(bind_rows(smoking_header, smoking_subrows),
                                    real_group_names)

# --- Number of exacerbations ---
exac_header  <- make_header_row(copd_data, "exacerbations_flag",
                                 "Number of exacerbations, n (%)",
                                 primary_level = ">0")
exac_subrows <- summarise_categorical_both(copd_data, "exacerbations_flag",
                                            level_labels = c(">0", "0"))
exac_block   <- fix_group_names(bind_rows(exac_header, exac_subrows),
                                 real_group_names)



# STEP 12 — ADD N ROW
# What this does: counts patients in each group and adds a Subjects n row
# at the top of the table. as.character() is essential — converts the
# integer count to text to match the character type of all other cells,
# preventing the bind_rows type mismatch error.


n_row <- data.frame(
  Variable         = "Subjects, n",
  Stat             = "",
  Group1           = as.character(
                       sum(copd_data$eos_class == GROUP_ORDER[1],
                           na.rm = TRUE)),
  Group2           = as.character(
                       sum(copd_data$eos_class == GROUP_ORDER[2],
                           na.rm = TRUE)),
  `p-value`        = "",
  Test             = "",
  stringsAsFactors = FALSE,
  check.names      = FALSE
)
names(n_row)[3:4] <- GROUP_ORDER


# STEP 13 — PULL INDIVIDUAL CONTINUOUS ROWS BY LABEL
# What this does: extracts each variable row from results_df by matching
# the display label. Uses LABELS["varname"] instead of retyping strings
# to guarantee exact match — eliminates the 0-observation issue.
# Run unique(results_df$Variable) if any row shows 0 observations.


# Verify all variables are present before pulling
cat("\nRow presence check:\n")
for (var in c(CONTINUOUS_MEAN_SD, CONTINUOUS_MEDIAN_IQR)) {
  present <- any(results_df$Variable == LABELS[var])
  cat(sprintf("%-25s : %s\n", var, if (present) "OK" else "MISSING — check LABELS"))
}

# Pull each row using LABELS key lookup
age_row        <- results_df[results_df$Variable == LABELS["age"], ]
bmi_row        <- results_df[results_df$Variable == LABELS["bmi"], ]
pack_row       <- results_df[results_df$Variable == LABELS["pack_years"], ]
cat_row        <- results_df[results_df$Variable == LABELS["cat_score"], ]
fev1_pct_row   <- results_df[results_df$Variable == LABELS["fev1_pct"], ]
fev1_fvc_row   <- results_df[results_df$Variable == LABELS["fev1_fvc"], ]
rv_tlc_row     <- results_df[results_df$Variable == LABELS["rv_tlc"], ]
tlco_row       <- results_df[results_df$Variable == LABELS["tlco"], ]
kco_row        <- results_df[results_df$Variable == LABELS["kco"], ]
eos_row        <- results_df[results_df$Variable == LABELS["eosinophils"], ]
neut_row       <- results_df[results_df$Variable == LABELS["neutrophils"], ]
leuk_row       <- results_df[results_df$Variable == LABELS["leukocytes"], ]
lymph_row      <- results_df[results_df$Variable == LABELS["lymphocytes"], ]
mono_row       <- results_df[results_df$Variable == LABELS["monocytes"], ]
crp_row        <- results_df[results_df$Variable == LABELS["crp"], ]
air_trap_row   <- results_df[results_df$Variable == LABELS["air_trapping"], ]
emph_row       <- results_df[results_df$Variable == LABELS["emphysema"], ]
mucus_plug_row <- results_df[results_df$Variable == LABELS["mucus_plugs"], ]
mucus_vol_row  <- results_df[results_df$Variable == LABELS["mucus_volume"], ]
pi10_row       <- results_df[results_df$Variable == LABELS["pi10"], ]
ax_row         <- results_df[results_df$Variable == LABELS["ax"], ]
x5_row         <- results_df[results_df$Variable == LABELS["x5"], ]

# Final check — all rows should show 1 observation
cat("\nObservation count per row (all should be 1):\n")
row_names <- c("age_row", "bmi_row", "pack_row", "cat_row",
               "fev1_pct_row", "fev1_fvc_row", "rv_tlc_row",
               "tlco_row", "kco_row", "eos_row", "neut_row",
               "leuk_row", "lymph_row", "mono_row", "crp_row",
               "air_trap_row", "emph_row", "mucus_plug_row",
               "mucus_vol_row", "pi10_row", "ax_row", "x5_row")

for (r in row_names) {
  n      <- nrow(get(r))
  status <- if (n == 1) "OK" else "*** PROBLEM ***"
  cat(sprintf("%-20s : %d  %s\n", r, n, status))
}


# STEP 14 — BUILD SECTION HEADERS AND ASSEMBLE FINAL TABLE
# (spirometry, diffusion, blood counts, HRCT, oscillometry), then assembles
# everything into the final results_df in your exact requested row order.


make_section_header <- function(label, col_names) {
  row           <- data.frame(matrix("", nrow = 1,
                                     ncol = length(col_names)))
  colnames(row) <- col_names
  row$Variable  <- label
  row
}

cn <- names(results_df)

section_spirometry   <- make_section_header(
                          "Post-bronchodilator spirometry", cn)
section_diffusion    <- make_section_header(
                          "Diffusion capacity", cn)
section_blood        <- make_section_header(
                          "Blood cell counts", cn)
section_hrct         <- make_section_header(
                          "High-Resolution Computed Tomography (HRCT)", cn)
section_oscillometry <- make_section_header(
                          "Impulse Oscillometry", cn)

# Assemble final table in exact requested order
results_df <- bind_rows(
  # Header
  n_row,

  # Demographics and clinical
  age_row,
  sex_block,
  bmi_row,
  pack_row,
  cat_row,
  smoking_block,
  ics_block,
  oral_block,
  exac_block,

  # Post-bronchodilator spirometry
  section_spirometry,
  fev1_pct_row,
  fev1_fvc_row,
  rv_tlc_row,

  # Diffusion capacity
  section_diffusion,
  tlco_row,
  kco_row,

  # Blood cell counts
  section_blood,
  eos_row,
  neut_row,
  leuk_row,
  lymph_row,
  mono_row,
  crp_row,

  # HRCT
  section_hrct,
  air_trap_row,
  emph_row,
  mucus_plug_row,
  mucus_vol_row,
  pi10_row,

  # Impulse oscillometry
  section_oscillometry,
  ax_row,
  x5_row
)

rownames(results_df) <- NULL

cat("\nFinal table preview:\n")
print(results_df, row.names = FALSE)



#STEP 15 EXPORT TO EXCEL


section_labels <- c(
  "Post-bronchodilator spirometry",
  "Diffusion capacity",
  "Blood cell counts",
  "High-Resolution Computed Tomography (HRCT)",
  "Impulse Oscillometry"
)

wb <- createWorkbook()
addWorksheet(wb, "COPD Final Baseline")
writeData(wb, "COPD Final Baseline", results_df)

# Blue bold header row (row 1 in Excel = column names)
addStyle(wb, "COPD Final Baseline",
  createStyle(fgFill         = "#2F5496",
              fontColour     = "#FFFFFF",
              textDecoration = "bold",
              halign         = "center"),
  rows = 1, cols = 1:ncol(results_df), gridExpand = TRUE)

# N row — bold grey (row 2 in Excel)
addStyle(wb, "COPD Final Baseline",
  createStyle(textDecoration = "bold", fgFill = "#D9D9D9"),
  rows = 2, cols = 1:ncol(results_df), gridExpand = TRUE)

# Section header rows — bold light grey
# +1 because Excel row 1 is the column header, so dataframe row 1 = Excel row 2
section_excel_rows <- which(results_df$Variable %in% section_labels) + 1
if (length(section_excel_rows) > 0) {
  addStyle(wb, "COPD Final Baseline",
    createStyle(textDecoration = "bold", fgFill = "#E8E8E8"),
    rows       = section_excel_rows,
    cols       = 1:ncol(results_df),
    gridExpand = TRUE)
}

# Alternating blue shading on data rows (skip section header rows)
for (i in seq(4, nrow(results_df) + 1, 2)) {
  if (!i %in% section_excel_rows) {
    addStyle(wb, "COPD Final Baseline",
      createStyle(fgFill = "#DCE6F1"),
      rows = i, cols = 1:ncol(results_df), gridExpand = TRUE)
  }
}

# Red bold for significant p-values (p < 0.05)
# sig_idx finds which dataframe rows have significant p-values
# +1 adjusts for Excel's header row offset
sig_idx <- which(grepl("^0\\.0[0-4]|<0\\.001", results_df$`p-value`))
if (length(sig_idx) > 0) {
  p_col <- which(names(results_df) == "p-value")
  addStyle(wb, "COPD Final Baseline",
    createStyle(fontColour     = "#C00000",
                textDecoration = "bold"),
    rows       = sig_idx + 1,
    cols       = p_col,
    gridExpand = TRUE)
}

# Auto-fit column widths
setColWidths(wb, "COPD Final Baseline",
             cols   = 1:ncol(results_df),
             widths = "auto")

# Save the workbook
saveWorkbook(wb, 
             "./data/results/Table1_COPDFinal_baseline.xlsx", overwrite = TRUE)
cat("Excel file saved: Table_COPDFinal_baseline.xlsx\n")




