
#BASELINE CHARACTERISTICS TABLE — CONTROL GROUP (Eosinophilic vs Non-Eosinophilic)

#PURPOSE:
#Builds a "Table 1"-style baseline characteristics table comparing the control
#group split into two subgroups by blood eosinophil count:
#   - Eosinophilic:      blood eosinophils >= 0.3 x10^9/L
#   - Non-Eosinophilic:  blood eosinophils <  0.3 x10^9/L
#
#WHAT THE SCRIPT DOES:
#   1.Loads the control group data and calculates BMI from weight/height.
#   2.Renames raw dataset columns to short, consistent internal names.
#   3.Classifies each subject as Eosinophilic / Non-Eosinophilic and flags
#      subjects with >0 exacerbations vs 0.
#   4.Runs Shapiro-Wilk tests to check normality of continuous variables.
#   5.Summarizes each variable by group:
#        - Normally distributed continuous variables -> Mean +/- SD, t-test
#        - Non-normal continuous variables            -> Median [IQR], Mann-Whitney U
#        - Categorical variables                       -> n (%), Fisher's Exact test
#      Categorical variables (sex, ICS use, oral steroids, smoking status,
#      exacerbations) get a header row (n/%, p-value) plus indented sub-rows
#      for each category level.
#   6.Assembles all rows into one final table, organized into sections:
#      demographics, post-bronchodilator spirometry, diffusion capacity,
#      blood cell counts, HRCT, and impulse oscillometry — with a subject
#      count (n) row at the top.
#   7.Exports the final table to a formatted Excel file and
#      significant p-values (<0.05) highlighted in red.

# OUTPUT:
#   ./data/results/Table1_Control_Final_baseline.xlsx

# STEP 0 — INSTALL AND LOAD PACKAGES


#packages <- c("readxl", "dplyr", "huxtable", "openxlsx")

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
invisible(lapply(packages, install_if_missing))

library(readxl)
library(dplyr)
library(huxtable)
library(openxlsx)


# STEP 1 — LOAD Control Clinical DATA AND CHECK COLUMN NAMES


control_data <- read_excel(
"data/input/control_Group.xlsx"
)


cat("Rows:", nrow(control_data), "| Columns:", ncol(control_data), "\n")
cat("Column names:\n")
print(names(control_data))
control_data <- control_data[, -617]
##Calculate BMI
control_data <- control_data %>%
  mutate(bmi =crf_kg/height^2)
summary(control_data$bmi)



# STEP 2 — RENAME COLUMNS TO INTERNAL SHORT NAMES       


control_data <- control_data %>%
  rename(
    # Demographics
    age             = "crf_age",
    sex             = "crf_gender",
    bmi             = "bmi",             
    
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
    kco            = "post_KCO_pctpred",
    tlco           = "post_TLCO_pctpred",
    
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

cat("Renamed columns:\n")
print(names(control_data))






# STEP 3 — CLASSIFY EOSINOPHILIC / NON-EOSINOPHILIC
# Threshold: >= 0.3 x10^9/L = Eosinophilic
#            <  0.3 x10^9/L = Non-Eosinophilic


control_data <- control_data %>%
  mutate(eos_class = if_else(eosinophils >= 0.3,
                              "Control-Eosinophilic",
                              "Control-Non_Eosinophilic"))

cat("\nGroup sizes:\n")
print(table(control_data$eos_class))



# STEP 4 — SET GROUP ORDER AS FACTOR


GROUP_ORDER <- c("Control-Eosinophilic", "Control-Non_Eosinophilic")

control_data <- control_data %>%
  mutate(eos_class = factor(eos_class, levels = GROUP_ORDER))



# STEP 5 — CREATE EXACERBATIONS FLAG (>0 vs 0)


control_data <- control_data %>%
  mutate(exacerbations_flag = if_else(exacerbations > 0, ">0", "0"))

cat("\nExacerbations flag distribution:\n")
print(table(control_data$exacerbations_flag))



# STEP 6 — DEFINE VARIABLE LISTS AND DISPLAY LABELS


# Normally distributed continuous variables — reported as Mean +/- SD
CONTINUOUS_MEAN_SD <- c("age", "bmi", "fev1_pct", "fev1_fvc",
                         "rv_tlc", "tlco", "kco", "pi10", "cat_score")

# Non-normally distributed continuous variables — reported as Median [IQR]
CONTINUOUS_MEDIAN_IQR <- c("pack_years", "eosinophils", "neutrophils",
                            "leukocytes", "lymphocytes", "monocytes", "crp",
                            "air_trapping", "emphysema", "mucus_plugs",
                            "mucus_volume", "ax", "x5")

# Categorical variables with sub-rows
CATEGORICAL_WITH_SUBROWS <- c("sex", "ics_use", "oral_steroids",
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
  tlco               = "TLCO %pred",
  kco                = "KCO %pred",
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


all_continuous <- c(CONTINUOUS_MEAN_SD, CONTINUOUS_MEDIAN_IQR)

cat("\nShapiro-Wilk Normality Test Results:\n")
cat(strrep("-", 65), "\n")

for (var in all_continuous) {
  if (!var %in% names(control_data)) next
  x <- na.omit(control_data[[var]])
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
fmt_p <- function(p) {
  if (is.na(p)) return("—")
  if (p < 0.001) return("<0.001")
  formatC(p, digits = 3, format = "f")
}

# 9b — Mean +/- SD summary with independent t-test
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

# 9c — Median [IQR] summary with Mann-Whitney U test
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

# 9d — Detect positive category (Male/Female, Yes/No, 1/0)
positive_level <- function(x) {
  vals <- unique(na.omit(as.character(x)))
  pos  <- c("1", "yes", "y", "true", "male", "m")
  hit  <- vals[tolower(vals) %in% pos]
  if (length(hit) > 0) return(hit[1])
  sort(vals)[length(vals)]
}

# 9e — Categorical header row with real count/% and Fisher's Exact p-value
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

# 9f — Categorical sub-rows (indented, blank p-value)
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
fix_group_names <- function(block, real_group_names) {
  names(block)[names(block) == "Group1"] <- real_group_names[1]
  names(block)[names(block) == "Group2"] <- real_group_names[2]
  block
}



# STEP 10 — BUILD CONTINUOUS VARIABLE ROWS


results_list <- list()

# Mean +/- SD variables
for (var in CONTINUOUS_MEAN_SD) {
  if (!var %in% names(control_data)) {
    cat("Skipping (not found):", var, "\n"); next
  }
  row_df              <- as.data.frame(t(summarise_mean_sd(control_data, var)),
                                       stringsAsFactors = FALSE)
  row_df$Variable     <- LABELS[var]
  results_list[[var]] <- row_df
}

# Median [IQR] variables
for (var in CONTINUOUS_MEDIAN_IQR) {
  if (!var %in% names(control_data)) {
    cat("Skipping (not found):", var, "\n"); next
  }
  row_df              <- as.data.frame(t(summarise_median_iqr(control_data, var)),
                                       stringsAsFactors = FALSE)
  row_df$Variable     <- LABELS[var]
  results_list[[var]] <- row_df
}

# Stack all continuous rows into one table
results_df <- bind_rows(results_list) %>%
  select(Variable, Stat, Group1, Group2, `p-value`, Test)

# Rename Group1/Group2 to actual group names
names(results_df)[names(results_df) == "Group1"] <- GROUP_ORDER[1]
names(results_df)[names(results_df) == "Group2"] <- GROUP_ORDER[2]

# Store real group column names for use in categorical blocks
real_group_names <- setdiff(names(results_df),
                             c("Variable", "Stat", "p-value", "Test"))



# STEP 11 — BUILD CATEGORICAL BLOCKS (HEADER + SUB-ROWS)


# Check actual category labels in your data
cat("\nUnique values — sex:           "); print(unique(control_data$sex))
cat("Unique values — ics_use:       "); print(unique(control_data$ics_use))
cat("Unique values — oral_steroids: "); print(unique(control_data$oral_steroids))
cat("Unique values — crf_smoking:   "); print(unique(control_data$crf_smoking))
cat("Unique values — exacerbations_flag: "); print(unique(control_data$exacerbations_flag))

# --- Sex (Male / Female) ---
sex_header  <- make_header_row(control_data, "sex",
                                "Sex, n (%)", primary_level = "Male")
sex_subrows <- summarise_categorical_both(control_data, "sex",
                                           level_labels = c("Male", "Female"))
sex_block   <- fix_group_names(bind_rows(sex_header, sex_subrows),
                                real_group_names)

# --- ICS use (Yes / No) ---
ics_header  <- make_header_row(control_data, "ics_use",
                                "ICS use, n (%)", primary_level = "Yes")
ics_subrows <- summarise_categorical_both(control_data, "ics_use",
                                           level_labels = c("Yes", "No"))
ics_block   <- fix_group_names(bind_rows(ics_header, ics_subrows),
                                real_group_names)

# --- Oral steroids (Yes / No) ---
oral_header  <- make_header_row(control_data, "oral_steroids",
                                 "Oral steroids, n (%)", primary_level = "Yes")
oral_subrows <- summarise_categorical_both(control_data, "oral_steroids",
                                            level_labels = c("Yes", "No"))
oral_block   <- fix_group_names(bind_rows(oral_header, oral_subrows),
                                 real_group_names)

# --- Smoking status (crf_smoking) ---
# *** Common values: "Current smoker, Ex-smoker"
smoking_header  <- make_header_row(control_data, "crf_smoking",
                                    "Smoking status, n (%)",
                                    primary_level = "Current smoker")
smoking_subrows <- summarise_categorical_both(control_data, "crf_smoking",
                                               level_labels = c("Current smoker",
                                                                 "Ex-smoker"))
smoking_block   <- fix_group_names(bind_rows(smoking_header, smoking_subrows),
                                    real_group_names)

# --- Number of exacerbations (>0 / 0) ---
exac_header  <- make_header_row(control_data, "exacerbations_flag",
                                 "Number of exacerbations, n (%)",
                                 primary_level = ">0")
exac_subrows <- summarise_categorical_both(control_data, "exacerbations_flag",
                                            level_labels = c(">0", "0"))
exac_block   <- fix_group_names(bind_rows(exac_header, exac_subrows),
                                 real_group_names)



# STEP 12 — ADD N ROW


n_row <- data.frame(
  Variable         = "Subjects, n",
  Stat             = "",
  Group1           = as.character(sum(control_data$eos_class == GROUP_ORDER[1],
                                      na.rm = TRUE)),
  Group2           = as.character(sum(control_data$eos_class == GROUP_ORDER[2],
                                      na.rm = TRUE)),
  `p-value`        = "",
  Test             = "",
  stringsAsFactors = FALSE,
  check.names      = FALSE
)
names(n_row)[3:4] <- GROUP_ORDER


cat("Columns in results_df after step 10:\n")
print(names(results_df))



# STEP 13 — PULL INDIVIDUAL CONTINUOUS ROWS BY DISPLAY LABEL


age_row        <- results_df[results_df$Variable == "Age, years", ]
bmi_row        <- results_df[results_df$Variable == "BMI, kg/m2", ]
pack_row       <- results_df[results_df$Variable == "Pack-years", ]
cat_row        <- results_df[results_df$Variable == "CAT total score", ]
fev1_pct_row   <- results_df[results_df$Variable == "FEV1% pred", ]
fev1_fvc_row   <- results_df[results_df$Variable == "FEV1/FVC%", ]
rv_tlc_row     <- results_df[results_df$Variable == "RV/TLC %pred",] 
tlco_row       <- results_df[results_df$Variable == "TLCO %pred",] 
kco_row        <- results_df[results_df$Variable == "KCO %pred",]
eos_row        <- results_df[results_df$Variable == "Eosinophils x10^9/L",]
neut_row       <- results_df[results_df$Variable == "Neutrophils x10^9/L", ]
leuk_row       <- results_df[results_df$Variable == "Leukocytes x10^9/L", ]
lymph_row      <- results_df[results_df$Variable == "Lymphocytes x10^9/L", ]
mono_row       <- results_df[results_df$Variable == "Monocytes x10^9/L", ]
crp_row        <- results_df[results_df$Variable == "CRP, mg/L", ]
air_trap_row   <- results_df[results_df$Variable == "Air trapping, %", ]
emph_row       <- results_df[results_df$Variable == "Emphysema, % (LAA < -950 HU)", ]
mucus_plug_row <- results_df[results_df$Variable == "Mucus plugs count", ]
mucus_vol_row  <- results_df[results_df$Variable == "Mucus volume, mL", ]
pi10_row       <- results_df[results_df$Variable == "Pi10, mm", ]
ax_row         <- results_df[results_df$Variable == "AX, kPa.s/L", ]
x5_row         <- results_df[results_df$Variable == "X5, kPa.s/L", ]

names(results_list)


# STEP 14 — BUILD SECTION HEADERS AND ASSEMBLE FINAL TABLE


make_section_header <- function(label, col_names) {
  row           <- data.frame(matrix("", nrow = 1, ncol = length(col_names)))
  colnames(row) <- col_names
  row$Variable  <- label
  row
}

cn <- names(results_df)

section_spirometry   <- make_section_header("Post-bronchodilator spirometry", cn)
section_diffusion    <- make_section_header("Diffusion capacity", cn)
section_blood        <- make_section_header("Blood cell counts", cn)
section_hrct         <- make_section_header(
                          "High-Resolution Computed Tomography (HRCT)", cn)
section_oscillometry <- make_section_header("Impulse Oscillometry", cn)

# Assemble in final order matching your requested table layout
results_df <- bind_rows(

  # Header
  n_row,

  # Demographics
  age_row,
  sex_block,
  bmi_row,
  pack_row,
  cat_row,

  # Categorical clinical variables
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


#STEP — EXPORT TO EXCEL


wb <- createWorkbook()
addWorksheet(wb, "Control Baseline")
writeData(wb, "Control Baseline", results_df)

# Blue bold header
addStyle(wb, "Control Baseline",
  createStyle(fgFill         = "#2F5496",
              fontColour     = "#FFFFFF",
              textDecoration = "bold",
              halign         = "center"),
  rows = 1, cols = 1:ncol(results_df), gridExpand = TRUE)

# N row — bold grey
addStyle(wb, "Control Baseline",
  createStyle(textDecoration = "bold", fgFill = "#D9D9D9"),
  rows = 2, cols = 1:ncol(results_df), gridExpand = TRUE)

# Section header rows — bold light grey
section_excel_rows <- which(results_df$Variable %in% LABELS) + 1
if (length(section_excel_rows) > 0) {
  addStyle(wb, "Control Baseline",
    createStyle(textDecoration = "bold", fgFill = "#E8E8E8"),
    rows       = section_excel_rows,
    cols       = 1:ncol(results_df),
    gridExpand = TRUE)
}


# Alternating shading on data rows
for (i in seq(4, nrow(results_df) + 1, 2)) {
  if (!i %in% section_excel_rows) {
    addStyle(wb, "Control Baseline",
      createStyle(fgFill = "#DCE6F1"),
      rows = i, cols = 1:ncol(results_df), gridExpand = TRUE)
  }
}

# Red bold for significant p-values
sig_idx <- which(grepl("^0\\.0[0-4]|<0\\.001", results_df$`p-value`))
if (length(sig_idx) > 0) {
  p_col2 <- which(names(results_df) == "p-value")
  addStyle(wb, "Control Baseline",
    createStyle(fontColour = "#C00000", textDecoration = "bold"),
    rows       = sig_idx + 1,
    cols       = p_col2,
    gridExpand = TRUE)
}

setColWidths(wb, "Control Baseline",
             cols   = 1:ncol(results_df),
             widths = "auto")


saveWorkbook(
  wb,
  "./data/results/Table1_ControlFinal_baseline.xlsx",
  overwrite = TRUE
)

cat("Excel file saved: Table1_ControlFinal_baseline.xlsx\n")











