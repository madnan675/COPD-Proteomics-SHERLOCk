#Merged Clinical and proteomics Data for downstream analysis

#PURPOSE:
# Merges baseline clinical data with Olink proteomics (NPX) measurements
# into a single patient-level dataset, ready for downstream differential
# expression analysis (Stage 3).
#
# WHAT IT DOES:
#   1. Loads clinical data and selects/renames the columns needed
#      downstream (ID, disease classification, eosinophils, and the
#      DEP covariates: age, sex, smoking status).
#   2. Removes patients with missing eosinophil values, since group
#      assignment depends on this.
#   3. Loads the proteomics CSV (proteins as rows, patients as columns)
#      and transposes it to patients as rows, proteins as columns.
#   4. Matches patients between the two datasets by Subject_ID and
#      inner-joins them, keeping only patients present in both.
#   5. Classifies each patient as Eosinophilic / Non-Eosinophilic
#      (blood eosinophils >= 0.3 x10^9/L).
#   6. Splits the merged dataset into COPD-only and Control-only subsets.
#   7. Saves all three (full merged, COPD-only, Control-only) as sheets
#      in one Excel file.

Install and Load Packages-----
install.packages("tidyverse")   # data wrangling and pipes
install.packages("readxl")      # reading Excel (.xlsx) files
install.packages("writexl")     # saving results back to Excel 

library(tidyverse)
library(readxl)
library(writexl)
library(dplyr)

#load baseline Clinical data
clinical_raw <- read_excel("data/input/clinical_data_with classification.xlsx")


dim(clinical_raw)        # tells you rows (patients) x columns (variables)
colnames(clinical_raw)   # prints ALL column names — verify they look right
head(clinical_raw, 3)    # shows first 3 rows to confirm data loaded correctly

clinical_clean <- clinical_raw %>%
  dplyr::select(
    Subject_ID      = "class_incl_study_id",          
    classification  = "Group",            
    eosinophils     = "blood_eusinophils",
    age             = "crf_age",
    sex             = "crf_gender",
    smoking_status  = "crf_smoking",
  )

sum(duplicated(clinical_clean$Subject_ID)) #to check theres any duplicate id or not 
table(clinical_clean$classification)  #how many copd and control patients



#Remove rows with NA in eosinophils column
clinical_no_na <- clinical_clean %>%
  filter(!is.na(eosinophils))

#remove rows with NA in eosinophils column and save the removed Subject_IDs
removed_ids <- clinical_clean %>%
  filter(is.na(eosinophils)) %>%
  pull(Subject_ID)

#load Proteomics data 
olink_raw <- read.csv("data/input/proteomics_patients_x_2395proteins.csv")

dim(olink_raw)


protein_names <- olink_raw[[1]]   #save protein names before transposing
olink_values <- olink_raw[ , -1]  #Remove the protein name column
olink_transposed <- as.data.frame(t(olink_values)) #Transpose flip rows and columns
colnames(olink_transposed) <- protein_names    #Assign protein names as column headers
olink_transposed$Subject_ID <- rownames(olink_transposed) 

rownames(olink_transposed)



olink_transposed <- olink_transposed[
  !rownames(olink_transposed) %in% c("SEO234", "SEO266", "SEO296", "SEO413"),    #remove those samples whose blood Eosinophils are NA 
]

c("SEO234", "SEO266", "SEO296", "SEO413") %in% rownames(olink_transposed) #to confirm these samples have been removed


olink_transposed$Subject_ID <- rownames(olink_transposed)
olink_transposed <- olink_transposed[, c("Subject_ID", setdiff(names(olink_transposed), "Subject_ID"))]

head(clinical_no_na$Subject_ID, 10)             #check the subject_IDs match
head(olink_transposed$Subject_ID, 10)

common_ids <- intersect(clinical_no_na$Subject_ID, olink_transposed$Subject_ID)
length(common_ids)

merged_data <- inner_join(clinical_no_na,
                          olink_transposed,
                          by = "Subject_ID")
dim(merged_data)  #to verify merged data
table(merged_data$classification)

#define eosinophilic Groups
merged_data$Eosinophils_Status <-ifelse(
  merged_data$eosinophils >=0.3,
  "Eosinophilic",
  "Non_Eosinophilic"
)




merged_data <- merged_data %>%
  relocate(Eosinophils_Status, .before = 7) 

merged_data <- merged_data %>%
  dplyr::rename(Eos_Group = Eosinophils_Status)


table(merged_data$classification, merged_data$Eos_Group)    #Check how many patients fall into each group

#SPLIT INTO TWO ANALYSIS DATASETS
Copd_proteomics <- merged_data %>%
  filter(classification == "COPD")

control_proteomics <- merged_data %>%
  filter(classification == "Control")

write_xlsx(
  list(
    "Full_Merged"    = merged_data,    # all patients, all variables
    "COPD_only"      = Copd_proteomics,      # COPD patients only
    "Control_only"   = control_proteomics    # Control patients only
  ),
  path = "./data/results/Merged_Clinical_Proteomics.xlsx"
)

cat("Files saved successfully.\n")
cat("Total merged patients:", nrow(merged_data), "\n")
cat("COPD patients:", nrow(Copd_proteomics), "\n")
cat("Control patients:", nrow(control_proteomics), "\n")
