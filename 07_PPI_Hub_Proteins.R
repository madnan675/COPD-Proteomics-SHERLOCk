#PROTEIN-PROTEIN INTERACTION (PPI) NETWORK HUB PROTEIN ANALYSIS
#PURPOSE:
  # Loads a STRING database interaction network (exported as TSV) and
  # identifies "hub" proteins those with the most connections (highest
  # degree) in the network then annotates them with their COPD
  # differential expression results (logFC, adjusted p-value).
  #
  # WHAT IT DOES:
  #   1. Loads STRING interactions and counts how often each protein
  #      appears across both interaction columns (its degree/connectivity).
  #   2. Ranks proteins by degree to surface the top hub proteins.
  #   3. Joins in logFC and FDR from the COPD significant results, so hub
  #      proteins can be cross-referenced against differential expression.
  #   4. Saves the annotated hub protein table to Excel.
  #
  # OUTPUT:
  #   ./data/results/PPI_Hub_Proteins.xlsx

library(dplyr)
library(readxl)
library(writexl)

#Load the String database interactions
string_int <- read.delim("data/input/string_interactions.tsv",
                         header = TRUE,
                         stringsAsFactors = FALSE)


colnames(string_int)

#count apperances across both columns
hub_counts <- c(string_int$X.node1, string_int$node2) %>%
  table() %>%
  as.data.frame() %>%
  dplyr::rename(Protein = 1, Degree = 2) %>%
  dplyr::arrange(desc(Degree))

#Top 10 hub proteins
print(head(hub_counts, 10))


#Add Fold Change and FDR from the DEP results

copd_sig<- read_excel(
  "./data/input/DEP_Results_Stage3.xlsx",
  sheet = "COPD_Significant"
)

hub_annotated <- hub_counts %>%
  left_join(copd_sig %>% dplyr::select(Protein, logFC, adj.P.Val),
            by = "Protein") %>%
  dplyr::arrange(desc(Degree))


print(head(hub_annotated, 10))

write_xlsx(
  list("Hub_Proteins" = hub_annotated),
  path = "./data/results/PPI_Hub_Proteins.xlsx"
)
cat("Saved: PPI_Hub_Proteins.xlsx")



