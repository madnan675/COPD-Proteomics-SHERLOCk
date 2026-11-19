#PATHWAY ENRICHMENT ANALYSIS-COPD Data Only 
#COPD DATA ONLY

#PURPOSE:
  # Takes the differentially expressed proteins identified in Stage 3
  # (COPD Eosinophilic vs Non-Eosinophilic) and tests whether they cluster
  # into known biological pathways, using over-representation analysis
  # against GO Biological Process and KEGG databases.
  #
  # THREE ANALYSES:
  #   1. GO Biological Process — COPD significant proteins (FDR < 0.20)
  #   2. KEGG Pathways         — COPD significant proteins (FDR < 0.20)
  #   3. GO Biological Process — COPD-specific proteins (disease-specific,
  #                              not shared with Control's top 200)
  #
  # All three use the full protein panel (COPD_All_Proteins) as the
  # statistical background/universe, and test enrichment at p.adjust < 0.05.
  #
  # INPUT:
  #   ./data/input/DEP_Results_Stage3.xlsx
  #   (sheets: COPD_Significant, COPD_All_Proteins, COPD_Specific)
  #
  # OUTPUT:
  #   ./data/results/Pathway_Enrichment_Results_COPD.xlsx  (GO/KEGG tables)
  #   Dotplots, protein-pathway network plots, and a fold-change-coloured
  #   network plot, saved as matching .pdf/.png pairs.


#STEP 1 — INSTALL AND LOAD PACKAGES

# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("enrichplot")
# 
# install.packages("tidyverse")
# install.packages("readxl")
# install.packages("writexl")
# install.packages("ggplot2")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(tidyverse)
library(readxl)
library(writexl)
library(ggplot2)


#STEP 2 LOAD EXCEL Data COPD Significant from DEP Results Stage 3 



copd_sig_df <- read_excel(
  "./data/input/DEP_Results_Stage3.xlsx",
  sheet = "COPD_Significant"
)


all_proteins_df <- read_excel(
  "./data/input/DEP_Results_Stage3.xlsx",
  sheet = "COPD_All_Proteins"
)

copd_specific_df <- read_excel(
  "./data/input/DEP_Results_Stage3.xlsx",
  sheet = "COPD_Specific"
)



cat("COPD significant: ", nrow(copd_sig_df), "\n")
cat("Background:       ", nrow(all_proteins_df), "\n")
cat("COPD specific:    ", nrow(copd_specific_df), "\n")

head(copd_sig_df$Protein, 10)



#STEP 3 — EXTRACT PROTEIN NAME VECTORS


proteins_copd_sig      <- as.character(copd_sig_df$Protein)
proteins_copd_specific <- as.character(copd_specific_df$Protein)
proteins_background    <- as.character(all_proteins_df$Protein)

cat("COPD significant:  ", length(proteins_copd_sig), "\n")
cat("COPD specific:     ", length(proteins_copd_specific), "\n")
cat("Background:        ", length(proteins_background), "\n")



#STEP 4 — CONVERT GENE SYMBOLS TO ENTREZ IDs


map_sig <- bitr(
  proteins_copd_sig,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

unmapped_sig <- setdiff(proteins_copd_sig, map_sig$SYMBOL)
print(unmapped_sig)

map_specific <- bitr(
  proteins_copd_specific,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

map_background <- bitr(
  proteins_background,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

cat("\n=== MAPPING SUMMARY ===\n")
cat("COPD significant: ", length(proteins_copd_sig),
    "input →", nrow(map_sig), "mapped\n")
cat("COPD specific:    ", length(proteins_copd_specific),
    "input →", nrow(map_specific), "mapped\n")
cat("Background:       ", length(proteins_background),
    "input →", nrow(map_background), "mapped\n")



#STEP 5 — ANALYSIS 1: GO BIOLOGICAL PROCESS
#COPD SIGNIFICANT (217 proteins)


go_bp_copd <- enrichGO(
  gene          = map_sig$ENTREZID,
  universe      = map_background$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  minGSSize     = 5,
  maxGSSize     = 500,
  readable      = TRUE
)

n_go_bp <- nrow(go_bp_copd@result %>% filter(p.adjust < 0.05))
cat("\nGO BP significant terms:", n_go_bp, "\n")

if(n_go_bp > 0) {
  go_bp_copd@result %>%
    filter(p.adjust < 0.05) %>%
    select(Description, Count, GeneRatio, p.adjust) %>%
    mutate(p.adjust = round(p.adjust, 5)) %>%
    arrange(p.adjust) %>%
    head(15) %>%
    print()
}



#STEP 6 — ANALYSIS 2: KEGG PATHWAYS
#COPD SIGNIFICANT (217 proteins)


kegg_copd <- enrichKEGG(
  gene          = map_sig$ENTREZID,
  universe      = map_background$ENTREZID,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  minGSSize     = 5,
  maxGSSize     = 500
)

kegg_copd <- setReadable(kegg_copd,
                          OrgDb   = org.Hs.eg.db,
                          keyType = "ENTREZID")

n_kegg <- nrow(kegg_copd@result %>% filter(p.adjust < 0.05))
cat("\nKEGG significant pathways:", n_kegg, "\n")

if(n_kegg > 0) {
  kegg_copd@result %>%
    filter(p.adjust < 0.05) %>%
    select(Description, Count, GeneRatio, p.adjust) %>%
    mutate(p.adjust = round(p.adjust, 5)) %>%
    arrange(p.adjust) %>%
    head(15) %>%
    print()
}



#STEP 7 — ANALYSIS 3: GO BIOLOGICAL PROCESS
#COPD SPECIFIC (200 proteins)


go_bp_specific <- enrichGO(
  gene          = map_specific$ENTREZID,
  universe      = map_background$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  minGSSize     = 5,
  maxGSSize     = 500,
  readable      = TRUE
)

n_go_specific <- nrow(go_bp_specific@result %>% filter(p.adjust < 0.05))
cat("\nGO BP COPD specific significant terms:", n_go_specific, "\n")

if(n_go_specific > 0) {
  go_bp_specific@result %>%
    filter(p.adjust < 0.05) %>%
    select(Description, Count, GeneRatio, p.adjust) %>%
    mutate(p.adjust = round(p.adjust, 5)) %>%
    arrange(p.adjust) %>%
    head(15) %>%
    print()
}


#STEP 8 — DOTPLOT: GO BP COPD SIGNIFICANT


if(n_go_bp > 0) {

  p1 <- dotplot(
    go_bp_copd,
    showCategory = 20,
    font.size    = 9,
    title        = "GO Biological Process\nCOPD: Eosinophilic vs Non-Eosinophilic\n(FDR < 0.20, n = 217 proteins)",
    color        = "p.adjust",
    orderBy      = "GeneRatio"
  ) +
  scale_colour_gradient(
    low  = "#B2182B",
    high = "#4393C3",
    name = "FDR"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title   = element_text(face = "bold", size = 11,
                                 hjust = 0.5),
    axis.text.y  = element_text(size = 8),
    axis.text.x  = element_text(size = 9)
  )

  print(p1)

  ggsave("./data/results/Dotplot_GO_BP_COPD_Significant.pdf",
         plot = p1, width = 11, height = 10, dpi = 300)
  ggsave("./data/results/Dotplot_GO_BP_COPD_Significant.png",
         plot = p1, width = 11, height = 10, dpi = 300)

  cat("Saved: Dotplot_GO_BP_COPD_Significant.pdf and .png\n")

} else {
  cat("No significant GO BP terms at FDR 0.05\n")
  cat("Try changing pvalueCutoff to 0.10 in Step 5\n")
}



#STEP 9 — DOTPLOT: KEGG COPD SIGNIFICANT


if(n_kegg > 0) {

  p2 <- dotplot(
    kegg_copd,
    showCategory = 20,
    font.size    = 9,
    title        = "KEGG Pathway Analysis\nCOPD: Eosinophilic vs Non-Eosinophilic\n(FDR < 0.20, n = 217 proteins)",
    color        = "p.adjust",
    orderBy      = "GeneRatio"
  ) +
  scale_colour_gradient(
    low  = "#B2182B",
    high = "#4393C3",
    name = "FDR"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title   = element_text(face = "bold", size = 11,
                                 hjust = 0.5),
    axis.text.y  = element_text(size = 8),
    axis.text.x  = element_text(size = 9)
  )

  print(p2)

  ggsave("./data/results/Dotplot_KEGG_COPD_Significant.pdf",
         plot = p2, width = 11, height = 10, dpi = 300)
  ggsave("./data/results/Dotplot_KEGG_COPD_Significant.png",
         plot = p2, width = 11, height = 10, dpi = 300)

  cat("Saved: Dotplot_KEGG_COPD_Significant.pdf and .png\n")

} else {
  cat("No significant KEGG pathways at FDR 0.05\n")
}



#STEP 10 — DOTPLOT: GO BP COPD SPECIFIC


if(n_go_specific > 0) {

  p3 <- dotplot(
    go_bp_specific,
    showCategory = 20,
    font.size    = 9,
    title        = "GO Biological Process — COPD-Specific Proteins\nEosinophilic vs Non-Eosinophilic\n(n = 200 proteins, disease-specific only)",
    color        = "p.adjust",
    orderBy      = "GeneRatio"
  ) +
  scale_colour_gradient(
    low  = "#B2182B",
    high = "#4393C3",
    name = "FDR"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title   = element_text(face = "bold", size = 11,
                                 hjust = 0.5),
    axis.text.y  = element_text(size = 8),
    axis.text.x  = element_text(size = 9)
  )

  print(p3)

  ggsave("./data/results/Dotplot_GO_BP_COPD_Specific.pdf",
         plot = p3, width = 11, height = 10, dpi = 300)
  ggsave("./data/results/Dotplot_GO_BP_COPD_Specific.png",
         plot = p3, width = 11, height = 10, dpi = 300)

  cat("Saved: Dotplot_GO_BP_COPD_Specific.pdf and .png\n")

} else {
  cat("No significant GO BP terms for COPD specific at FDR 0.05\n")
}



#STEP 11 — SAVE RESULTS TO EXCEL


write_xlsx(
  list(
    "COPD_GO_BP"          = go_bp_copd@result %>%
                              arrange(p.adjust),
    "COPD_KEGG"           = kegg_copd@result %>%
                              arrange(p.adjust),
    "COPDspecific_GO_BP"  = go_bp_specific@result %>%
                              arrange(p.adjust)
  ),
  path = "./data/results/Pathway_Enrichment_Results_COPD.xlsx"
)

cat("\nResults saved: Pathway_Enrichment_Results_COPD.xlsx\n")

#Network plot — GO BP COPD significant

if (nrow(go_bp_copd@result %>% filter(p.adjust < 0.05)) > 0) {
  
  p_network <- cnetplot(
    go_bp_copd,
    showCategory = 8,
    color_edge = "grey70",
    size_edge = 0.8,
    node_label = "all"
  ) +
    labs(
      title = "Protein-Pathway Network\nCOPD: Eosinophilic vs Non-Eosinophilic"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 14,
        hjust = 0.5
      )
    )
  
  ggsave(
    "./data/results/Network_GO_BP_COPD.pdf",
    plot = p_network,
    width = 14,
    height = 12,
    dpi = 300
  )
  
  ggsave(
    "./data/results/Network_GO_BP_COPD.png",
    plot = p_network,
    width = 14,
    height = 12,
    dpi = 300
  )
  
  cat("Saved: Network_GO_BP_COPD\n")
}


# Network plot — KEGG COPD
if (nrow(kegg_copd@result %>% filter(p.adjust < 0.05)) > 0) {
  
  p_network_kegg <- cnetplot(
    kegg_copd,
    showCategory = 5,
    color_edge = "grey70",
    size_edge = 0.8,
    node_label = "all"
  ) +
    labs(
      title = "Protein-Pathway Network — KEGG\nCOPD: Eosinophilic vs Non-Eosinophilic"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 14,
        hjust = 0.5
      )
    )
  
  ggsave(
    "./data/results/Network_KEGG_COPD.pdf",
    plot = p_network_kegg,
    width = 14,
    height = 12,
    dpi = 300
  )
  
  ggsave(
    "./data/results/Network_KEGG_COPD.png",
    plot = p_network_kegg,
    width = 14,
    height = 12,
    dpi = 300
  )
  
  cat("Saved: Network_KEGG_COPD\n")
}
#Network pathway ______________________________________________________________
# First create a named vector of logFC values for your proteins
# This colours each protein node by its fold change direction

logfc_vector <- copd_sig_df$logFC
names(logfc_vector) <- copd_sig_df$Protein

# Network plot with fold change colouring
logfc_vector <- copd_sig_df$logFC
names(logfc_vector) <- copd_sig_df$Protein

# Network plot with fold-change colouring
p_net_fc <- cnetplot(
  go_bp_copd,
  showCategory = 8,
  foldChange = logfc_vector,
  node_label = "all",
  color_category = "#E8A838",
  color_edge = "grey70",
  size_edge = 0.8,
  layout = igraph::layout_with_kk
) +
  scale_colour_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Log2FC"
  ) +
  labs(
    title = "Protein-Pathway Network",
    subtitle = "COPD Eosinophilic vs Non-Eosinophilic\nNode colour = Log2 Fold Change"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5,
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      size = 11,
      hjust = 0.5,
      colour = "grey40",
      margin = margin(b = 10)
    ),
    legend.position = "right",
    legend.title = element_text(
      size = 10,
      face = "bold"
    ),
    legend.text = element_text(size = 9),
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    plot.margin = margin(20, 20, 20, 20)
  )

ggsave(
  "./data/results/Network_COPD_FoldChange.pdf",
  plot = p_net_fc,
  width = 16,
  height = 14,
  dpi = 300
)

ggsave(
  "./data/results/Network_COPD_FoldChange.png",
  plot = p_net_fc,
  width = 16,
  height = 14,
  dpi = 300,
  bg = "white"
)

cat("Saved: Network_COPD_FoldChange\n")
