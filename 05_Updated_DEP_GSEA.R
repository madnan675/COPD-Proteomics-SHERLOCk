
#BOXPLOTS, VOLCANO PLOT & GSEA — COPD Eosinophilic Analysis

#PURPOSE:
# Visualises the Stage 3 differential expression results and runs
# Gene Set Enrichment Analysis (GSEA) on the full ranked protein list.
#
# WHAT IT DOES:
#   1. Boxplots of the top 5 most significant COPD proteins (NPX by
#      Eosinophilic/Non-Eosinophilic group).
#   2. Volcano plot of all COPD proteins (FDR < 0.20, |logFC| > 0.5),
#      with the top 25 proteins labelled.
#   3. Ranks all 2,395 proteins by signed -log10(p-value) and maps gene
#      symbols to Entrez IDs for GSEA.
#   4. Runs GSEA against GO Biological Process and KEGG pathways using
#      the full ranked list (no significance cutoff on the input).
#   5. Produces dotplots and per-pathway enrichment plots for the
#      significant GO terms (p.adjust < 0.05).
#   6. Saves all GSEA results to Excel.
# Input:
#      Merged Clinical proteomics and DEP results for GSEA
# OUTPUT:
#     ./data/results/Boxplots_COPD_Top5.png, 
#     ./data/results/Volcano_COPD_Eos_vs_NonEos.png, 
#     ./data/results/GSEA_Dotplot_GO_BP.png,
#     ./data/results/GSEA_Enrichment_hydrogen_peroxide_catabolic_pr.png, 
#     ./data/results/GSEA_Enrichment_hydrogen_peroxide_metabolic_pr.png,
#     ./data/results/GSEA_Enrichment_protein_containing_complex_rem.png
#     ./data/results/GSEA_Results.xlsx




#STEP 1 — LOAD PACKAGES


library(tidyverse)
library(readxl)
library(writexl)
library(ggplot2)
library(ggrepel)
library(limma)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")


library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(fgsea)



#load merged clinical data for boxplots
merged_data  <- read_excel("data/input/Merged_Clinical_Proteomics.xlsx",
                           sheet = "Full_Merged")

copd_data    <- read_excel("data/input/Merged_Clinical_Proteomics.xlsx",
                           sheet = "COPD_only")

control_data <- read_excel("data/input/Merged_Clinical_Proteomics.xlsx",
                           sheet = "Control_only")

results_copd <- read_excel("data/results/DEP_Results_Stage3.xlsx",
                           sheet = "COPD_All_Proteins")
sig_copd     <- read_excel("data/results/DEP_Results_Stage3.xlsx",
                           sheet = "COPD_Significant")

#STEP 2 — BOXPLOTS FOR TOP 5 COPD PROTEINS

#Top 5 most significant COPD proteins by FDR
#Shows NPX distribution in Eosinophilic vs Non-Eosinophilic
#Using COPD patients only

# Identify protein columns
# Adjust column number if your clinical columns differ
all_proteins <- colnames(merged_data)[8:ncol(merged_data)]

#Convert protein columns to numeric
copd_data <- copd_data %>%
  mutate(across(all_of(all_proteins), as.numeric))

# Top 5 proteins
top5 <- sig_copd %>%
  arrange(adj.P.Val) %>%
  head(5) %>%
  pull(Protein)

cat("\nTop 5 COPD proteins for boxplots:\n")
print(top5)

library(ggplot2)
# Build statistics label for each panel
stats_labels <- sig_copd %>%
  filter(Protein %in% top5) %>%
  mutate(
    panel_label = paste0(
      Protein, "\n",
      "FDR = ", formatC(adj.P.Val, format = "e", digits = 2),
      "  logFC = ", round(logFC, 2)
    )
  ) %>%
  dplyr::select(Protein, panel_label)

#Reshape to long format for ggplot
copd_long <- copd_data %>%
  dplyr::select(Subject_ID, Eos_Group, all_of(top5)) %>%
  pivot_longer(
    cols      = all_of(top5),
    names_to  = "Protein",
    values_to = "NPX"
  ) %>%
  left_join(stats_labels, by = "Protein") %>%
  mutate(
    Eos_Group = factor(Eos_Group,
                       levels = c("Eosinophilic",
                                  "Non_Eosinophilic"))
  )

# Draw boxplots
p_box <- ggplot(copd_long,
                aes(x = Eos_Group, y = NPX, fill = Eos_Group)) +
  
  geom_boxplot(
    outlier.shape = NA,
    alpha         = 0.70,
    width         = 0.50,
    colour        = "grey25"
  ) +
  
  geom_jitter(
    width  = 0.15,
    size   = 1.8,
    alpha  = 0.65,
    colour = "grey20"
  ) +
  
  facet_wrap(~ panel_label,
             scales = "free_y",
             ncol   = 3) +
  
  scale_fill_manual(values = c(
    "Eosinophilic"     = "#C0392B",
    "Non_Eosinophilic" = "#2471A3"
  )) +
  
  scale_x_discrete(labels = c(
    "Eosinophilic"     = paste0("Eosinophilic\n(n=28)"),
    "Non_Eosinophilic" = paste0("Non-Eosinophilic\n(n=147)")
  )) +
  
  labs(
    title    = "COPD: Top 5 Most Significant Proteins",
    subtitle = "FDR < 0.20  |  Adjusted for age, sex, smoking status",
    x        = NULL,
    y        = "NPX (Normalised Protein Expression)",
    fill     = "Group"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    strip.text       = element_text(size = 8.5, face = "bold"),
    strip.background = element_rect(fill  = "grey95",
                                    colour = "grey70"),
    axis.text.x      = element_text(angle = 15, hjust = 1,
                                    size  = 9),
    legend.position  = "none",
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    panel.grid.minor = element_blank()
  )

print(p_box)

ggsave("Boxplots_COPD_Top5.pdf",
       plot = p_box, width = 13, height = 9, dpi = 300)
ggsave("./data/results/Boxplots_COPD_Top5.png",
       plot = p_box, width = 13, height = 9, dpi = 300,
       bg   = "white")


cat("Saved: Boxplots_COPD_Top5.pdf/.png\n")



# STEP 2B — VOLCANO PLOT: COPD Eosinophilic vs Non-Eosinophilic

# X-axis: log2 fold change (positive = higher in Eosinophilic)
# Y-axis: -log10(FDR-adjusted p-value)
# Dashed lines mark thresholds: FDR = 0.20, |log2FC| = 0.5
# Top 25 proteins by FDR are labelled

volcano_data <- results_copd %>%
  mutate(
    neg_log10_fdr = -log10(adj.P.Val),
    significant = case_when(
      adj.P.Val < 0.20 & logFC >  0 ~ "Up in Eosinophilic",
      adj.P.Val < 0.20 & logFC < 0 ~ "Down in Eosinophilic",
      TRUE                            ~ "Not Significant"
    )
  )

# Counts for the subtitle — computed from the data, not hardcoded
n_total <- nrow(volcano_data)
n_up    <- sum(volcano_data$significant == "Up in Eosinophilic")
n_down  <- sum(volcano_data$significant == "Down in Eosinophilic")

# Top 25 proteins by FDR — labelled only, not a separate filter of the plot
top25_labels <- volcano_data %>%
  arrange(adj.P.Val) %>%
  head(25)

p_volcano <- ggplot(volcano_data,
                    aes(x = logFC, y = neg_log10_fdr, colour = significant)) +
  
  geom_point(alpha = 0.65, size = 1.8) +
  
  geom_hline(yintercept = -log10(0.20),
             linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = c(-0.5, 0.5),
             linetype = "dashed", colour = "grey40") +
  
  annotate("text", x = Inf, y = -log10(0.20),
           label = "FDR = 0.20",
           hjust = 1.1, vjust = -0.6,
           size = 3.5, colour = "grey30") +
  
  geom_text_repel(
    data         = top25_labels,
    aes(label = Protein),
    size         = 3,
    max.overlaps = Inf,
    box.padding  = 0.4,
    segment.size = 0.3,
    show.legend  = FALSE
  ) +
  
  scale_colour_manual(
    breaks = c("Down in Eosinophilic", "Not Significant", "Up in Eosinophilic"),
    values = c(
      "Up in Eosinophilic"   = "#C0392B",
      "Down in Eosinophilic" = "#2471A3",
      "Not Significant"      = "grey70"
    )
  ) +
  
  labs(
    title    = "COPD: Eosinophilic vs Non-Eosinophilic",
    subtitle = paste0(
      "Total proteins: ", format(n_total, big.mark = ","),
      "  |  Upregulated: ", n_up,
      "  |  Downregulated: ", n_down,
      "  |  FDR < 0.20, |logFC| > 0.5  |  Top 25 labelled"
    ),
    x        = expression(Log[2]~"Fold Change  (Eosinophilic vs Non-Eosinophilic)"),
    y        = expression(-Log[10]~"(FDR-adjusted P-value)"),
    colour   = NULL
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 15),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

print(p_volcano)

ggsave("./data/results/Volcano_COPD_Eos_vs_NonEos.pdf",
       plot = p_volcano, width = 12, height = 9, dpi = 300)
ggsave("./data/results/Volcano_COPD_Eos_vs_NonEos.png",
       plot = p_volcano, width = 12, height = 9, dpi = 300, bg = "white")

cat("Saved: Volcano_COPD_Eos_vs_NonEos.pdf/.png\n")


#STEP3 — RANKED PROTEIN LIST FOR GSEA

#GSEA uses ALL proteins ranked by significance
#Standard ranking metric: signed -log10(p-value)
#positive = upregulated in Eosinophilic
#negative = downregulated in Eosinophilic

#This metric captures both statistical significance (p-value)and direction (sign of logFC) simultaneously

copd_ranked <- results_copd %>%
  arrange(P.Value) %>%
  mutate(
    #Ranking metric: direction × -log10(p-value)
    rank_metric = sign(logFC) * (-log10(P.Value))
  )

# Create named numeric vector: names = protein, value = rank
ranked_vector <- copd_ranked$rank_metric
names(ranked_vector) <- copd_ranked$Protein

# Sort in DECREASING order (required by fgsea and clusterProfiler)
ranked_vector <- sort(ranked_vector, decreasing = TRUE)

#CHECK ---
cat("\n=== GSEA RANKED LIST ===\n")
cat("Total proteins in ranked list: ", length(ranked_vector), "\n")
cat("Top 10 (most upregulated):\n")
print(head(ranked_vector, 10))
cat("Bottom 10 (most downregulated):\n")
print(tail(ranked_vector, 10))



#STEP 4 — CONVERT PROTEIN NAMES TO ENTREZ IDs FOR GSEA

library(org.Hs.eg.db)
# Convert names of ranked vector to Entrez IDs
gene_map <- bitr(
  names(ranked_vector),
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

AnnotationDbi::keytypes(org.Hs.eg.db)

cat("Proteins mapped for GSEA: ", nrow(gene_map),
    "out of", length(ranked_vector), "\n")

#Filter ranked vector to only mapped proteins and rename with Entrez IDs
ranked_entrez <- ranked_vector[names(ranked_vector) %in%
                                 gene_map$SYMBOL]

#Replace gene symbol names with Entrez ID names
gene_map_match <- gene_map[match(names(ranked_entrez),
                                 gene_map$SYMBOL), ]
names(ranked_entrez) <- gene_map_match$ENTREZID


cat("Final GSEA ranked list length: ", length(ranked_entrez), "\n")



#STEP 5 — RUN GSEA: GO BIOLOGICAL PROCESS

# gseGO runs GSEA against GO Biological Process database
# Uses your full ranked protein list — no cutoff
# nPerm = number of permutations for significance testing
# minGSSize / maxGSSize = gene set size limits
# seed for reproducibility


library(clusterProfiler)
library(org.Hs.eg.db)
set.seed(123)
gsea_go_bp <- gseGO(
  geneList      = ranked_entrez,   # full ranked list
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = FALSE,
  by            = "fgsea",         # uses fgsea algorithm internally
  eps           = 0,
  seed          = TRUE    
)


exists("gsea_go_bp")
#Convert IDs to readable gene symbols in output
gsea_go_bp <- setReadable(gsea_go_bp,
                          OrgDb   = org.Hs.eg.db,
                          keyType = "ENTREZID")

library(dplyr)
n_gsea_go <- nrow(gsea_go_bp@result %>% filter(p.adjust < 0.05))
cat("\nGSEA GO BP significant terms: ", n_gsea_go, "\n")

if(n_gsea_go > 0) {
  gsea_go_bp@result %>%
    filter(p.adjust < 0.05) %>%
    dplyr::select(Description, setSize, enrichmentScore,
           NES, p.adjust) %>%
    mutate(
      NES       = round(NES, 3),
      p.adjust  = formatC(p.adjust, format = "e", digits = 2)
    ) %>%
    arrange(desc(NES)) %>%
    head(20) %>%
    print()
}



#STEP 6 — RUN GSEA: KEGG PATHWAYS



library(clusterProfiler)

set.seed(123)
gsea_kegg <- gseKEGG(
  geneList      = ranked_entrez,
  organism      = "hsa",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = FALSE,
  eps           = 0,
  seed          =TRUE
)

library(org.Hs.eg.db)
gsea_kegg <- setReadable(gsea_kegg,
                         OrgDb   = org.Hs.eg.db,
                         keyType = "ENTREZID")

n_gsea_kegg <- nrow(gsea_kegg@result %>% filter(p.adjust < 0.05))
cat("\nGSEA KEGG significant pathways: ", n_gsea_kegg, "\n")




#STEP 7 — GSEA VISUALISATION: DOTPLOT

#NES = Normalised Enrichment Score
#Positive NES = pathway enriched in UPREGULATED proteins
#Negative NES = pathway enriched in DOWNREGULATED proteins
library(ggplot2)

if(n_gsea_go > 0) {
  
  p_gsea_go <- dotplot(
    gsea_go_bp,
    showCategory = 20,
    split        = ".sign",   # split by activated vs suppressed
    font.size    = 9,
    title        = "GSEA — GO Biological Process\nCOPD: Eosinophilic vs Non-Eosinophilic\n(All 2,395 proteins ranked by significance)"
  ) +
    facet_grid(. ~ .sign) +
    scale_colour_gradient(
      low  = "#B2182B",
      high = "#4393C3",
      name = "FDR"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title   = element_text(face = "bold", size = 11,
                                  hjust = 0.5),
      axis.text.y  = element_text(size = 8),
      strip.text   = element_text(face = "bold", size = 10)
    )
  
  ggsave(
    "./data/results/GSEA_Dotplot_GO_BP.png",
         plot = p_gsea_go, 
         width = 14, 
         height = 10,
         dpi = 300,
         bg = "white"
    )
  
  cat("Saved: GSEA_Dotplot_GO_BP.pdf/.png\n")
}
  

#STEP 8 - GSEA ENRICHMENT PLOT FOR TOP PATHWAYS

# The classic GSEA enrichment plot shows the running
# enrichment score across your ranked protein list
# One plot per pathway — shows where in the ranked list
# the pathway proteins cluster

library(enrichplot)
library(ggplot2)
library(ggplot2)
library(tidyverse)
#Output location
FIG_DIR <- "./data/results"          
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

if(n_gsea_go > 0) {
  
  #Check which GO terms have valid NES scores
  #Filter out any terms with NA in NES or enrichmentScore
  valid_go_terms <- gsea_go_bp@result %>%
    filter(p.adjust < 0.05) %>%
    filter(!is.na(NES)) %>%
    filter(!is.na(enrichmentScore)) %>%
    arrange(desc(abs(NES)))
  
  cat("Valid GO terms for enrichment plots:", nrow(valid_go_terms), "\n")
  
  if(nrow(valid_go_terms) > 0) {
    
    #Take top 3 valid terms
    top_terms <- valid_go_terms %>%
      head(3) %>%
      pull(ID)    # use GO ID not Description to avoid name issues
    
    cat("Plotting enrichment for:\n")
    print(valid_go_terms %>% head(3) %>%
            dplyr::select(ID, Description, NES, p.adjust))
    
    #Plot each term safely with tryCatch
    for(term_id in top_terms) {
      
      term_desc <- gsea_go_bp@result %>%
        filter(ID == term_id) %>%
        pull(Description)
      
      tryCatch({
        
        p_enrich <- gseaplot2(
          gsea_go_bp,
          geneSetID    = term_id,    # use ID not Description
          title        = paste0("GSEA: ", term_desc),
          color        = "#C0392B",
          pvalue_table = TRUE,
          ES_geom      = "line"
        )
        
        clean_name <- gsub("[^a-zA-Z0-9]", "_",
                           substr(term_desc, 1, 30))
      
        
        ggsave(
          file.path(FIG_DIR,paste0("GSEA_Enrichment_", clean_name, ".png")),
          plot   = p_enrich,
          width  = 10,
          height = 7,
          dpi    = 300,
          bg     = "white"
        )
        
        cat("Saved enrichment plot:", term_desc, "\n")
        
      }, error = function(e) {
        cat("Could not plot:", term_desc, "\n")
        cat("Error:", conditionMessage(e), "\n")
        cat("Skipping to next term.\n\n")
      })
    }
    
  } else {
    cat("No valid GO terms with complete NES scores found.\n")
  }
}

library(dplyr)
valid_go_terms %>%
  head(3) %>%
  dplyr::select(Description, NES, pvalue, p.adjust) %>%
  mutate(
    NES      = round(NES, 3),
    pvalue   = formatC(pvalue, format = "e", digits = 2),
    p.adjust = formatC(p.adjust, format = "e", digits = 3)
  ) %>%
  as.data.frame() %>%
  print()


#Step 9 - SAVE ALL GSEA RESULTS TO EXCEL

library(writexl)

write_xlsx(
  list(
    "GSEA_GO_BP" = gsea_go_bp@result %>%
      dplyr::select(
        Description,
        setSize,
        enrichmentScore,
        NES,
        pvalue,
        p.adjust,
        core_enrichment
      ) %>%
      dplyr::arrange(p.adjust)
  ),
  path = "./data/results/GSEA_Results.xlsx"
)

cat("\nGSEA results saved: GSEA_Results.xlsx\n")

