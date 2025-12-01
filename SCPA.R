
library(SCPA)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(statmod)
library(dplyr)
library(msigdbr)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Check and convert gene IDs if needed
print("=== Converting Ensembl IDs to Gene Symbols ===")
current_genes <- rownames(GSE164897)
print(paste("Current genes in Seurat object:", length(current_genes)))
print("Sample genes:")
print(head(current_genes, 10))

# Check if Ensembl IDs
is_ensembl <- grepl("^ENSG", head(current_genes, 100))
ensembl_pct <- sum(is_ensembl) / length(is_ensembl) * 100

if (ensembl_pct > 50) {
  print(paste("Detected", round(ensembl_pct, 1), "% Ensembl IDs - converting to gene symbols..."))
  
  # Join layers first if needed
  tryCatch({
    GSE164897[["RNA"]] <- JoinLayers(GSE164897[["RNA"]])
    print("Joined RNA layers")
  }, error = function(e) {
    print("Layers already joined or not needed")
  })
  
  # Map Ensembl to symbols
  ensembl_ids <- rownames(GSE164897)
  gene_symbols <- mapIds(x = org.Hs.eg.db,
                         keys = ensembl_ids,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
  
  # Make unique and remove NAs
  gene_symbols <- make.unique(as.character(gene_symbols))
  
  # Update rownames
  rownames(GSE164897) <- gene_symbols
  
  # Remove NA genes
  GSE164897 <- GSE164897[!is.na(rownames(GSE164897)), ]
  
  print(paste("Conversion complete. Genes after conversion:", nrow(GSE164897)))
  print("New gene names:")
  print(head(rownames(GSE164897), 10))
} else {
  print("Gene symbols already present - no conversion needed")
}

DimPlot(GSE164897, split.by = "treatment") +
  theme(aspect.ratio = 1)



untreated <- seurat_extract(GSE164897,
                            meta1 = "treatment", value_meta1 = "untreated")

vem_cob <- seurat_extract(GSE164897,
                          meta1 = "treatment", value_meta1 = "vem_cob")

vem_tram <- seurat_extract(GSE164897,
                           meta1 = "treatment", value_meta1 = "vem_tram")

vem <- seurat_extract(GSE164897,
                      meta1 = "treatment", value_meta1 = "Vemurafenib")

# Debug: Check structure of extracted samples
print("Checking extracted sample structure:")
print(paste("untreated class:", class(untreated)))
print(paste("untreated dimensions:", paste(dim(untreated), collapse = " x ")))
print(paste("vem class:", class(vem)))
print(paste("vem dimensions:", paste(dim(vem), collapse = " x ")))

# Debug: Check gene names in Seurat object
print("\n=== Checking gene names in GSE164897 ===")
seurat_genes <- rownames(GSE164897)
print(paste("Total genes in Seurat object:", length(seurat_genes)))
print("First 10 genes in Seurat object:")
print(head(seurat_genes, 10))

# Check if they're Ensembl IDs or gene symbols
is_ensembl <- grepl("^ENSG", head(seurat_genes, 100))
print(paste("Percentage Ensembl IDs:", round(sum(is_ensembl) / length(is_ensembl) * 100, 1), "%"))

# Load Hallmark pathways from msigdbr  
hallmark_df <- msigdbr(species = "Homo sapiens", 
                       category = "H")

print(paste("\nTotal Hallmark gene sets retrieved:", length(unique(hallmark_df$gs_name))))

# Use SCPA's format_pathways helper function to convert msigdbr output
pathways <- format_pathways(hallmark_df)

print(paste("Formatted", length(pathways), "Hallmark pathways for SCPA"))
print("First 10 genes in first pathway:")
print(head(pathways[[1]], 10))

# Check overlap between Seurat genes and pathway genes
all_pathway_genes <- unique(unlist(pathways))
overlap <- intersect(seurat_genes, all_pathway_genes)
print(paste("\nGene overlap check:"))
print(paste("  Genes in pathways:", length(all_pathway_genes)))
print(paste("  Genes in Seurat object:", length(seurat_genes)))
print(paste("  Overlapping genes:", length(overlap)))
print(paste("  Overlap percentage:", round(length(overlap) / length(seurat_genes) * 100, 2), "%"))

# Run SCPA comparison  
print("\n=== Running SCPA comparison ===")
untreated_vs_vem <- compare_pathways(samples = list(untreated, vem), 
                                     pathways = pathways, 
                                     parallel = TRUE, 
                                     cores = 4)

print("\n=== SCPA Analysis Complete ===")
print(paste("Results dimensions:", paste(dim(untreated_vs_vem), collapse = " x ")))
print("\nTop 10 pathways by qval:")
print(head(untreated_vs_vem[order(untreated_vs_vem$qval), ], 10))

# Save results
write.csv(untreated_vs_vem, "SCPA_untreated_vs_vem.csv", row.names = FALSE)
print("\n✓ Saved results to SCPA_untreated_vs_vem.csv")

# =============================================================================
# VISUALIZE SCPA RESULTS
# =============================================================================

print("\n=== Creating Visualizations ===")

# 1. Pathway enrichment plot (bar plot of top pathways)
top_n <- 20
top_pathways <- untreated_vs_vem %>%
  arrange(qval) %>%
  head(top_n) %>%
  mutate(Pathway_clean = gsub("HALLMARK_", "", Pathway),
         Pathway_clean = gsub("_", " ", Pathway_clean))

p1 <- ggplot(top_pathways, aes(x = reorder(Pathway_clean, -FC), y = FC, fill = qval < 0.05)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#95A5A6"),
                    labels = c("TRUE" = "Significant (q < 0.05)", "FALSE" = "Not significant"),
                    name = "Significance") +
  labs(title = "SCPA: Untreated vs Vemurafenib",
       subtitle = paste("Top", top_n, "Pathways by q-value"),
       x = "Pathway",
       y = "Fold Change (FC)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11),
        axis.text.y = element_text(size = 9))

ggsave("SCPA_BarPlot_untreated_vs_vem.pdf", p1, width = 12, height = 10)
ggsave("SCPA_BarPlot_untreated_vs_vem.png", p1, width = 12, height = 10, dpi = 300)
print("✓ Saved bar plot")

# 2. Volcano plot (FC vs -log10(qval))
untreated_vs_vem <- untreated_vs_vem %>%
  mutate(Pathway_clean = gsub("HALLMARK_", "", Pathway),
         Pathway_clean = gsub("_", " ", Pathway_clean),
         neg_log_qval = -log10(qval + 1e-300),
         Significant = qval < 0.05)

p2 <- ggplot(untreated_vs_vem, aes(x = FC, y = neg_log_qval, 
                                    color = Significant, 
                                    size = abs(FC))) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  scale_color_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#95A5A6"),
                     labels = c("TRUE" = "Significant (q < 0.05)", "FALSE" = "Not significant")) +
  labs(title = "SCPA Volcano Plot: Untreated vs Vemurafenib",
       x = "Fold Change (FC)",
       y = "-Log10(q-value)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave("SCPA_VolcanoPlot_untreated_vs_vem.pdf", p2, width = 10, height = 8)
ggsave("SCPA_VolcanoPlot_untreated_vs_vem.png", p2, width = 10, height = 8, dpi = 300)
print("✓ Saved volcano plot")

# 3. Dot plot of significant pathways
sig_pathways <- untreated_vs_vem %>%
  filter(qval < 0.05) %>%
  arrange(desc(abs(FC))) %>%
  head(20)

if (nrow(sig_pathways) > 0) {
  p3 <- ggplot(sig_pathways, aes(x = FC, y = reorder(Pathway_clean, FC))) +
    geom_point(aes(size = -log10(qval), color = FC)) +
    scale_color_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C",
                         midpoint = 0, name = "Fold Change") +
    scale_size_continuous(name = "-Log10(q-value)") +
    labs(title = "Significant Pathways (q < 0.05)",
         subtitle = "Untreated vs Vemurafenib",
         x = "Fold Change",
         y = "Pathway") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14),
          axis.text.y = element_text(size = 9))
  
  ggsave("SCPA_DotPlot_Significant.pdf", p3, width = 12, height = 10)
  ggsave("SCPA_DotPlot_Significant.png", p3, width = 12, height = 10, dpi = 300)
  print("✓ Saved dot plot")
} else {
  print("No significant pathways found (q < 0.05)")
}

print("\n=== SCPA Visualization Complete ===")
print(paste("Total pathways analyzed:", nrow(untreated_vs_vem)))
print(paste("Significant pathways (q < 0.05):", sum(untreated_vs_vem$qval < 0.05)))

