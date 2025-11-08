# =============================================================================
# GENE SET ENRICHMENT ANALYSIS (GSEA)
# Perform GSEA on differentially expressed genes between treatment groups
# Seurat object: GSE164897
# =============================================================================

library(Seurat)
library(fgsea)
library(ggplot2)
library(dplyr)
library(msigdbr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(DT)
library(pheatmap)

print("═══════════════════════════════════════════════════════════════════════════")
print("   GENE SET ENRICHMENT ANALYSIS (GSEA)")
print("═══════════════════════════════════════════════════════════════════════════")

# Check if GSE164897 object exists
if (!exists("GSE164897")) {
  stop("ERROR: GSE164897 Seurat object not found. Please load it first.")
}

# =============================================================================
# GENE ID CONVERSION: Ensembl IDs to Gene Symbols (if needed)
# =============================================================================

print("\n═══════════════════════════════════════════════════════════════════════════")
print("   GENE ID CONVERSION: Ensembl IDs → Gene Symbols")
print("═══════════════════════════════════════════════════════════════════════════")

current_genes <- rownames(GSE164897)
sample_genes <- head(current_genes, 100)
is_ensembl <- grepl("^ENSG", sample_genes, ignore.case = TRUE)
ensembl_ratio <- sum(is_ensembl) / length(is_ensembl)

print(paste("Total genes:", length(current_genes)))
print(paste("Ensembl ID pattern detected in", round(ensembl_ratio * 100, 1), "% of sample genes"))

# Create mapping from Ensembl IDs to Gene Symbols
if (ensembl_ratio > 0.5) {
  print("Detected Ensembl IDs - creating mapping to gene symbols...")
  
  # Strip version suffix if present
  ens_ids <- sub("\\..*$", "", current_genes)
  
  # Map Ensembl IDs to Gene Symbols
  tryCatch({
    gene_mapping <- AnnotationDbi::select(org.Hs.eg.db,
                                         keys = unique(ens_ids),
                                         keytype = "ENSEMBL",
                                         columns = "SYMBOL")
    
    # Remove NAs and duplicates
    gene_mapping <- gene_mapping[!is.na(gene_mapping$SYMBOL), ]
    gene_mapping <- gene_mapping[!duplicated(gene_mapping$ENSEMBL), ]
    
    # Create named vector
    ens_to_symbol <- setNames(gene_mapping$SYMBOL, gene_mapping$ENSEMBL)
    
    print(paste("Successfully mapped", length(ens_to_symbol), "Ensembl IDs to gene symbols"))
    
  }, error = function(e) {
    print(paste("Error in gene ID mapping:", e$message))
    stop("Cannot proceed without gene symbol mapping for GSEA")
  })
  
} else {
  print("Rownames appear to be gene symbols already - no conversion needed")
  ens_to_symbol <- setNames(current_genes, current_genes)
}

# =============================================================================
# LOAD HALLMARK GENE SETS
# =============================================================================

print("\n═══════════════════════════════════════════════════════════════════════════")
print("   LOADING HALLMARK GENE SETS")
print("═══════════════════════════════════════════════════════════════════════════")

# Fetch ALL hallmark gene sets from MSigDB
print("Fetching ALL hallmark gene sets from MSigDB...")
tryCatch({
  msigdb_hallmark <- tryCatch({
    msigdbr(species = "Homo sapiens", collection = "H")
  }, error = function(e) {
    msigdbr(species = "Homo sapiens", category = "H")
  })
  
  # Get all unique hallmark gene set names
  all_hallmark_names <- unique(msigdb_hallmark$gs_name)
  print(paste("Found", length(all_hallmark_names), "hallmark gene sets"))
  
  # Create gene set list for fgsea (named list of character vectors)
  pathways <- split(msigdb_hallmark$gene_symbol, msigdb_hallmark$gs_name)
  
  print(paste("\nSuccessfully loaded", length(pathways), "hallmark pathways"))
  print(paste("Pathway sizes range from", min(sapply(pathways, length)), 
              "to", max(sapply(pathways, length)), "genes"))
  
}, error = function(e) {
  print(paste("Error: Failed to fetch hallmark gene sets:", e$message))
  stop("Cannot proceed without hallmark gene sets. Please check your internet connection and msigdbr package.")
})

# =============================================================================
# DIFFERENTIAL EXPRESSION ANALYSIS
# =============================================================================

print("\n═══════════════════════════════════════════════════════════════════════════")
print("   DIFFERENTIAL EXPRESSION ANALYSIS")
print("═══════════════════════════════════════════════════════════════════════════")

if (!"treatment" %in% colnames(GSE164897@meta.data)) {
  stop("ERROR: Treatment column not found in Seurat object metadata.")
}

# Set identity to treatment
Idents(GSE164897) <- "treatment"

# Get treatment groups
treatments <- unique(GSE164897$treatment)
treatments <- treatments[!is.na(treatments)]
print(paste("Treatment groups:", paste(treatments, collapse = ", ")))

# Define comparisons: (ident.1 vs ident.2)
comparisons <- list(
  list(ident1 = "untreated", ident2 = "Vemurafenib", name = "Untreated_vs_Vemurafenib"),
  list(ident1 = "Vemurafenib", ident2 = "vem_tram", name = "Vemurafenib_vs_vem_tram"),
  list(ident1 = "vem_tram", ident2 = "vem_cob", name = "vem_tram_vs_vem_cob"),
  list(ident1 = "untreated", ident2 = "vem_cob", name = "Untreated_vs_vem_cob"),
  list(ident1 = "untreated", ident2 = "vem_tram", name = "Untreated_vs_vem_tram"),
  list(ident1 = "Vemurafenib", ident2 = "vem_cob", name = "Vemurafenib_vs_vem_cob")
)

# Store GSEA results
gsea_results <- list()

# =============================================================================
# PERFORM GSEA FOR EACH COMPARISON
# =============================================================================

print("\n═══════════════════════════════════════════════════════════════════════════")
print("   RUNNING GSEA")
print("═══════════════════════════════════════════════════════════════════════════")

for (comp in comparisons) {
  
  print(paste("\n", paste(rep("=", 60), collapse = "")))
  print(paste("Comparing", comp$ident1, "vs", comp$ident2))
  print(paste(paste(rep("=", 60), collapse = "")))
  
  # Find differentially expressed genes
  tryCatch({
    print("Finding differentially expressed genes...")
    markers <- FindMarkers(GSE164897,
                          ident.1 = comp$ident1,
                          ident.2 = comp$ident2,
                          test.use = "wilcox",
                          logfc.threshold = 0,
                          min.pct = 0.1,
                          verbose = FALSE)
    
    # Add gene names
    markers$gene_ensembl <- rownames(markers)
    
    # Convert Ensembl IDs to gene symbols if needed
    if (ensembl_ratio > 0.5) {
      marker_ens_no_version <- sub("\\..*$", "", markers$gene_ensembl)
      markers$gene <- ens_to_symbol[marker_ens_no_version]
      # Fill missing symbols with Ensembl IDs
      markers$gene[is.na(markers$gene)] <- markers$gene_ensembl[is.na(markers$gene)]
    } else {
      markers$gene <- markers$gene_ensembl
    }
    
    # Handle NA values
    markers$p_val_adj[is.na(markers$p_val_adj)] <- 1
    markers$avg_log2FC[is.na(markers$avg_log2FC)] <- 0
    
    # Remove genes without valid symbols
    markers <- markers[!is.na(markers$gene) & markers$gene != "", ]
    
    if (nrow(markers) == 0) {
      warning(paste("No valid genes found for", comp$name, "- skipping"))
      next
    }
    
    # Prepare ranked gene list for GSEA
    # Rank by -log10(p-value) * sign(log2FC) or by log2FC
    # Using a combined metric: signed log p-value
    markers$rank_score <- sign(markers$avg_log2FC) * (-log10(markers$p_val_adj + 1e-300))
    
    # Create named vector: gene symbol -> rank score
    ranks <- setNames(markers$rank_score, markers$gene)
    ranks <- ranks[!is.na(ranks) & is.finite(ranks)]
    
    print(paste("Ranked", length(ranks), "genes for GSEA"))
    
    # Filter pathways to only include genes present in ranked list
    pathways_filtered <- lapply(pathways, function(pathway) {
      pathway[pathway %in% names(ranks)]
    })
    
    # Remove empty pathways
    pathways_filtered <- pathways_filtered[sapply(pathways_filtered, length) >= 10]
    
    print(paste("Using", length(pathways_filtered), "pathways with at least 10 genes"))
    
    # Run GSEA
    print("Running GSEA (this may take a few minutes)...")
    set.seed(42)
    fgseaRes <- fgsea(pathways = pathways_filtered,
                      stats = ranks,
                      minSize = 10,
                      maxSize = 500,
                      nperm = 10000)
    
    # Sort by normalized enrichment score (NES)
    fgseaRes <- fgseaRes[order(fgseaRes$NES, decreasing = TRUE), ]
    
    print(paste("GSEA completed. Found", nrow(fgseaRes), "pathways"))
    
    # Print top enriched pathways
    print("\nTop 10 enriched pathways (positive NES):")
    top_enriched <- head(fgseaRes[fgseaRes$NES > 0, ], 10)
    if (nrow(top_enriched) > 0) {
      print(top_enriched[, c("pathway", "NES", "pval", "padj")])
    }
    
    print("\nTop 10 depleted pathways (negative NES):")
    top_depleted <- head(fgseaRes[fgseaRes$NES < 0, ], 10)
    if (nrow(top_depleted) > 0) {
      print(top_depleted[, c("pathway", "NES", "pval", "padj")])
    }
    
    # Add comparison name
    fgseaRes$comparison <- comp$name
    
    # Store results
    gsea_results[[comp$name]] <- fgseaRes
    
    # Save individual results
    write.csv(fgseaRes, paste0("GSEA_", comp$name, ".csv"), row.names = FALSE)
    print(paste("✓ Saved GSEA results for", comp$name))
    
  }, error = function(e) {
    print(paste("Error in", comp$name, "comparison:", e$message))
  })
}

# =============================================================================
# COMBINE AND SUMMARIZE RESULTS
# =============================================================================

print("\n═══════════════════════════════════════════════════════════════════════════")
print("   SUMMARIZING GSEA RESULTS")
print("═══════════════════════════════════════════════════════════════════════════")

if (length(gsea_results) > 0) {
  # Combine all results
  all_gsea_results <- do.call(rbind, gsea_results)
  
  # Save combined results
  write.csv(all_gsea_results, "GSEA_All_Comparisons.csv", row.names = FALSE)
  print("✓ Saved combined GSEA results")
  
  # Identify significantly enriched/depleted pathways (padj < 0.05)
  sig_pathways <- all_gsea_results[all_gsea_results$padj < 0.05, ]
  write.csv(sig_pathways, "GSEA_Significant_Pathways.csv", row.names = FALSE)
  print(paste("✓ Saved", nrow(sig_pathways), "significantly enriched/depleted pathways"))
  
  # Summary by comparison
  print("\nSignificant pathways per comparison:")
  for (comp_name in names(gsea_results)) {
    sig_count <- sum(gsea_results[[comp_name]]$padj < 0.05, na.rm = TRUE)
    print(paste("  ", comp_name, ":", sig_count, "significant pathways"))
  }
}

# =============================================================================
# VISUALIZE GSEA RESULTS
# =============================================================================

print("\n═══════════════════════════════════════════════════════════════════════════")
print("   VISUALIZING GSEA RESULTS")
print("═══════════════════════════════════════════════════════════════════════════")

for (comp in comparisons) {
  comp_name <- comp$name
  
  if (!comp_name %in% names(gsea_results)) {
    next
  }
  
  fgseaRes <- gsea_results[[comp_name]]
  
  # 1. Bar plot of top enriched/depleted pathways
  print(paste("Creating visualizations for", comp_name, "..."))
  
  # Get top pathways (by absolute NES)
  top_pathways <- fgseaRes[order(abs(fgseaRes$NES), decreasing = TRUE), ]
  top_pathways <- head(top_pathways, 20)
  
  # Clean pathway names
  top_pathways$pathway_clean <- gsub("HALLMARK_", "", top_pathways$pathway)
  top_pathways$pathway_clean <- gsub("_", " ", top_pathways$pathway_clean)
  
  # Bar plot
  p_bar <- ggplot(top_pathways, aes(x = reorder(pathway_clean, NES), y = NES, fill = NES > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#D62728", "FALSE" = "#2CA02C"),
                      labels = c("TRUE" = "Enriched", "FALSE" = "Depleted"),
                      name = "Direction") +
    labs(title = paste("GSEA Results:", comp$ident1, "vs", comp$ident2),
         subtitle = "Top 20 Pathways by |NES|",
         x = "Pathway",
         y = "Normalized Enrichment Score (NES)") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 12),
          axis.text.y = element_text(size = 9))
  
  ggsave(paste0("GSEA_BarPlot_", comp_name, ".pdf"), p_bar, width = 12, height = 10)
  ggsave(paste0("GSEA_BarPlot_", comp_name, ".png"), p_bar, width = 12, height = 10, dpi = 300)
  
  # 2. Scatter plot: NES vs -log10(padj)
  p_scatter <- ggplot(fgseaRes, aes(x = NES, y = -log10(padj + 1e-300), 
                                     color = padj < 0.05, size = size)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("TRUE" = "#D62728", "FALSE" = "#95A5A6"),
                       labels = c("TRUE" = "Significant (padj < 0.05)", "FALSE" = "Not significant"),
                       name = "Significance") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    labs(title = paste("GSEA Results:", comp$ident1, "vs", comp$ident2),
         x = "Normalized Enrichment Score (NES)",
         y = "-Log10(Adjusted P-value)") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14))
  
  ggsave(paste0("GSEA_ScatterPlot_", comp_name, ".pdf"), p_scatter, width = 10, height = 8)
  ggsave(paste0("GSEA_ScatterPlot_", comp_name, ".png"), p_scatter, width = 10, height = 8, dpi = 300)
  
  # 3. Enrichment plots for top pathways
  if (nrow(fgseaRes) > 0) {
    # Get the ranked gene list for this comparison
    markers <- FindMarkers(GSE164897,
                          ident.1 = comp$ident1,
                          ident.2 = comp$ident2,
                          test.use = "wilcox",
                          logfc.threshold = 0,
                          min.pct = 0.1,
                          verbose = FALSE)
    
    if (ensembl_ratio > 0.5) {
      marker_ens_no_version <- sub("\\..*$", "", rownames(markers))
      markers$gene <- ens_to_symbol[marker_ens_no_version]
      markers$gene[is.na(markers$gene)] <- rownames(markers)[is.na(markers$gene)]
    } else {
      markers$gene <- rownames(markers)
    }
    
    markers$rank_score <- sign(markers$avg_log2FC) * (-log10(markers$p_val_adj + 1e-300))
    ranks <- setNames(markers$rank_score, markers$gene)
    ranks <- ranks[!is.na(ranks) & is.finite(ranks)]
    
    # Select top 5 enriched and top 5 depleted pathways for enrichment plots
    top_enriched_pathways <- head(fgseaRes[fgseaRes$NES > 0 & fgseaRes$padj < 0.05, ], 5)
    top_depleted_pathways <- head(fgseaRes[fgseaRes$NES < 0 & fgseaRes$padj < 0.05, ], 5)
    
    # Create enrichment plots
    if (nrow(top_enriched_pathways) > 0) {
      pdf(paste0("GSEA_EnrichmentPlots_", comp_name, ".pdf"), width = 12, height = 10)
      for (i in 1:min(5, nrow(top_enriched_pathways))) {
        pathway_name <- top_enriched_pathways$pathway[i]
        plotEnrichment(pathways_filtered[[pathway_name]], ranks) +
          labs(title = paste("Enrichment Plot:", gsub("HALLMARK_", "", pathway_name)),
               subtitle = paste("NES =", round(top_enriched_pathways$NES[i], 3),
                              ", padj =", format(top_enriched_pathways$padj[i], scientific = TRUE)))
      }
      dev.off()
      print(paste("  ✓ Saved enrichment plots for", comp_name))
    }
  }
  
  print(paste("✓ Saved visualizations for", comp_name))
}

# =============================================================================
# HEATMAP OF PATHWAY ENRICHMENT ACROSS COMPARISONS
# =============================================================================

if (length(gsea_results) > 1) {
  print("\nCreating heatmap of pathway enrichment across comparisons...")
  
  # Create matrix: pathways x comparisons
  all_pathways <- unique(unlist(lapply(gsea_results, function(x) x$pathway)))
  
  nes_matrix <- matrix(NA, nrow = length(all_pathways), ncol = length(gsea_results))
  rownames(nes_matrix) <- all_pathways
  colnames(nes_matrix) <- names(gsea_results)
  
  for (comp_name in names(gsea_results)) {
    fgseaRes <- gsea_results[[comp_name]]
    nes_matrix[fgseaRes$pathway, comp_name] <- fgseaRes$NES
  }
  
  # Filter to pathways that are significant in at least one comparison
  sig_in_any <- apply(nes_matrix, 1, function(x) {
    any(!is.na(x) & abs(x) > 0)
  })
  nes_matrix_filtered <- nes_matrix[sig_in_any, , drop = FALSE]
  
  # Clean pathway names for display
  rownames(nes_matrix_filtered) <- gsub("HALLMARK_", "", rownames(nes_matrix_filtered))
  rownames(nes_matrix_filtered) <- gsub("_", " ", rownames(nes_matrix_filtered))
  
  if (nrow(nes_matrix_filtered) > 0) {
    pdf("GSEA_Heatmap_All_Comparisons.pdf", width = 10, height = max(8, nrow(nes_matrix_filtered) * 0.3))
    pheatmap(nes_matrix_filtered,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             color = colorRampPalette(c("blue", "white", "red"))(100),
             main = "GSEA NES Heatmap Across Comparisons",
             fontsize = 8,
             fontsize_row = 7,
             fontsize_col = 10,
             show_rownames = TRUE,
             show_colnames = TRUE,
             na_col = "grey90")
    dev.off()
    
    png("GSEA_Heatmap_All_Comparisons.png", width = 10, height = max(8, nrow(nes_matrix_filtered) * 0.3), 
        units = "in", res = 300)
    pheatmap(nes_matrix_filtered,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             color = colorRampPalette(c("blue", "white", "red"))(100),
             main = "GSEA NES Heatmap Across Comparisons",
             fontsize = 8,
             fontsize_row = 7,
             fontsize_col = 10,
             show_rownames = TRUE,
             show_colnames = TRUE,
             na_col = "grey90")
    dev.off()
    
    print("✓ Saved GSEA heatmap across comparisons")
  }
}

print("\n═══════════════════════════════════════════════════════════════════════════")
print("   GSEA ANALYSIS COMPLETED SUCCESSFULLY!")
print("═══════════════════════════════════════════════════════════════════════════")
print("\nGSEA results have been saved to CSV files and visualizations to PDF/PNG files.")

