# =============================================================================
# DIFFERENTIAL GENE EXPRESSION ANALYSIS
# =============================================================================
# Comprehensive DEG analysis for melanoma drug resistance
# Identifies genes driving resistance across treatments and cell types
# =============================================================================

library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

# Load reproducibility configuration
source("/projects/vanaja_lab/satya/SCPA/Reproduce.R")
set_analysis_seed("DEG_analysis")

cat("\n=============================================================================\n")
cat("DIFFERENTIAL GENE EXPRESSION ANALYSIS\n")
cat("=============================================================================\n\n")

# Check if GSE164897 exists
if (!exists("GSE164897")) {
  cat("Loading GSE164897 dataset...\n")
  source("GSE164897.R")
}

# Create output directory
dir.create("DEG_Results", showWarnings = FALSE)

# =============================================================================
# 1. TREATMENT COMPARISONS - OVERALL
# =============================================================================

cat("\n=== 1. TREATMENT vs UNTREATED COMPARISONS ===\n")

treatments <- c("Vemurafenib", "vem_cob", "vem_tram")
treatment_degs <- list()

for (trt in treatments) {
  cat(paste("\nAnalyzing:", trt, "vs untreated\n"))
  
  # Set identity to treatment
  Idents(GSE164897) <- "treatment"
  
  # Find markers
  degs <- FindMarkers(
    GSE164897,
    ident.1 = trt,
    ident.2 = "untreated",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    verbose = FALSE
  )
  
  # Add gene names and statistics
  degs <- degs %>%
    rownames_to_column("gene") %>%
    mutate(
      treatment = trt,
      significant = p_val_adj < 0.05 & abs(avg_log2FC) > 0.5,
      direction = ifelse(avg_log2FC > 0, "Up", "Down")
    ) %>%
    arrange(p_val_adj)
  
  treatment_degs[[trt]] <- degs
  
  # Summary
  n_up <- sum(degs$avg_log2FC > 0.5 & degs$p_val_adj < 0.05)
  n_down <- sum(degs$avg_log2FC < -0.5 & degs$p_val_adj < 0.05)
  
  cat(paste("  Upregulated genes:", n_up, "\n"))
  cat(paste("  Downregulated genes:", n_down, "\n"))
  cat(paste("  Total significant DEGs:", n_up + n_down, "\n"))
  
  # Save results
  write.csv(degs, 
            file.path("DEG_Results", paste0("DEGs_", trt, "_vs_untreated.csv")),
            row.names = FALSE)
}

# =============================================================================
# 2. CELL TYPE-SPECIFIC DEGs
# =============================================================================

cat("\n=== 2. CELL TYPE-SPECIFIC DEGs ===\n")

cell_types <- unique(GSE164897$celltype)
celltype_degs <- list()

for (ct in cell_types) {
  cat(paste("\nAnalyzing cell type:", ct, "\n"))
  
  # Subset to this cell type
  ct_subset <- subset(GSE164897, celltype == ct)
  
  if (ncol(ct_subset) < 50) {
    cat(paste("  Skipping", ct, "- too few cells\n"))
    next
  }
  
  Idents(ct_subset) <- "treatment"
  
  ct_degs <- list()
  
  for (trt in treatments) {
    # Check if both groups have enough cells
    n_trt <- sum(ct_subset$treatment == trt)
    n_untreated <- sum(ct_subset$treatment == "untreated")
    
    if (n_trt < 10 || n_untreated < 10) {
      cat(paste("  Skipping", trt, "- insufficient cells\n"))
      next
    }
    
    tryCatch({
      degs <- FindMarkers(
        ct_subset,
        ident.1 = trt,
        ident.2 = "untreated",
        test.use = "wilcox",
        logfc.threshold = 0.25,
        min.pct = 0.1,
        verbose = FALSE
      )
      
      degs <- degs %>%
        rownames_to_column("gene") %>%
        mutate(
          celltype = ct,
          treatment = trt,
          significant = p_val_adj < 0.05 & abs(avg_log2FC) > 0.5
        )
      
      ct_degs[[trt]] <- degs
      
      n_sig <- sum(degs$significant)
      cat(paste("  ", trt, ":", n_sig, "significant DEGs\n"))
      
    }, error = function(e) {
      cat(paste("  Error in", trt, ":", e$message, "\n"))
    })
  }
  
  if (length(ct_degs) > 0) {
    celltype_degs[[ct]] <- bind_rows(ct_degs)
    
    # Save cell type-specific results
    write.csv(celltype_degs[[ct]],
              file.path("DEG_Results", paste0("DEGs_", gsub(" ", "_", ct), "_all_treatments.csv")),
              row.names = FALSE)
  }
}

# =============================================================================
# 3. RESISTANCE STAGE DEGs (Early vs Late Pseudotime)
# =============================================================================

cat("\n=== 3. RESISTANCE STAGE DEGs (Early vs Late) ===\n")

if ("pseudotime_clinical" %in% colnames(GSE164897@meta.data)) {
  
  Idents(GSE164897) <- "pseudotime_clinical"
  
  stage_degs <- FindMarkers(
    GSE164897,
    ident.1 = "Resistant",
    ident.2 = "Sensitive",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    verbose = FALSE
  )
  
  stage_degs <- stage_degs %>%
    rownames_to_column("gene") %>%
    mutate(
      significant = p_val_adj < 0.05 & abs(avg_log2FC) > 0.5,
      direction = ifelse(avg_log2FC > 0, "Resistant", "Sensitive")
    ) %>%
    arrange(p_val_adj)
  
  n_resistant <- sum(stage_degs$avg_log2FC > 0.5 & stage_degs$p_val_adj < 0.05)
  n_sensitive <- sum(stage_degs$avg_log2FC < -0.5 & stage_degs$p_val_adj < 0.05)
  
  cat(paste("  Resistant-enriched genes:", n_resistant, "\n"))
  cat(paste("  Sensitive-enriched genes:", n_sensitive, "\n"))
  
  write.csv(stage_degs,
            file.path("DEG_Results", "DEGs_Resistant_vs_Sensitive.csv"),
            row.names = FALSE)
  
} else {
  cat("  Pseudotime stages not found. Run Psudotime.R first.\n")
  stage_degs <- NULL
}

# =============================================================================
# 4. CELL TYPE MARKERS
# =============================================================================

cat("\n=== 4. CELL TYPE MARKER GENES ===\n")

Idents(GSE164897) <- "celltype"

celltype_markers <- FindAllMarkers(
  GSE164897,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.5,
  test.use = "wilcox",
  verbose = FALSE
)

celltype_markers <- celltype_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE)

cat(paste("  Found markers for", length(unique(celltype_markers$cluster)), "cell types\n"))

write.csv(celltype_markers,
          file.path("DEG_Results", "CellType_Markers.csv"),
          row.names = FALSE)

# =============================================================================
# 5. VISUALIZATIONS
# =============================================================================

cat("\n=== 5. CREATING VISUALIZATIONS ===\n")

# --- 5.1 Volcano Plots for Each Treatment ---
cat("  Creating volcano plots...\n")

for (trt in names(treatment_degs)) {
  degs <- treatment_degs[[trt]]
  
  # Label top genes
  top_genes <- degs %>%
    filter(significant) %>%
    arrange(desc(abs(avg_log2FC))) %>%
    head(20)
  
  p <- ggplot(degs, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    geom_text_repel(data = top_genes, aes(label = gene), 
                    size = 3, max.overlaps = 20, box.padding = 0.5) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#E74C3C"),
                      labels = c("Not Sig", "Significant")) +
    labs(title = paste("DEGs:", trt, "vs Untreated"),
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value",
         color = NULL) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          panel.border = element_rect(fill = NA, color = "black"),
          legend.position = "bottom")
  
  ggsave(file.path("DEG_Results", paste0("Volcano_", trt, "_vs_untreated.pdf")),
         p, width = 10, height = 8)
  
  print(p)
}

# --- 5.2 Heatmap of Top DEGs ---
cat("  Creating DEG heatmap...\n")

# Get top DEGs from each treatment
top_degs_per_treatment <- lapply(treatment_degs, function(df) {
  df %>%
    filter(significant) %>%
    arrange(desc(abs(avg_log2FC))) %>%
    head(30) %>%
    pull(gene)
})

all_top_genes <- unique(unlist(top_degs_per_treatment))

if (length(all_top_genes) > 0) {
  # Get expression matrix
  expr_matrix <- GetAssayData(GSE164897, layer = "data")[all_top_genes, ]
  
  # Average by treatment
  avg_expr <- sapply(c("untreated", treatments), function(trt) {
    cells <- colnames(GSE164897)[GSE164897$treatment == trt]
    rowMeans(expr_matrix[, cells, drop = FALSE])
  })
  
  # Scale for visualization
  avg_expr_scaled <- t(scale(t(avg_expr)))
  
  # Remove rows with NA/Inf values
  valid_rows <- apply(avg_expr_scaled, 1, function(x) all(is.finite(x)))
  avg_expr_scaled <- avg_expr_scaled[valid_rows, , drop = FALSE]
  
  if (nrow(avg_expr_scaled) > 0) {
    # Create heatmap
    pheatmap(avg_expr_scaled,
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             color = colorRampPalette(c("blue", "white", "red"))(100),
             breaks = seq(-2, 2, length.out = 101),
             main = "Top DEGs Across Treatments",
             fontsize = 10,
             fontsize_row = 6,
             angle_col = 45,
             filename = file.path("DEG_Results", "Heatmap_Top_DEGs.pdf"),
             width = 8,
             height = 12)
  }
}

# --- 5.3 Resistance Signature Heatmap ---
if (!is.null(stage_degs)) {
  cat("  Creating resistance signature heatmap...\n")
  
  resistance_genes <- stage_degs %>%
    filter(significant, avg_log2FC > 0.5) %>%
    arrange(desc(avg_log2FC)) %>%
    head(50) %>%
    pull(gene)
  
  if (length(resistance_genes) > 0) {
    expr_matrix <- GetAssayData(GSE164897, layer = "data")[resistance_genes, ]
    
    # Average by pseudotime stage
    avg_expr <- sapply(c("Sensitive", "Adaptive", "Resistant"), function(stage) {
      cells <- colnames(GSE164897)[GSE164897$pseudotime_clinical == stage]
      if (length(cells) > 0) {
        rowMeans(expr_matrix[, cells, drop = FALSE])
      } else {
        rep(NA, length(resistance_genes))
      }
    })
    
    avg_expr_scaled <- t(scale(t(avg_expr)))
    
    # Remove rows with NA/Inf values
    valid_rows <- apply(avg_expr_scaled, 1, function(x) all(is.finite(x)))
    avg_expr_scaled <- avg_expr_scaled[valid_rows, , drop = FALSE]
    
    if (nrow(avg_expr_scaled) > 0) {
      pheatmap(avg_expr_scaled,
               cluster_rows = TRUE,
               cluster_cols = FALSE,
               color = colorRampPalette(c("blue", "white", "red"))(100),
               breaks = seq(-2, 2, length.out = 101),
               main = "Resistance Gene Signature Across Stages",
               fontsize = 10,
               fontsize_row = 6,
               filename = file.path("DEG_Results", "Heatmap_Resistance_Signature.pdf"),
               width = 6,
               height = 12)
    }
  }
}

# --- 5.4 Cell Type Marker Heatmap ---
cat("  Creating cell type marker heatmap...\n")

top_markers <- celltype_markers %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC) %>%
  pull(gene) %>%
  unique()

if (length(top_markers) > 0) {
  expr_matrix <- GetAssayData(GSE164897, layer = "data")[top_markers, ]
  
  # Average by cell type
  cell_types_ordered <- sort(unique(GSE164897$celltype))
  avg_expr <- sapply(cell_types_ordered, function(ct) {
    cells <- colnames(GSE164897)[GSE164897$celltype == ct]
    rowMeans(expr_matrix[, cells, drop = FALSE])
  })
  
  avg_expr_scaled <- t(scale(t(avg_expr)))
  
  # Remove rows with NA/Inf values
  valid_rows <- apply(avg_expr_scaled, 1, function(x) all(is.finite(x)))
  avg_expr_scaled <- avg_expr_scaled[valid_rows, , drop = FALSE]
  
  if (nrow(avg_expr_scaled) > 0) {
    pheatmap(avg_expr_scaled,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             color = colorRampPalette(c("blue", "white", "red"))(100),
             breaks = seq(-2, 2, length.out = 101),
             main = "Top Cell Type Markers",
             fontsize = 10,
             fontsize_row = 6,
             angle_col = 45,
             filename = file.path("DEG_Results", "Heatmap_CellType_Markers.pdf"),
             width = 10,
             height = 14)
  }
}

# =============================================================================
# 6. SUMMARY STATISTICS
# =============================================================================

cat("\n=== 6. SUMMARY STATISTICS ===\n")

# Create summary table
summary_df <- data.frame(
  Treatment = character(),
  Total_DEGs = integer(),
  Upregulated = integer(),
  Downregulated = integer(),
  stringsAsFactors = FALSE
)

for (trt in names(treatment_degs)) {
  degs <- treatment_degs[[trt]]
  total <- sum(degs$significant)
  up <- sum(degs$avg_log2FC > 0.5 & degs$p_val_adj < 0.05)
  down <- sum(degs$avg_log2FC < -0.5 & degs$p_val_adj < 0.05)
  
  summary_df <- rbind(summary_df, data.frame(
    Treatment = trt,
    Total_DEGs = total,
    Upregulated = up,
    Downregulated = down
  ))
}

print(summary_df)

write.csv(summary_df,
          file.path("DEG_Results", "DEG_Summary_Statistics.csv"),
          row.names = FALSE)

# =============================================================================
# 7. KNOWN RESISTANCE GENES
# =============================================================================

cat("\n=== 7. KNOWN RESISTANCE GENE ANALYSIS ===\n")

# Known melanoma resistance genes
known_resistance_genes <- c(
  "AXL", "NGFR", "EGFR", "PDGFRB", "MET", "IGF1R",
  "SOX10", "MITF", "SOX9", "FOXD3",
  "VIM", "TWIST1", "ZEB1", "SNAI2", "CDH2",
  "WNT5A", "JUN", "FOS", "NFKB1", "STAT3"
)

# Check which are present in our data
present_genes <- intersect(known_resistance_genes, rownames(GSE164897))

cat(paste("  Found", length(present_genes), "of", length(known_resistance_genes), 
          "known resistance genes\n"))

# Extract their expression across treatments
if (length(present_genes) > 0) {
  known_gene_expr <- GetAssayData(GSE164897, layer = "data")[present_genes, ]
  
  avg_expr <- sapply(c("untreated", treatments), function(trt) {
    cells <- colnames(GSE164897)[GSE164897$treatment == trt]
    rowMeans(known_gene_expr[, cells, drop = FALSE])
  })
  
  # Calculate fold changes
  fc_matrix <- log2(avg_expr[, treatments] / avg_expr[, "untreated"])
  
  # Remove rows with NA/Inf values
  valid_rows <- apply(fc_matrix, 1, function(x) all(is.finite(x)))
  fc_matrix <- fc_matrix[valid_rows, , drop = FALSE]
  
  if (nrow(fc_matrix) > 0) {
    # Create heatmap
    pheatmap(fc_matrix,
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             color = colorRampPalette(c("blue", "white", "red"))(100),
             breaks = seq(-2, 2, length.out = 101),
             main = "Known Resistance Genes - Log2 FC vs Untreated",
             fontsize = 10,
             angle_col = 45,
             filename = file.path("DEG_Results", "Heatmap_Known_Resistance_Genes.pdf"),
             width = 6,
             height = 8)
  }
  
  # Save expression data
  write.csv(avg_expr,
            file.path("DEG_Results", "Known_Resistance_Genes_Expression.csv"))
}

# =============================================================================
# 8. SAVE WORKSPACE
# =============================================================================

cat("\n=== 8. SAVING RESULTS ===\n")

# Save all DEG results
save(treatment_degs, celltype_degs, stage_degs, celltype_markers,
     file = file.path("DEG_Results", "All_DEG_Results.RData"))

# Save session info
save_session_info(file.path("DEG_Results", "session_info_DEG.txt"))

cat("\n=============================================================================\n")
cat("✓ DIFFERENTIAL GENE EXPRESSION ANALYSIS COMPLETE\n")
cat("=============================================================================\n")
cat(paste("Results saved in:", file.path(getwd(), "DEG_Results"), "\n"))
cat("\nKey outputs:\n")
cat("  - DEG tables for each treatment comparison\n")
cat("  - Cell type-specific DEGs\n")
cat("  - Resistance stage DEGs\n")
cat("  - Cell type marker genes\n")
cat("  - Volcano plots\n")
cat("  - Expression heatmaps\n")
cat("  - Known resistance gene analysis\n")
cat("=============================================================================\n\n")
