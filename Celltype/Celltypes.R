# Melanoma Cell State Annotation and Cluster Naming
# This script identifies melanoma cell states based on specific gene signatures (C1-C7)
# and renames Seurat clusters accordingly.

library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tibble)

# Load reproducibility configuration
if (file.exists("/projects/vanaja_lab/satya/SCPA/Reproduce.R")) {
  source("/projects/vanaja_lab/satya/SCPA/Reproduce.R")
  set_analysis_seed("melanoma_states")
}

# Check if GSE164897 object exists
if (!exists("GSE164897")) {
  if (file.exists("GSE164897.R")) {
    source("GSE164897.R")
  } else {
    stop("ERROR: GSE164897 Seurat object not found and GSE164897.R missing.")
  }
}

message("\n=== Melanoma Cell State Annotation based on C1-C7 Signatures ===")

# =============================================================================
# 1. DEFINE GENE SIGNATURES FROM SCREENSHOT
# =============================================================================

# Gene lists extracted from user-provided table
# Note: Aliases used where common (e.g., PTGS2 for COX2, POU5F1 for OCT4)
c_signatures <- list(
  Neural_Crest_like = c("MAJIN", "NME9", "CD74", "ERBB3", "RAMP1", "HLA-DRA", "SOX10", "KRT23"),
  
  High_Cycling_G1S = c("UBE2S", "KRT81", "KRT17", "HMGN2", "ANKRD1", "FOSL1", 
                       "KRT18", "RANBP1", "STMN1", "CENPU", "CDK1", "CDK4", "MYC"),
  
  Slow_Cycling_Stromal = c("CDKN1A", "KDM5B", "IL1B", "H1-2", "GSN", 
                           "FAM3C", "MTOR", "MDM2"),
  
  # C4 refers to C2 in table, but represents G2M. We include C2 genes plus standard G2M markers to differentiate
  High_Cycling_G2M = c("UBE2S", "CDK1", "TOP2A", "MKI67", "CCNB1", "CENPF", "NUSAP1", "UBE2C"),
  
  Translation = c("S100A6", "RPL7A", "RPL31", "RACK1", "EEF1A1"),
  
  Stressed_Pluripotency = c("UTS2", "PTGS2", "NTS", "MMP1", "EGR1", 
                            "SUCNR1", "PCSK6", "POU5F1", "SETD1A", "KLF4"), 
  # COX2->PTGS2, OCT4->POU5F1
  
  MAPK_Reactivating = c("EPO", "MARCKS", "TSPAN8", "GUCY1B1", "PDE6B", "DDX10", 
                        "BRAF", "WDR63", "SPANXD")
)

# Filter markers to those present in dataset
c_signatures_filtered <- lapply(c_signatures, function(x) intersect(x, rownames(GSE164897)))

message("Signature gene availability:")
for (sig in names(c_signatures_filtered)) {
  pct <- round(length(c_signatures_filtered[[sig]]) / length(c_signatures[[sig]]) * 100, 1)
  message(paste0("  ", sig, ": ", length(c_signatures_filtered[[sig]]), "/", length(c_signatures[[sig]]), " (", pct, "%)"))
}

# =============================================================================
# 2. FIND MOST EXPRESSED GENES (MARKERS) PER CLUSTER
# =============================================================================

message("\n=== Finding most expressed genes for each Seurat cluster ===")

# Ensure we are using seurat_clusters
Idents(GSE164897) <- "seurat_clusters"

# Find top markers for every cluster
all_markers <- FindAllMarkers(
  GSE164897,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.5,
  test.use = "wilcox",
  verbose = FALSE
)

# Print top 5 markers per cluster for verification
message("Top 5 markers per cluster:")
top_markers_df <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  select(cluster, gene, avg_log2FC)

print(as.data.frame(top_markers_df))

# =============================================================================
# 3. SCORE CLUSTERS AND ASSIGN C1-C7 NAMES
# =============================================================================

message("\n=== Scoring and Naming Clusters ===")

# Calculate Module Scores for each signature
# We clean up old score columns first to avoid duplication issues
old_scores <- grep("^C[1-7]_", colnames(GSE164897@meta.data))
if (length(old_scores) > 0) {
  GSE164897@meta.data <- GSE164897@meta.data[, -old_scores]
}

# Add scores
for (sig_name in names(c_signatures_filtered)) {
  if (length(c_signatures_filtered[[sig_name]]) > 0) {
    GSE164897 <- AddModuleScore(
      GSE164897,
      features = list(c_signatures_filtered[[sig_name]]),
      name = sig_name,
      seed = 42
    )
    
    # Fix Seurat's naming convention (it adds '1')
    current_col <- paste0(sig_name, "1")
    if (current_col %in% colnames(GSE164897@meta.data)) {
      GSE164897@meta.data[[sig_name]] <- GSE164897@meta.data[[current_col]]
      GSE164897@meta.data[[current_col]] <- NULL
    }
  }
}

# Determine the best matching signature for each Cluster
cluster_ids <- unique(GSE164897$seurat_clusters)
mapping_df <- data.frame(Cluster = character(), best_match = character(), score = numeric(), stringsAsFactors = FALSE)

signature_names <- names(c_signatures_filtered)
# Exclude "Translation" from the "best match" competition because it is often 
# highly expressed in all cells (ribosomal genes) and masks specific identities.
signature_names <- setdiff(signature_names, "Translation")

for (cl in cluster_ids) {
  # Get cells in this cluster
  cells_in_cluster <- WhichCells(GSE164897, idents = cl)
  
  # Calculate average score for each signature in this cluster
  avg_scores <- numeric(length(signature_names))
  names(avg_scores) <- signature_names
  
  for (sig in signature_names) {
    avg_scores[sig] <- mean(GSE164897@meta.data[cells_in_cluster, sig], na.rm = TRUE)
  }
  
  # Find max
  best_sig <- names(which.max(avg_scores))
  max_val <- max(avg_scores)
  
  mapping_df <- rbind(mapping_df, data.frame(Cluster = cl, best_match = best_sig, score = max_val))
}

# Simplify names for assignment
mapping_df$short_name <- gsub("_", " ", mapping_df$best_match)
mapping_df$short_name <- gsub("like", "like", mapping_df$short_name) # aesthetic cleanup

message("Cluster Annotation Results:")
print(mapping_df)

# Rename clusters in the Object
new_ids <- setNames(paste0(mapping_df$Cluster, ": ", mapping_df$short_name), mapping_df$Cluster)
GSE164897 <- RenameIdents(GSE164897, new_ids)

# Store this in metadata as well
GSE164897$C_State_Annotation <- Idents(GSE164897)

message("\nRenamed identities in Seurat object to match C1-C7 states.")

# =============================================================================
# 4. VISUALIZATION
# =============================================================================

# Heatmap of Signature Scores by Original Cluster
start_cluster_col <- which(colnames(GSE164897@meta.data) == "seurat_clusters")
score_cols_indices <- match(signature_names, colnames(GSE164897@meta.data))

# Valid score columns
valid_sigs <- signature_names[!is.na(score_cols_indices)]

if (length(valid_sigs) > 0) {
  avg_sig_scores <- GSE164897@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(across(all_of(valid_sigs), mean)) %>%
    column_to_rownames("seurat_clusters")
  
  # Normalize row-wise for heatmap (z-score)
  # avg_sig_scores_scaled <- t(scale(t(avg_sig_scores))) 
  
  pheatmap(avg_sig_scores,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           main = "Average C-State Scores per Cluster",
           color = colorRampPalette(c("white", "red"))(100),
           display_numbers = TRUE,
           number_format = "%.2f",
           fontsize = 10,
           angle_col = 45,
           filename = "Melanoma_C_State_Assignment_Heatmap.pdf",
           width = 8, height = 6)
  
  message("Saved assignment heatmap to Melanoma_C_State_Assignment_Heatmap.pdf")
}

# UMAP with new labels
p_umap <- DimPlot(GSE164897, reduction = "umap", label = TRUE, repel = TRUE, label.size = 3) +
  ggtitle("Clusters Annotated by C1-C7 Signatures") +
  theme(legend.position = "bottom", legend.direction = "vertical")

print(p_umap)
ggsave("Annotated_Clusters_UMAP.pdf", plot = p_umap, width = 10, height = 8)

# DotPlot of Original Signatures
# To verify, let's plot the top marker from each signature across the new clusters
top_genes_per_sig <- sapply(c_signatures_filtered, function(x) x[1]) # Pick first gene as proxy
top_genes_per_sig <- top_genes_per_sig[!is.na(top_genes_per_sig)]

p_dot <- DotPlot(GSE164897, features = unique(unlist(c_signatures_filtered))) +
  RotatedAxis() +
  ggtitle("C-State Marker Expression") +
  theme(axis.text.x = element_text(size = 8))

ggsave("C_State_Markers_DotPlot.pdf", plot = p_dot, width = 14, height = 8)

message("\n=== Analysis Complete ===\n")
