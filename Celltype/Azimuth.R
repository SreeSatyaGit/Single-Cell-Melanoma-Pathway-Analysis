# Reference-Based Cell Type Transfer
# Using GSE115978 (reference) to annotate GSE164897 (query)
# This transfers cell type labels from the annotated GSE115978 dataset to GSE164897

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)

message("\n=== Reference-Based Cell Type Transfer ===")
message("Reference: GSE115978 (pre-annotated)")
message("Query: GSE164897 (to be annotated)")

# ============================================================================
# STEP 1: Check Required Objects
# ============================================================================

message("\n=== Checking Required Objects ===")

# Check if GSE115978 exists (reference)
if (!exists("GSE115978")) {
  message("GSE115978 not found. Loading reference dataset...")
  source("Melanoma_Annotated/GES115978.R")
  message("GSE115978 loaded successfully")
} else {
  message("✓ GSE115978 found (reference dataset)")
}

# Check if GSE164897 exists (query)
if (!exists("GSE164897")) {
  message("GSE164897 not found. Loading query dataset...")
  source("GSE164897.R")
  message("GSE164897 loaded successfully")
} else {
  message("✓ GSE164897 found (query dataset)")
}

# ============================================================================
# STEP 2: Identify Reference Cell Type Column
# ============================================================================

message("\n=== Identifying Reference Cell Type Column ===")

# Find potential cell type columns in GSE115978
ref_metadata_cols <- colnames(GSE115978@meta.data)
message("Available columns in GSE115978:")
print(ref_metadata_cols)

# Look for cell type columns (excluding our melanoma_state)
celltype_candidates <- grep("celltype|cell_type|CellType|annotation|cluster|malignant|cell.type", 
                            ref_metadata_cols, 
                            value = TRUE, 
                            ignore.case = TRUE)

# Exclude our added columns
celltype_candidates <- setdiff(celltype_candidates, 
                               c("melanoma_state", "melanoma_state_score", 
                                 grep("_score", celltype_candidates, value = TRUE)))

message("\nPotential cell type columns:")
print(celltype_candidates)

# Select the reference cell type column
if (length(celltype_candidates) > 0) {
  ref_celltype_col <- celltype_candidates[1]
  message(paste("\nUsing reference column:", ref_celltype_col))
  
  # Show distribution
  message("\nReference cell type distribution:")
  print(table(GSE115978@meta.data[[ref_celltype_col]]))
} else {
  # Fallback to melanoma_state if no other annotation found
  ref_celltype_col <- "melanoma_state"
  message(paste("\nNo original annotations found. Using:", ref_celltype_col))
  message("\nReference cell type distribution:")
  print(table(GSE115978@meta.data[[ref_celltype_col]]))
}

# ============================================================================
# STEP 3: Prepare Reference Dataset
# ============================================================================

message("\n=== Preparing Reference Dataset ===")

# Ensure reference is processed
if (!"pca" %in% names(GSE115978@reductions)) {
  message("Processing reference dataset...")
  GSE115978 <- NormalizeData(GSE115978)
  GSE115978 <- FindVariableFeatures(GSE115978, nfeatures = 2000)
  GSE115978 <- ScaleData(GSE115978)
  GSE115978 <- RunPCA(GSE115978)
  GSE115978 <- RunUMAP(GSE115978, dims = 1:30)
  message("Reference processing complete")
} else {
  message("✓ Reference already processed")
}

# ============================================================================
# STEP 4: Prepare Query Dataset
# ============================================================================

message("\n=== Preparing Query Dataset ===")

# Ensure query is processed
if (!"pca" %in% names(GSE164897@reductions)) {
  message("Processing query dataset...")
  GSE164897 <- NormalizeData(GSE164897)
  GSE164897 <- FindVariableFeatures(GSE164897, nfeatures = 2000)
  GSE164897 <- ScaleData(GSE164897)
  GSE164897 <- RunPCA(GSE164897)
  message("Query processing complete")
} else {
  message("✓ Query already processed")
}

# ============================================================================
# STEP 5: Find Transfer Anchors
# ============================================================================

message("\n=== Finding Transfer Anchors ===")
message("This may take several minutes...")

# Find anchors between reference and query
transfer_anchors <- FindTransferAnchors(
  reference = GSE115978,
  query = GSE164897,
  dims = 1:30,
  reference.reduction = "pca",
  normalization.method = "LogNormalize"
)

message(paste("Found", nrow(transfer_anchors@anchors), "anchors"))

# ============================================================================
# STEP 6: Transfer Cell Type Labels
# ============================================================================

message("\n=== Transferring Cell Type Labels ===")

# Transfer data
predictions <- TransferData(
  anchorset = transfer_anchors,
  refdata = GSE115978@meta.data[[ref_celltype_col]],
  dims = 1:30
)

# Add predictions to query object
GSE164897$predicted_celltype <- predictions$predicted.id
GSE164897$prediction_score <- predictions$prediction.score.max

# Get prediction scores for each cell type
pred_score_cols <- grep("prediction.score", colnames(predictions), value = TRUE)
for (col in pred_score_cols) {
  GSE164897@meta.data[[col]] <- predictions[[col]]
}

message("Cell type transfer complete!")

# ============================================================================
# STEP 7: Evaluate Transfer Quality
# ============================================================================

message("\n=== Transfer Quality Metrics ===")

# Prediction score summary
message("\nPrediction Score Summary:")
message(paste("Mean prediction score:", round(mean(GSE164897$prediction_score), 3)))
message(paste("Median prediction score:", round(median(GSE164897$prediction_score), 3)))
message(paste("Min prediction score:", round(min(GSE164897$prediction_score), 3)))
message(paste("Max prediction score:", round(max(GSE164897$prediction_score), 3)))

# Cells with low confidence (score < 0.5)
low_conf <- sum(GSE164897$prediction_score < 0.5)
message(paste("\nCells with low confidence (score < 0.5):", low_conf, 
              paste0("(", round(100*low_conf/ncol(GSE164897), 1), "%)")))

# Predicted cell type distribution
message("\n=== Predicted Cell Type Distribution ===")
pred_table <- table(GSE164897$predicted_celltype)
print(pred_table)

# Distribution by treatment
message("\n=== Predicted Cell Types by Treatment ===")
pred_by_treatment <- table(GSE164897$predicted_celltype, GSE164897$treatment)
print(pred_by_treatment)

# ============================================================================
# STEP 8: Visualizations
# ============================================================================

message("\n=== Creating Visualizations ===")

# 1. UMAP with predicted cell types
p1 <- DimPlot(GSE164897, reduction = "umap", group.by = "predicted_celltype", 
              label = TRUE, repel = TRUE) +
  ggtitle("GSE164897 - Predicted Cell Types (from GSE115978)") +
  theme_minimal()
print(p1)

# 2. Split by treatment
p2 <- DimPlot(GSE164897, reduction = "umap", group.by = "predicted_celltype", 
              split.by = "treatment", label = TRUE, repel = TRUE, ncol = 2) +
  ggtitle("Predicted Cell Types by Treatment") +
  theme_minimal()
print(p2)

# 3. Prediction score distribution
p3 <- ggplot(GSE164897@meta.data, aes(x = prediction_score)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
  labs(title = "Prediction Score Distribution",
       x = "Prediction Score", y = "Number of Cells") +
  theme_minimal()
print(p3)

# 4. Prediction score by cell type
p4 <- ggplot(GSE164897@meta.data, aes(x = predicted_celltype, y = prediction_score)) +
  geom_violin(fill = "lightblue") +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "Prediction Score by Cell Type",
       x = "Predicted Cell Type", y = "Prediction Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p4)

# 5. Compare with melanoma states
if ("celltype" %in% colnames(GSE164897@meta.data)) {
  message("\n=== Comparing Predicted Types with Melanoma States ===")
  comparison <- table(GSE164897$predicted_celltype, GSE164897$celltype)
  print(comparison)
  
  # Heatmap of comparison
  library(pheatmap)
  comparison_prop <- prop.table(comparison, margin = 2) * 100
  
  pheatmap(comparison_prop,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           main = "Predicted Cell Type vs Melanoma State (%)",
           display_numbers = TRUE,
           number_format = "%.1f",
           fontsize = 10,
           angle_col = "45")
}

# 6. Feature plot of prediction scores
p5 <- FeaturePlot(GSE164897, features = "prediction_score", reduction = "umap") +
  scale_color_viridis_c() +
  ggtitle("Prediction Score on UMAP") +
  theme_minimal()
print(p5)

# ============================================================================
# STEP 9: Summary and Recommendations
# ============================================================================

message("\n=== Transfer Summary ===")
message(paste("Reference dataset:", "GSE115978"))
message(paste("Query dataset:", "GSE164897"))
message(paste("Reference annotation:", ref_celltype_col))
message(paste("Cells annotated:", ncol(GSE164897)))
message(paste("Cell types identified:", length(unique(GSE164897$predicted_celltype))))
message(paste("Mean prediction confidence:", round(mean(GSE164897$prediction_score), 3)))

message("\n=== Annotations Stored ===")
message("GSE164897$predicted_celltype - Transferred cell type labels")
message("GSE164897$prediction_score - Prediction confidence (0-1)")
message("GSE164897$prediction.score.* - Per-cell-type scores")

message("\n=== Recommendations ===")
if (mean(GSE164897$prediction_score) > 0.7) {
  message("✓ High confidence predictions (mean > 0.7)")
  message("  The transferred labels are reliable")
} else if (mean(GSE164897$prediction_score) > 0.5) {
  message("⚠ Moderate confidence predictions (mean 0.5-0.7)")
  message("  Consider filtering low-confidence cells")
} else {
  message("⚠ Low confidence predictions (mean < 0.5)")
  message("  Datasets may be too different for reliable transfer")
  message("  Consider using melanoma_state annotations instead")
}

message("\n=== Reference-Based Annotation Complete! ===")