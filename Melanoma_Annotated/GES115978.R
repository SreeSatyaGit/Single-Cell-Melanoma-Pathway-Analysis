# GSE115978 - Melanoma scRNA-seq Data Processing
# Pre-annotated melanoma dataset with cell type labels
# This script processes the GSE115978 dataset with existing cell annotations

library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(data.table)


message("=== GSE115978 Data Loading ===")
message(paste("Working directory:", getwd()))

# ============================================================================
# STEP 1: Load Data
# ============================================================================

message("\n=== Loading Count Matrix ===")
# Read counts (this may take a few minutes for large files)
counts <- fread("/projects/vanaja_lab/satya/SCPA/Melanoma_Annotated/GSE115978_counts.csv.gz", data.table = FALSE)
message(paste("Loaded counts matrix:", nrow(counts), "genes x", ncol(counts)-1, "cells"))

# First column is gene names
gene_names <- counts[, 1]
counts <- counts[, -1]
rownames(counts) <- gene_names

# Convert to sparse matrix for memory efficiency
library(Matrix)
counts <- Matrix(as.matrix(counts), sparse = TRUE)
message("Converted to sparse matrix")

message("\n=== Loading Cell Annotations ===")
# Read cell annotations
annotations <- fread("/projects/vanaja_lab/satya/SCPA/Melanoma_Annotated/GSE115978_cell.annotations.csv.gz", data.table = FALSE)
message(paste("Loaded annotations for", nrow(annotations), "cells"))
message("Available annotation columns:")
print(colnames(annotations))

# Display first few rows
message("\nFirst few annotation rows:")
print(head(annotations))

# ============================================================================
# STEP 2: Create Seurat Object
# ============================================================================

message("\n=== Creating Seurat Object ===")

# Ensure cell names match between counts and annotations
if ("cell" %in% colnames(annotations)) {
  rownames(annotations) <- annotations$cell
} else if ("Cell" %in% colnames(annotations)) {
  rownames(annotations) <- annotations$Cell
} else {
  # Use first column as cell names
  rownames(annotations) <- annotations[, 1]
}

# Match cells between counts and annotations
common_cells <- intersect(colnames(counts), rownames(annotations))
message(paste("Common cells between counts and annotations:", length(common_cells)))

# Subset to common cells
counts <- counts[, common_cells]
annotations <- annotations[common_cells, ]

# Create Seurat object
GSE115978 <- CreateSeuratObject(
  counts = counts,
  meta.data = annotations,
  project = "GSE115978",
  min.cells = 3,
  min.features = 200
)

message(paste("Created Seurat object with", ncol(GSE115978), "cells and", nrow(GSE115978), "genes"))

# ============================================================================
# STEP 3: Quality Control
# ============================================================================

message("\n=== Quality Control Metrics ===")

# Ensure nFeature_RNA and nCount_RNA exist
# These should be created by CreateSeuratObject, but if metadata was passed, they might be missing
if (!"nFeature_RNA" %in% colnames(GSE115978@meta.data)) {
  message("Adding nFeature_RNA and nCount_RNA columns...")
  GSE115978$nFeature_RNA <- colSums(GetAssayData(GSE115978, layer = "counts") > 0)
  GSE115978$nCount_RNA <- colSums(GetAssayData(GSE115978, layer = "counts"))
}

# Calculate mitochondrial percentage
GSE115978[["percent.mt"]] <- PercentageFeatureSet(GSE115978, pattern = "^MT-")

# QC metrics summary
message("\nQC Metrics Summary:")
message(paste("Median genes per cell:", median(GSE115978$nFeature_RNA)))
message(paste("Median UMIs per cell:", median(GSE115978$nCount_RNA)))
message(paste("Median % mitochondrial:", round(median(GSE115978$percent.mt, na.rm = TRUE), 2)))

# Visualize QC metrics
p_qc <- VlnPlot(GSE115978, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(p_qc)

# Feature scatter
p_scatter <- FeatureScatter(GSE115978, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm")
print(p_scatter)

# Optional: Filter cells based on QC (adjust thresholds as needed)
message("\n=== Filtering Cells (Optional) ===")
message("Current cell count:", ncol(GSE115978))

# Uncomment to apply filters:
# GSE115978 <- subset(GSE115978, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
# message("After filtering:", ncol(GSE115978))

# ============================================================================
# STEP 4: Normalization and Scaling
# ============================================================================

message("\n=== Normalization and Scaling ===")

GSE115978 <- NormalizeData(GSE115978, normalization.method = "LogNormalize", scale.factor = 10000)
message("Data normalized")

GSE115978 <- FindVariableFeatures(GSE115978, selection.method = "vst", nfeatures = 2000)
message("Variable features identified")

# Top 10 variable genes
top10 <- head(VariableFeatures(GSE115978), 10)
message("Top 10 variable genes:")
print(top10)

GSE115978 <- ScaleData(GSE115978)
message("Data scaled")

# ============================================================================
# STEP 5: Dimensionality Reduction
# ============================================================================

message("\n=== Dimensionality Reduction ===")

GSE115978 <- RunPCA(GSE115978, features = VariableFeatures(object = GSE115978))
message("PCA complete")

# Visualize PCA
print(DimPlot(GSE115978, reduction = "pca"))

# Elbow plot
print(ElbowPlot(GSE115978, ndims = 50))

# ============================================================================
# STEP 6: Clustering
# ============================================================================

message("\n=== Clustering ===")

GSE115978 <- FindNeighbors(GSE115978, dims = 1:30)
GSE115978 <- FindClusters(GSE115978, resolution = 0.8)
message("Clustering complete")

# ============================================================================
# STEP 7: UMAP
# ============================================================================

message("\n=== Running UMAP ===")

GSE115978 <- RunUMAP(GSE115978, dims = 1:30)
message("UMAP complete")

# Visualize by clusters
p_clusters <- DimPlot(GSE115978, reduction = "umap", label = TRUE) +
  ggtitle("GSE115978 - Seurat Clusters")
print(p_clusters)

# ============================================================================
# STEP 8: Visualize Existing Annotations
# ============================================================================

message("\n=== Visualizing Existing Cell Type Annotations ===")

# Find cell type column (common names)
celltype_cols <- grep("celltype|cell_type|CellType|annotation|cluster", 
                      colnames(GSE115978@meta.data), 
                      value = TRUE, ignore.case = TRUE)

message("Found potential cell type columns:")
print(celltype_cols)

# Visualize each annotation column
for (col in celltype_cols) {
  if (col %in% colnames(GSE115978@meta.data)) {
    message(paste("\nVisualizing:", col))
    print(table(GSE115978@meta.data[[col]]))
    
    p <- DimPlot(GSE115978, reduction = "umap", group.by = col, label = TRUE, repel = TRUE) +
      ggtitle(paste("GSE115978 -", col))
    print(p)
  }
}

# ============================================================================
# STEP 9: Add Melanoma Cell State Annotations
# ============================================================================

message("\n=== Adding Melanoma Cell State Annotations ===")

# Define melanoma-specific marker gene sets
melanoma_states <- list(
  Differentiated = c("MLANA", "TYR", "DCT", "PMEL", "TYRP1", "GPR143", "OCA2"),
  Undifferentiated = c("SOX10", "MITF", "PAX3", "SOX9", "FOXD3"),
  Invasive = c("AXL", "NGFR", "VIM", "TWIST1", "ZEB1", "SNAI2", "FN1", "CDH2"),
  Proliferative = c("MKI67", "TOP2A", "PCNA", "CDK1", "CCNB1", "CCNA2"),
  Resistant = c("EGFR", "PDGFRB", "AXL", "NGFR", "WNT5A", "JUN", "FOS"),
  Hypoxic = c("HIF1A", "VEGFA", "SLC2A1", "LDHA", "PDK1", "BNIP3"),
  ImmuneEvasion = c("CD274", "PDCD1LG2", "IDO1", "HAVCR2", "LAG3")
)

# Filter to available genes
melanoma_states_filtered <- lapply(melanoma_states, function(markers) {
  intersect(markers, rownames(GSE115978))
})
melanoma_states_filtered <- melanoma_states_filtered[sapply(melanoma_states_filtered, length) > 0]

message(paste("Calculating module scores for", length(melanoma_states_filtered), "melanoma states..."))

# Remove existing score columns
existing_scores <- grep("_score", colnames(GSE115978@meta.data), value = TRUE)
if (length(existing_scores) > 0) {
  GSE115978@meta.data <- GSE115978@meta.data[, !colnames(GSE115978@meta.data) %in% existing_scores]
}

# Calculate module scores
for (state in names(melanoma_states_filtered)) {
  if (length(melanoma_states_filtered[[state]]) > 0) {
    GSE115978 <- AddModuleScore(
      GSE115978,
      features = list(melanoma_states_filtered[[state]]),
      name = paste0(state, "_score"),
      seed = 42
    )
  }
}

# Rename score columns
score_cols <- grep("_score1$", colnames(GSE115978@meta.data), value = TRUE)
for (i in seq_along(names(melanoma_states_filtered))) {
  old_name <- paste0(names(melanoma_states_filtered)[i], "_score1")
  new_name <- paste0(names(melanoma_states_filtered)[i], "_score")
  if (old_name %in% colnames(GSE115978@meta.data)) {
    colnames(GSE115978@meta.data)[colnames(GSE115978@meta.data) == old_name] <- new_name
  }
}

score_cols <- paste0(names(melanoma_states_filtered), "_score")

# Assign dominant state
state_scores <- GSE115978@meta.data[, score_cols, drop = FALSE]
GSE115978$melanoma_state <- colnames(state_scores)[apply(state_scores, 1, which.max)]
GSE115978$melanoma_state <- gsub("_score", "", GSE115978$melanoma_state)
GSE115978$melanoma_state_score <- apply(state_scores, 1, max)

message("\n=== Melanoma State Distribution ===")
print(table(GSE115978$melanoma_state))

# Visualize melanoma states
p_melanoma <- DimPlot(GSE115978, reduction = "umap", group.by = "melanoma_state", 
                      label = TRUE, repel = TRUE) +
  ggtitle("GSE115978 - Melanoma Cell States")
print(p_melanoma)

# ============================================================================
# STEP 10: Summary and Save
# ============================================================================

message("\n=== Processing Complete ===")
message(paste("Final cell count:", ncol(GSE115978)))
message(paste("Final gene count:", nrow(GSE115978)))
message(paste("Clusters identified:", length(unique(GSE115978$seurat_clusters))))
message(paste("Melanoma states identified:", length(unique(GSE115978$melanoma_state))))

# Display available metadata
message("\n=== Available Metadata Columns ===")
print(colnames(GSE115978@meta.data))

# Save the processed object (optional)
# saveRDS(GSE115978, "GSE115978_processed.rds")
# message("Saved processed object to: GSE115978_processed.rds")

message("\n=== GSE115978 Processing Complete! ===")
message("Seurat object: GSE115978")
message("Use DimPlot(GSE115978, group.by = 'melanoma_state') to visualize cell states")
