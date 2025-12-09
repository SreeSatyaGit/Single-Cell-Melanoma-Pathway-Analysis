# scType Cell Type Annotation
# Automated marker-based annotation for melanoma dataset

library(Seurat)
library(dplyr)
library(ggplot2)

# Install required packages
if (!requireNamespace("HGNChelper", quietly = TRUE))
  install.packages("HGNChelper")

if (!requireNamespace("openxlsx", quietly = TRUE))
  install.packages("openxlsx")

library(HGNChelper)
library(openxlsx)

# Check if GSE164897 exists
if (!exists("GSE164897")) {
  stop("ERROR: GSE164897 Seurat object not found. Please run GSE164897.R first.")
}

message("\n=== scType Cell Type Annotation ===")

# Load scType functions
message("Loading scType functions...")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# Load gene set database
message("Loading scType database...")
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"

# Download and check available tissues
temp_db <- tempfile(fileext = ".xlsx")
download.file(db_, temp_db, mode = "wb", quiet = TRUE)

# Read all sheets to see available tissues
available_sheets <- getSheetNames(temp_db)
message("\nAvailable tissue types in scType database:")
print(available_sheets)

# Try different tissue options for melanoma
tissue_options <- c("Skin", "Immune system", "Melanoma")
tissue <- NULL

for (tissue_opt in tissue_options) {
  if (tissue_opt %in% available_sheets) {
    tissue <- tissue_opt
    message(paste("\nUsing tissue type:", tissue))
    break
  }
}

if (is.null(tissue)) {
  # Default to first available tissue
  tissue <- available_sheets[1]
  message(paste("\nMelanoma-specific tissue not found. Using:", tissue))
  message("Note: Results may not be optimal for melanoma cells")
}

# Prepare gene sets
message("\nPreparing gene sets...")
gs_list <- tryCatch({
  gene_sets_prepare(db_, tissue)
}, error = function(e) {
  message("Error loading gene sets for tissue:", tissue)
  message("Trying alternative approach...")
  
  # Manual gene set preparation as fallback
  db_data <- read.xlsx(temp_db, sheet = tissue)
  
  # Create gene sets manually
  gs_positive <- list()
  gs_negative <- list()
  
  if (nrow(db_data) > 0) {
    for (i in 1:nrow(db_data)) {
      cell_type <- db_data[i, "cellName"]
      genes_pos <- unlist(strsplit(as.character(db_data[i, "geneSymbolmore1"]), ","))
      genes_neg <- unlist(strsplit(as.character(db_data[i, "geneSymbolmore2"]), ","))
      
      genes_pos <- genes_pos[genes_pos != "" & !is.na(genes_pos)]
      genes_neg <- genes_neg[genes_neg != "" & !is.na(genes_neg)]
      
      if (length(genes_pos) > 0) {
        gs_positive[[cell_type]] <- genes_pos
      }
      if (length(genes_neg) > 0) {
        gs_negative[[cell_type]] <- genes_neg
      }
    }
  }
  
  list(gs_positive = gs_positive, gs_negative = gs_negative)
})

message(paste("Loaded", length(gs_list$gs_positive), "cell type gene sets"))

# Check if we have gene sets
if (length(gs_list$gs_positive) == 0) {
  stop("No gene sets loaded. Please check the tissue type or database.")
}

# Ensure data is scaled
message("\nPreparing expression data...")
if (!"scale.data" %in% names(GSE164897@assays$RNA@layers)) {
  message("Scaling data...")
  GSE164897 <- ScaleData(GSE164897, features = rownames(GSE164897))
}

# Calculate scType scores
message("\nCalculating scType scores...")
es.max <- tryCatch({
  sctype_score(
    scRNAseqData = as.matrix(GetAssayData(GSE164897, layer = "scale.data")),
    scaled = TRUE,
    gs = gs_list$gs_positive,
    gs2 = gs_list$gs_negative
  )
}, error = function(e) {
  message("Error with scale.data, trying data layer...")
  sctype_score(
    scRNAseqData = as.matrix(GetAssayData(GSE164897, layer = "data")),
    scaled = FALSE,
    gs = gs_list$gs_positive,
    gs2 = gs_list$gs_negative
  )
})

message("scType scoring complete!")

# Assign cell types per cluster
message("\nAssigning cell types to clusters...")
cL_results <- do.call("rbind", lapply(unique(GSE164897$seurat_clusters), function(cl){
  es.max.cl <- sort(rowSums(es.max[, GSE164897$seurat_clusters == cl, drop = FALSE]), decreasing = TRUE)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, stringsAsFactors = FALSE), 10)
}))

# Get top cell type per cluster
sctype_scores <- cL_results %>% 
  group_by(cluster) %>% 
  top_n(n = 1, wt = scores)

# Add to Seurat object
GSE164897$celltype_sctype <- sctype_scores$type[match(GSE164897$seurat_clusters, sctype_scores$cluster)]

# Summary
message("\n=== scType Results ===")
message("Cell type assignments by cluster:")
print(sctype_scores)

message("\nCell type distribution:")
print(table(GSE164897$celltype_sctype))

message("\nCell types by treatment:")
print(table(GSE164897$celltype_sctype, GSE164897$treatment))

# Visualizations
message("\n=== Creating Visualizations ===")

# UMAP with scType annotations
p1 <- DimPlot(GSE164897, group.by = "celltype_sctype", label = TRUE, repel = TRUE) +
  ggtitle(paste("scType Cell Type Annotation (", tissue, ")", sep = "")) +
  theme_minimal()
print(p1)

# Split by treatment
p2 <- DimPlot(GSE164897, group.by = "celltype_sctype", split.by = "treatment",
              label = TRUE, repel = TRUE, ncol = 2) +
  ggtitle("scType Cell Types by Treatment") +
  theme_minimal()
print(p2)

# Score heatmap
message("\nCreating score heatmap...")
top_types <- head(unique(sctype_scores$type), 20)
if (length(top_types) > 0) {
  score_matrix <- es.max[top_types, , drop = FALSE]
  
  # Average scores by cluster
  avg_scores <- sapply(unique(GSE164897$seurat_clusters), function(cl) {
    rowMeans(score_matrix[, GSE164897$seurat_clusters == cl, drop = FALSE])
  })
  colnames(avg_scores) <- unique(GSE164897$seurat_clusters)
  
  library(pheatmap)
  pheatmap(avg_scores,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           main = "scType Scores by Cluster",
           fontsize = 10)
}

message("\n=== scType Annotation Complete ===")
message("Cell type labels stored in: GSE164897$celltype_sctype")
message(paste("Tissue database used:", tissue))

# Clean up
unlink(temp_db)