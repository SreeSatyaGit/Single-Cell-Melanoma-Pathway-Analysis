library(data.table)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(org.Hs.eg.db)
# Point to your directory
base_dir <- "/projects/vanaja_lab/satya/SCTMR"  


files <- list.files(
  base_dir,
  pattern = "^GSM\\d+_raw_counts_probe\\d+\\.txt\\.gz$",
  full.names = TRUE
)

# Function: read gzipped counts and make Seurat object
make_seurat_from_file <- function(path) {
  # Read file
  dt <- fread(path)
  
  # First column = gene IDs
  genes <- dt[[1]]
  dt[, 1 := NULL]   # drop gene column
  
  # Convert to matrix
  mat <- as.matrix(dt)
  rownames(mat) <- genes
  
  # Ensure unique colnames (cells)
  sample_id <- sub("\\.txt\\.gz$", "", basename(path))
  colnames(mat) <- paste0(sample_id, "_", make.unique(colnames(mat)))
  
  # Make sparse matrix
  mat <- Matrix(mat, sparse = TRUE)
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = mat,
    project = sample_id,
    min.cells = 3,
    min.features = 200
  )
  
  return(seurat_obj)
}


seurat_list <- lapply(files, make_seurat_from_file)
names(seurat_list) <- basename(files)

GSE164897 <- Reduce(function(x, y) merge(x = x, y = y), seurat_list)


# Treatment mapping
map <- c(
  "GSM5022595" = "Vemurafenib",
  "GSM5022596" = "untreated",
  "GSM5022597" = "vem_cob",
  "GSM5022598" = "vem_tram"
)

# Extract GSM ID from orig.ident (e.g., "GSM5022595_raw_counts_probe3" -> "GSM5022595")
gsm_ids <- sub("^(GSM\\d+).*$", "\\1", as.character(GSE164897$orig.ident))

# Add treatment column based on extracted GSM ID
GSE164897$treatment <- unname(map[gsm_ids])

# Make it an ordered factor
GSE164897$treatment <- factor(
  GSE164897$treatment,
  levels = c("untreated", "Vemurafenib", "vem_cob", "vem_tram")
)

# Check sample sizes per treatment group
cell_counts <- table(GSE164897$treatment, useNA = "ifany")
print("=== Cell counts per treatment group (before downsampling) ===")
print(cell_counts)

# Find minimum sample size across all groups - this ensures all groups can be equal
non_zero_counts <- cell_counts[cell_counts > 0 & !is.na(names(cell_counts))]
if (length(non_zero_counts) == 0) {
  stop("ERROR: All treatment groups have 0 cells after mapping.")
}
target_size <- min(non_zero_counts)
print(paste("\nTarget sample size (minimum across all groups):", target_size))
print(paste("This will downsample all groups to", target_size, "cells"))

# Downsample ALL treatment groups to match the minimum size
set.seed(42)  # For reproducibility
cells_to_keep <- c()

for (trt in levels(GSE164897$treatment)) {
  cells_in_group <- colnames(GSE164897)[GSE164897$treatment == trt & !is.na(GSE164897$treatment)]
  
  if (length(cells_in_group) > target_size) {
    # Downsample if group is larger than target
    sampled_cells <- sample(cells_in_group, size = target_size, replace = FALSE)
    cells_to_keep <- c(cells_to_keep, sampled_cells)
    print(paste("Downsampled", trt, "from", length(cells_in_group), "to", target_size, "cells"))
  } else if (length(cells_in_group) == target_size) {
    # Keep all if exactly at target size
    cells_to_keep <- c(cells_to_keep, cells_in_group)
    print(paste("Kept", trt, "at", length(cells_in_group), "cells (already at target size)"))
  } else if (length(cells_in_group) > 0) {
    # For groups smaller than target, keep all (but this shouldn't happen with min)
    cells_to_keep <- c(cells_to_keep, cells_in_group)
    print(paste("WARNING: Kept", trt, "at", length(cells_in_group), "cells (smaller than target - this will cause unequal sizes)"))
  } else {
    print(paste("Warning:", trt, "has 0 cells - skipping"))
  }
}

print(paste("Total cells to keep:", length(cells_to_keep)))

# Ensure cells_to_keep matches actual colnames
all_cells <- colnames(GSE164897)
cells_to_keep <- intersect(cells_to_keep, all_cells)
print(paste("After matching with colnames, cells to keep:", length(cells_to_keep)))

if (length(cells_to_keep) == 0) {
  stop("ERROR: No cells match after intersection. Check treatment mapping.")
}

# Subset Seurat object - use Seurat v5 compatible method
# Extract counts and metadata, then recreate object
print("Subsetting Seurat object...")

# Get layer names
layer_names <- Layers(GSE164897)
print(paste("Available layers:", paste(layer_names, collapse = ", ")))

# Extract counts - handle multiple layers in Seurat v5
# First, join layers to get a single counts matrix
tryCatch({
  # Join all layers into a single counts matrix
  GSE164897[["RNA"]] <- JoinLayers(GSE164897[["RNA"]])
  print("Joined layers successfully")
}, error = function(e) {
  print(paste("JoinLayers warning (may already be joined):", e$message))
})

# Now extract counts - try different methods
counts_subset <- NULL

# Method 1: Try GetAssayData after joining layers
tryCatch({
  counts_subset <- GetAssayData(GSE164897, slot = "counts", assay = "RNA")[, cells_to_keep, drop = FALSE]
  print("Extracted counts using GetAssayData")
}, error = function(e) {
  print(paste("GetAssayData failed:", e$message))
})

# Method 2: Try LayerData if multiple layers still exist
if (is.null(counts_subset)) {
  tryCatch({
    if (length(layer_names) > 0) {
      # Get all counts layers
      counts_layers <- grep("counts", layer_names, value = TRUE)
      if (length(counts_layers) > 0) {
        # Extract from first layer and subset
        counts_subset <- LayerData(GSE164897, layer = counts_layers[1], assay = "RNA")[, cells_to_keep, drop = FALSE]
        print(paste("Extracted counts from layer:", counts_layers[1]))
      }
    }
  }, error = function(e) {
    print(paste("LayerData failed:", e$message))
  })
}

# Method 3: Direct access as last resort
if (is.null(counts_subset)) {
  tryCatch({
    counts_subset <- GSE164897[["RNA"]]$counts[, cells_to_keep, drop = FALSE]
    print("Extracted counts using direct access")
  }, error = function(e) {
    stop(paste("All methods to extract counts failed:", e$message, 
               "\nPlease check that cells_to_keep matches colnames(GSE164897)"))
  })
}

# Verify counts columns match cells_to_keep
if (!identical(colnames(counts_subset), cells_to_keep)) {
  print("WARNING: Counts columns don't match cells_to_keep order - reordering...")
  counts_subset <- counts_subset[, cells_to_keep, drop = FALSE]
}

# Get metadata - ensure rows match cells_to_keep exactly
meta_subset <- GSE164897@meta.data[cells_to_keep, , drop = FALSE]

# Verify metadata rows match counts columns
if (nrow(meta_subset) != ncol(counts_subset)) {
  stop(paste("ERROR: Metadata rows (", nrow(meta_subset), ") don't match counts columns (", ncol(counts_subset), ")"))
}

# Ensure treatment column exists and is preserved
if (!"treatment" %in% colnames(meta_subset)) {
  stop("ERROR: treatment column missing from metadata!")
}

# Recreate Seurat object with subsetted data
GSE164897 <- CreateSeuratObject(counts = counts_subset, meta.data = meta_subset)
print("Seurat object recreated with downsampled cells")

# Ensure treatment factor is preserved
GSE164897$treatment <- factor(
  GSE164897$treatment,
  levels = c("untreated", "Vemurafenib", "vem_cob", "vem_tram")
)


print("Splitting object into layers for integration...")
GSE164897_split <- SplitObject(GSE164897, split.by = "orig.ident")
print(paste("Split into", length(GSE164897_split), "objects"))

# Merge back - this will create layers in Seurat v5
GSE164897 <- merge(GSE164897_split[[1]], y = GSE164897_split[-1])
print("Merged back into single object with layers")

# Verify equal sample sizes
cell_counts_after <- table(GSE164897$treatment)
print("\n=== Cell counts per treatment group (after downsampling) ===")
print(cell_counts_after)

# Check if all groups have equal sizes
unique_counts <- unique(as.numeric(cell_counts_after))
if (length(unique_counts) == 1) {
  print(paste("✓ SUCCESS: All treatment groups have equal size:", unique_counts[1], "cells"))
} else {
  print("WARNING: Treatment groups still have unequal sizes!")
  print("Expected size:", target_size)
  print("Actual sizes:", paste(paste(names(cell_counts_after), cell_counts_after, sep = "="), collapse = ", "))
}


GSE164897 <- NormalizeData(GSE164897)
GSE164897 <- FindVariableFeatures(GSE164897)
GSE164897 <- ScaleData(GSE164897)
GSE164897 <- RunPCA(GSE164897)



GSE164897 <- FindNeighbors(GSE164897, dims = 1:30, reduction = "pca")
GSE164897 <- FindClusters(GSE164897, resolution = 2, cluster.name = "unintegrated_clusters")


GSE164897 <- RunUMAP(GSE164897, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(GSE164897, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters"))


GSE164897 <- IntegrateLayers(object = GSE164897, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                               verbose = FALSE)
GSE164897 <- IntegrateLayers(object = GSE164897, method = CCAIntegration, 
                                    orig.reduction = "pca", new.reduction = "integrated.cca",
                                    assay = "RNA",
                                    verbose = FALSE)
  

# re-join layers after integration
GSE164897[["RNA"]] <- JoinLayers(GSE164897[["RNA"]])

GSE164897 <- FindNeighbors(GSE164897, reduction = "integrated.cca", dims = 1:30)
GSE164897 <- FindClusters(GSE164897, resolution = 1)

GSE164897 <- RunUMAP(GSE164897, dims = 1:30, reduction = "integrated.cca")

# Visualization (treatment column already added during downsampling)
DimPlot(GSE164897, reduction = "umap", group.by = c("seurat_clusters","treatment"))





