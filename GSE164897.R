library(data.table)
library(Matrix)
library(Seurat)

# Point to your directory
base_dir <- "/projects/vanaja_lab/satya/R/SCTMR"  # <-- change this


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


# run standard anlaysis workflow
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

# re-join layers after integration
GSE164897[["RNA"]] <- JoinLayers(GSE164897[["RNA"]])

GSE164897 <- FindNeighbors(GSE164897, reduction = "integrated.cca", dims = 1:30)
GSE164897 <- FindClusters(GSE164897, resolution = 1)

GSE164897 <- RunUMAP(GSE164897, dims = 1:30, reduction = "integrated.cca")
# Visualization
map <- c(
  "GSM5022595" = "Vemurafenib",
  "GSM5022596" = "untreated",
  "GSM5022597" = "vem_cob",
  "GSM5022598" = "vem_tram"
)

# Add new meta column (e.g., 'treatment')
GSE164897$treatment <- unname(map[as.character(GSE164897$orig.ident)])

# Make it an ordered factor (optional but recommended)
GSE164897$treatment <- factor(
  GSE164897$treatment,
  levels = c("untreated", "Vemurafenib", "vem_cob", "vem_tram")
)

DimPlot(GSE164897, reduction = "umap", group.by = c("seurat_clusters","treatment"))
