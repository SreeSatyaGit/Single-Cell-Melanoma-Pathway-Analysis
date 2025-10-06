library(vegan)
library(dplyr)
library(devtools)
library(oposSOM)
library(Seurat)
library(SeuratData)
library(patchwork)
library(biomaRt)
library(googledrive)


# Combined Seurat object setup
combined_seurat <- merge(
  x = s,
  y = list(v),
  add.cell.ids = c("S", "RV"),
  project = "CombinedDataset"
)


combined_seurat <- NormalizeData(combined_seurat)
combined_seurat <- FindVariableFeatures(combined_seurat)
combined_seurat <- ScaleData(combined_seurat)
combined_seurat <- RunPCA(combined_seurat)

combined_seurat <- IntegrateLayers(
  object = combined_seurat,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)

combined_seurat <- FindNeighbors(combined_seurat, reduction = "integrated.cca", dims = 1:30)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.1, cluster.name = "cca_clusters")
combined_seurat <- RunUMAP(combined_seurat, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

p1 <- DimPlot(combined_seurat, reduction = "umap.cca", group.by = c("orig.ident", "cca_clusters"))
p2 <- DimPlot(combined_seurat, reduction = "umap.cca", split.by = "cca_clusters", group.by = "orig.ident") 
p3 <- DimPlot(combined_seurat, reduction = "umap.cca", split.by = "orig.ident", group.by = "cca_clusters") 

p1
# Identify most common and highly expressed genes in each cluster
# Ensure clustering has been performed
combined_seurat <- FindClusters(combined_seurat, resolution = 0.1)

# Calculate average expression per cluster

average_expression <- AverageExpression(combined_seurat, group.by = "orig.ident")

# Extract average RNA expression data
average_rna <- average_expression$RNA

# Convert to a data frame
average_rna_df <- as.data.frame(as.matrix(average_rna))

# Find top 10 highly expressed genes per cluster
top_genes_per_cluster <- lapply(
  colnames(average_rna_df),
  function(cluster) {
    cluster_data <- average_rna_df[, cluster, drop = FALSE]  # Select data for the cluster
    cluster_data <- data.frame(gene = rownames(cluster_data), expression = cluster_data[, 1])  # Add gene names
    sorted_genes <- cluster_data[order(-cluster_data$expression), ]  # Sort by expression
    head(sorted_genes, 10)  # Return top 10 genes
  }
)
names(top_genes_per_cluster) <- colnames(average_rna_df)

# View results
for (cluster_id in names(top_genes_per_cluster)) {
  cat("Cluster", cluster_id, "Top Genes:\n")
  print(top_genes_per_cluster[[cluster_id]])
  cat("\n\n")
}

#oposSOM Analysis

# Create a new opossom environment
env <- opossom.new(list(
  dataset.name = "SandRVResutls",
  dim.1stLvlSom = 40
))

# Example: Assign normalized data to the environment
env$indata <- as.matrix(log(average_rna + 1))  # Add 1 to avoid log(0)




group.info <- data.frame( 
  group.labels = c("Senstive",
                   "RV"),
  
  group.colors = c(
    "red2",
    "brown"),
  
  row.names=colnames(average_rna_df))
# Run the oposSOM pipeline
opossom.run(env)




# Alpha diversity measures
metadata <- combined_seurat@meta.data
cluster_counts <- table(metadata$seurat_clusters)
species_richness <- length(cluster_counts)
proportions <- cluster_counts / sum(cluster_counts)
shannon_index <- -sum(proportions * log(proportions))
simpson_index <- 1 - sum(proportions^2)
abundance_matrix <- as.matrix(cluster_counts)
shannon_index_vegan <- diversity(abundance_matrix, index = "shannon")
simpson_index_vegan <- diversity(abundance_matrix, index = "simpson")
