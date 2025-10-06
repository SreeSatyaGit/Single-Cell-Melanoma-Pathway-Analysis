library(future)
library(Biobase)

combined_BlueGreen <- merge(
  x = clusterBlueSeurat,
  y = list(clusterGreenSeurat),
  add.cell.ids = c("Blue", "Grenn"),
  project = "CombinedBlueandGreen"
)


combined_BlueGreen <- NormalizeData(combined_BlueGreen)
combined_BlueGreen <- FindVariableFeatures(combined_BlueGreen)
combined_BlueGreen <- ScaleData(combined_BlueGreen)
combined_BlueGreen <- RunPCA(combined_BlueGreen)

combined_BlueGreen <- IntegrateLayers(
  object = combined_BlueGreen,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)

combined_BlueGreen <- FindNeighbors(combined_BlueGreen, reduction = "integrated.cca", dims = 1:30)
combined_BlueGreen <- FindClusters(combined_BlueGreen, resolution = 0.1, cluster.name = "cca_clusters")
combined_BlueGreen <- RunUMAP(combined_BlueGreen, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

p1 <- DimPlot(combined_BlueGreen, reduction = "umap.cca", group.by = c("orig.ident", "cca_clusters"))
p2 <- DimPlot(combined_BlueGreen, reduction = "umap.cca", split.by = "cca_clusters", group.by = "orig.ident") 
p3 <- DimPlot(combined_BlueGreen, reduction = "umap.cca", split.by = "orig.ident", group.by = "cca_clusters") 

p1
# Identify most common and highly expressed genes in each cluster
# Ensure clustering has been performed
combined_BlueGreen <- FindClusters(combined_BlueGreen, resolution = 0.1)

# Calculate average expression per cluster

average_expression_Blue <- AverageExpression(combined_BlueGreen, group.by = "orig.ident")

# Extract average RNA expression data
average_rna_Blue <- average_expression_Blue$RNA

# Convert to a data frame
average_rna_Blue_df <- as.data.frame(as.matrix(average_rna_Blue))

# Find top 10 highly expressed genes per cluster
top_genes_per_cluster <- lapply(
  colnames(average_rna_Blue_df),
  function(cluster) {
    cluster_data <- average_rna_Blue_df[, cluster, drop = FALSE]  # Select data for the cluster
    cluster_data <- data.frame(gene = rownames(cluster_data), expression = cluster_data[, 1])  # Add gene names
    sorted_genes <- cluster_data[order(-cluster_data$expression), ]  # Sort by expression
    head(sorted_genes, 10)  # Return top 10 genes
  }
)
names(average_rna_Blue_df) <- colnames(average_rna_Blue_df)

# View results
for (cluster_id in names(top_genes_per_cluster)) {
  cat("Cluster", cluster_id, "Top Genes:\n")
  print(top_genes_per_cluster[[cluster_id]])
  cat("\n\n")
}

#oposSOM Analysis
# Create a new opossom env_blueironment

env <- opossom.new(list(
  dataset.name = "BlueGrennSOM",
  dim.1stLvlSom = 40
))

# Example: Assign normalized data to the env_blueironment
env$indata <- as.matrix(log(average_rna_Blue_df + 1))

# Run the oposSOM pipeline
opossom.run(env)
