library(SCPA)
library(Seurat)
library(tidyverse)
library(magrittr)
library(dyno)
library(ComplexHeatmap)
library(circlize)
library(msigdbr)

source("/projects/vanaja_lab/satya/SCPA/Reproduce.R")
var_genes <- VariableFeatures(GSE164897)

counts <- GetAssayData(GSE164897, layer = "counts", assay = "RNA")
expression <- GetAssayData(GSE164897, layer = "data", assay = "RNA")

counts_var <- t(as.matrix(counts[var_genes, ]))
expression_var <- t(as.matrix(expression[var_genes, ]))

print("Wrapping data for dyno...")
dataset <- wrap_expression(
  counts = counts_var,
  expression = expression_var
)


print("Inferring trajectory using Slingshot...")

model <- infer_trajectory(dataset, method = "slingshot", parameters = list(cluster_method = "clara"))

# Calculate Pseudotime
print("Calculating pseudotime...")
pseudotime <- calculate_pseudotime(model)


GSE164897$pseudotime <- pseudotime[colnames(GSE164897)]

print("Plotting pseudotime and treatment...")
p1 <- plot_dimred(
  model, 
  "pseudotime", 
  pseudotime = pseudotime, 
  hex_cells = FALSE,
  plot_trajectory = TRUE, 
  size_cells = 1, 
  alpha_cells = 0.8
) + 
  ggtitle("Pseudotime") +
  theme(aspect.ratio = 1)

p2 <- plot_dimred(
  model, 
  grouping = GSE164897$treatment,
  plot_trajectory = TRUE, 
  size_cells = 1, 
  alpha_cells = 0.8
) + 
  ggtitle("Treatment") +
  theme(aspect.ratio = 1)

# Combine plots (requires patchwork, usually loaded with Seurat)
p1 + p2


plot_dimred(model, 
            grouping = group_onto_nearest_milestones(model), 
            hex_cells = F,
            plot_trajectory = T, 
            size_cells = 1, alpha_cells = 0.8) + 
  theme(aspect.ratio = 1)


mile_group <- data.frame(group_onto_nearest_milestones(model)) %>%
  set_colnames("milestone") %>%
  rownames_to_column("cell")


GSE164897$milestone <- mile_group$milestone



cd4_pseudo <- list()
for (i in 1:max(mile_group$milestone)) {
  cd4_pseudo[[i]] <- seurat_extract(GSE164897, meta1 = "milestone", value_meta1 = i)
}


# Load Hallmark Pathways
hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- format_pathways(hallmark_df)

# Run SCPA Analysis across Pseudotime Milestones
print("Running SCPA across milestones...")
pathActivation <- compare_pathways(samples = cd4_pseudo, 
                                   pathways = pathways, 
                                   parallel = TRUE, 
                                   cores = 4) # Reduced cores for safety


pathActivation <- as.data.frame(pathActivation)

pathActivation <- pathActivation %>%
  data.frame() %>%
  select(Pathway, qval) %>%
  column_to_rownames("Pathway")

# Handle Inf values for color scale (replace Inf with max finite value)
qvals <- as.matrix(pathActivation)
max_val <- max(qvals[is.finite(qvals)])
qvals[is.infinite(qvals)] <- max_val

col_hm <- colorRamp2(colors = c("white", "red"), breaks = c(0, max_val))

Heatmap(t(qvals),
        name = "Qvalue",
        col = col_hm,
        border = TRUE,
        rect_gp = gpar(col = "white", lwd = 0.1),
        show_column_dend = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 8))



plot_rank(scpa_out = pathActivation, 
          pathway = "HALLMARK_FATTY_ACID_METABOLISM",
          base_point_size = 2.5,
          highlight_point_size = 3)


# === Visualize Pathway Activity Across Clusters ===

print("Calculating per-cell pathway scores...")

# Filter pathways to only include genes present in the Seurat object
seurat_features <- rownames(GSE164897)

# Check if pathways is a data frame (SCPA format) or list
if (is.data.frame(pathways)) {
  pathways_list <- split(pathways$Genes, pathways$Pathway)
} else {
  pathways_list <- pathways
}

# Filter to genes present in dataset
pathways_filtered <- lapply(pathways_list, function(genes) {
  g_vec <- if(is.data.frame(genes) || is.list(genes)) unlist(genes) else genes
  base::intersect(as.character(g_vec), as.character(seurat_features))
})

# Remove empty pathways
pathways_filtered <- pathways_filtered[sapply(pathways_filtered, length) > 0]
print(paste("Retained", length(pathways_filtered), "pathways with overlapping genes"))

# Calculate pathway scores
GSE164897 <- AddModuleScore(GSE164897, features = pathways_filtered, name = "PathwayScore")

# Rename columns to pathway names
meta <- GSE164897@meta.data
renamed_count <- 0
for(i in seq_along(pathways_filtered)) {
  old_name <- paste0("PathwayScore", i)
  new_name <- names(pathways_filtered)[i]
  
  # Check both old name exists and new name is valid
  if(old_name %in% colnames(meta) && !is.null(new_name) && nchar(new_name) > 0) {
    colnames(meta)[colnames(meta) == old_name] <- new_name
    renamed_count <- renamed_count + 1
  }
}
GSE164897@meta.data <- meta
print(paste("Successfully renamed", renamed_count, "pathway columns"))

# Find top 5 most variable pathways across clusters
print("Identifying top variable pathways across clusters...")
pathway_names <- names(pathways_filtered)

# Verify pathway names exist in metadata
pathway_names_in_meta <- pathway_names[pathway_names %in% colnames(GSE164897@meta.data)]
print(paste("Found", length(pathway_names_in_meta), "pathway columns in metadata"))

if(length(pathway_names_in_meta) == 0) {
  print("ERROR: No pathway columns found in metadata. Using PathwayScore columns instead.")
  pathway_names_in_meta <- grep("^PathwayScore", colnames(GSE164897@meta.data), value = TRUE)
}

cluster_means <- GSE164897@meta.data %>%
  select(seurat_clusters, all_of(pathway_names_in_meta)) %>%
  group_by(seurat_clusters) %>%
  summarise(across(everything(), mean))

pathway_vars <- cluster_means %>%
  select(-seurat_clusters) %>%
  summarise(across(everything(), var)) %>%
  t() %>%
  as.data.frame() %>%
  setNames("variance") %>%
  arrange(desc(variance))

top_pathways <- rownames(pathway_vars)
# Remove any NA values
top_pathways <- top_pathways[!is.na(top_pathways)]
print(paste("Plotting", length(top_pathways), "pathways in groups of 5..."))

# Create FeaturePlots for all pathways, 5 per combined plot
if(length(top_pathways) > 0) {
  print("Creating FeaturePlots...")
  
  # Split into groups of 5
  pathway_groups <- split(top_pathways, ceiling(seq_along(top_pathways) / 5))
  
  # Create a plot for each group
  for(group_idx in seq_along(pathway_groups)) {
    group_pathways <- pathway_groups[[group_idx]]
    
    # Create a mapping from PathwayScore names to actual pathway names
    pathway_names_for_plot <- character(length(group_pathways))
    
    for(i in seq_along(group_pathways)) {
      score_col <- group_pathways[i]
      
      if(grepl("^PathwayScore", score_col)) {
        # Extract the index number
        idx <- as.numeric(gsub("PathwayScore", "", score_col))
        # Get the actual pathway name from the pathways list
        if(!is.na(idx) && idx <= length(pathways)) {
          # pathways is a list of tibbles, extract the Pathway name from the first row
          pathway_name <- pathways[[idx]]$Pathway[1]
          # Check if pathway name is valid (not NULL, not empty)
          if(!is.null(pathway_name) && length(pathway_name) > 0 && !is.na(pathway_name) && nchar(pathway_name) > 0) {
            pathway_names_for_plot[i] <- pathway_name
          } else {
            pathway_names_for_plot[i] <- score_col  # Fallback to score name
          }
        } else {
          pathway_names_for_plot[i] <- score_col
        }
      } else {
        pathway_names_for_plot[i] <- score_col
      }
    }
    
    plot_list <- lapply(1:length(group_pathways), function(i) {
      score_col <- group_pathways[i]
      pathway_name <- pathway_names_for_plot[i]
      
      FeaturePlot(GSE164897, features = score_col, reduction = "umap") + 
        scale_color_viridis_c(option = "C") + 
        ggtitle(gsub("HALLMARK_", "", pathway_name)) +
        theme(aspect.ratio = 1, 
              plot.title = element_text(size = 9))
    })
    
    # Combine plots for this group
    library(patchwork)
    combined_plot <- wrap_plots(plot_list, ncol = 3)
    print(combined_plot)
    print(paste("Plotted group", group_idx, "of", length(pathway_groups)))
  }
} else {
  print("No valid pathways found for plotting")
}
