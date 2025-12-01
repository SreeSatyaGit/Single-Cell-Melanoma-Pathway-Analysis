library(SCPA)
library(Seurat)
library(msigdbr)
library(magrittr)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)


hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")


c2_df <- msigdbr(species = "Homo sapiens", category = "C2")



target_subcats <- c("CP:KEGG", "CP:REACTOME", "CP:BIOCARTA", "CP:WIKIPATHWAYS", "CP:PID")
c2_filtered <- c2_df %>% 
  dplyr::filter(gs_subcat %in% target_subcats)

all_pathways_df <- bind_rows(hallmark_df, c2_filtered)


pathways <- format_pathways(all_pathways_df)


split_treatment <- SplitObject(GSE164897, split.by = "treatment")



# Initialize lists
Pvem_tram <- list(); Pvem_cob <- list(); Pvem <- list()

clusters <- unique(GSE164897$seurat_clusters)
print(paste("Processing", length(clusters), "clusters..."))

for (i in clusters) {
  print(paste("Analyzing cluster:", i))
  
  cluster_obj <- subset(GSE164897, seurat_clusters == i)
  

  if (ncol(cluster_obj) < 20) {
    print(paste("Skipping cluster", i, "- too few cells"))
    next
  }
  
  tryCatch({

    untreated_sub <- seurat_extract(cluster_obj, meta1 = "treatment", value_meta1 = "untreated")
    vem_sub <- seurat_extract(cluster_obj, meta1 = "treatment", value_meta1 = "Vemurafenib")
    vem_cob_sub <- seurat_extract(cluster_obj, meta1 = "treatment", value_meta1 = "vem_cob")
    vem_tram_sub <- seurat_extract(cluster_obj, meta1 = "treatment", value_meta1 = "vem_tram")
    

    idx <- as.character(i)

    Pvem_tram[[idx]] <- compare_pathways(list(untreated_sub, vem_tram_sub), pathways, parallel = TRUE, cores = 4)
    Pvem_cob[[idx]] <- compare_pathways(list(untreated_sub, vem_cob_sub), pathways, parallel = TRUE, cores = 4) 
    Pvem[[idx]] <- compare_pathways(list(untreated_sub, vem_sub), pathways, parallel = TRUE, cores = 4)
    
  }, error = function(e) {
    print(paste("Error in cluster", i, ":", e$message))
  })
}


# Helper function to get Signed Q-values (sign(FC) * -log10(qval))
get_signed_scores <- function(scpa_out, name) {
  df <- list()
  for (i in names(scpa_out)) {
    res <- scpa_out[[i]]
    # Calculate signed score: Positive = Upregulated, Negative = Downregulated
    # Add small epsilon to qval to avoid Inf
    res$score <- sign(res$FC) * -log10(res$qval + 1e-50)
    
    df[[i]] <- res %>% dplyr::select(Pathway, score)
  }
  
  col_names <- names(df)
  for (i in 1:length(df)) {
    df[[i]] <- setNames(df[[i]], c("pathway", paste(name, col_names[[i]], sep = "_")))
  }
  return(df)
}


scpa_results <- Reduce(full_join, c(get_signed_scores(Pvem_tram, "Vemtram"),
                                    get_signed_scores(Pvem_cob, "Vemcob"),
                                    get_signed_scores(Pvem, "Vem")))


all_data <- scpa_results %>%
  column_to_rownames("pathway")


all_data <- all_data[grep("HALLMARK_", rownames(all_data)), ]


all_data[is.na(all_data)] <- 0
mat <- as.matrix(all_data)


rownames(mat) <- gsub("HALLMARK_", "", rownames(mat))
rownames(mat) <- gsub("REACTOME_", "", rownames(mat))
rownames(mat) <- gsub("KEGG_", "", rownames(mat))
rownames(mat) <- gsub("BIOCARTA_", "", rownames(mat))
rownames(mat) <- gsub("WP_", "", rownames(mat))
rownames(mat) <- gsub("_", " ", rownames(mat)) 


row_var <- apply(mat, 1, var)
top_pathways <- names(sort(row_var, decreasing = TRUE))[1:20]
mat_subset <- mat[top_pathways, , drop = FALSE]


col_names <- colnames(mat)
treatment <- sapply(strsplit(col_names, "_"), `[`, 1)
cluster <- sapply(strsplit(col_names, "_"), `[`, 2)


n_clusters <- length(unique(cluster))
cluster_colors <- setNames(rainbow(n_clusters), unique(cluster))
treatment_colors <- c("Vem" = "#E74C3C", "Vemcob" = "#3498DB", "Vemtram" = "#2ECC71")

col_an <- HeatmapAnnotation(
  Treatment = treatment,
  Cluster = cluster,
  col = list(Treatment = treatment_colors, Cluster = cluster_colors),
  gp = gpar(col = "white", lwd = 0.05),
  simple_anno_size = unit(3, "mm")
)

print(paste("Plotting signed heatmap for top", nrow(mat_subset), "pathways"))


Heatmap(mat_subset,
        name = "Signed -log10(qval)",
        top_annotation = col_an,
        col = colorRamp2(c(-10, 0, 10), c("blue", "white", "red")), 
        show_row_names = TRUE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 9),
        column_names_gp = gpar(fontsize = 8),
        rect_gp = gpar(col = "white", lwd = 0.5))

apply(mat, 1, var) %>%
  as.data.frame() %>% 
  setNames("variance") %>%
  arrange(desc(variance)) %>% 
  rownames_to_column("pathway") %>%
  ggplot(aes(reorder(pathway, variance), variance)) +
  geom_point(shape = 21, cex = 3, fill = "royalblue2", color = 'black', stroke = 0.2) +
  scale_x_discrete(expand = c(0.04, 0.04)) +
  labs(title = "Pathway Variance Across Conditions",
       x = "Pathway", 
       y = "Variance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"))





cluster_id <- "0"  
if(cluster_id %in% names(Pvem_tram)) {
  plot_data <- Pvem_tram[[cluster_id]] %>%
    mutate(color = case_when(FC > 5 & adjPval < 0.01 ~ '#6dbf88',
                             FC < 5 & FC > -5 & adjPval < 0.01 ~ '#84b0f0',
                             FC < -5 & adjPval < 0.01 ~ 'mediumseagreen',
                             TRUE ~ 'black')) # Added TRUE catch-all
  
  aa_path <- plot_data %>% 
    filter(grepl(pattern = "reactome_arachi", ignore.case = TRUE, x = Pathway))
  
  print(paste("Plotting Volcano for Cluster", cluster_id))
  
  top_labels <- plot_data %>% 
    arrange(desc(abs(FC))) %>% 
    head(5)
    
  label_points <- bind_rows(top_labels, aa_path) %>% distinct()
  
  # Ensure ggrepel is available
  if (!require("ggrepel", quietly = TRUE)) install.packages("ggrepel")
  library(ggrepel)

  ggplot(plot_data, aes(-FC, qval)) +
    geom_vline(xintercept = c(-5, 5), linetype = "dashed", col = 'black', lwd = 0.3) +
    geom_point(cex = 2.6, shape = 21, fill = plot_data$color, stroke = 0.3) +
    geom_point(data = aa_path, shape = 21, cex = 2.8, fill = "orangered2", color = "black", stroke = 0.3) +
    # Add labels
    geom_text_repel(data = label_points, aes(label = Pathway), size = 3, max.overlaps = 20, box.padding = 0.5) +
    # Dynamic limits to ensure all points are shown
    scale_x_continuous(expand = expansion(mult = 0.1)) +
    scale_y_continuous(expand = expansion(mult = 0.1)) +
    xlab("Enrichment") +
    ylab("Qval") +
    ggtitle(paste("Volcano Plot: Cluster", cluster_id)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          aspect.ratio = 1)
} else {
  print(paste("Cluster", cluster_id, "not found in results"))
}
