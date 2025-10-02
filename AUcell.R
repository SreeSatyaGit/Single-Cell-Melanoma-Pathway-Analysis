# AUCell analysis
library(AUCell)
library(ggplot2)

cells_rankings <- AUCell_buildRankings(E2, nCores = 4, plotStats = FALSE)
cells_AUC <- AUCell_calcAUC(pathways, cells_rankings)

# AUC matrix: pathways x cells
auc_matrix <- getAUC(cells_AUC)

# top active pathways globally
rowMeans(auc_matrix) |> sort(decreasing = TRUE) |> head(10)
pathway_means <- rowMeans(auc_matrix, na.rm = TRUE)

# sort and take top 10
pathway_means <- rowMeans(auc_matrix, na.rm = TRUE)
top10 <- sort(pathway_means, decreasing = TRUE)[1:10]

op <- par(mar = c(5, 14, 3, 1))  # wider left margin for labels
barplot(
  rev(as.numeric(top10)),
  names.arg = rev(names(top10)),
  horiz = TRUE,
  las = 1,
  xlab = "Mean AUCell score",
  col = "gray70",
  main = "Top 10 Highly Expressed Pathways (AUCell)"
)
par(op)

