library(Seurat)
library(msigdbr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(fgsea)
library(sctransform)
library(glmGamPoi)
# assume your Seurat object is called 'seu'
# Extract current gene names (Ensembl IDs, maybe with versions)
ens <- rownames(GSE164897)

# 1) Strip version suffix if present
ens_novers <- sub("\\..*$", "", ens)  # ENSG00000100867.10 → ENSG00000100867

# 2) Map ENSG → SYMBOL
map <- AnnotationDbi::select(org.Hs.eg.db,
                             keys = unique(ens_novers),
                             keytype = "ENSEMBL",
                             columns = "SYMBOL")

# Drop NAs and keep unique mappings
map <- map[!is.na(map$SYMBOL), c("ENSEMBL","SYMBOL")]
map <- map[!duplicated(map$ENSEMBL), ]   # remove duplicates

# 3) Replace rownames of Seurat object
# Match ENS IDs (no version) to SYMBOL
ens2sym <- setNames(map$SYMBOL, map$ENSEMBL)
newnames <- ens2sym[ens_novers]

# Only keep genes that successfully mapped
keep <- !is.na(newnames)
GSE164897 <- GSE164897[keep, ]
rownames(GSE164897) <- newnames[keep]



GSE164897 <- SCTransform(GSE164897, verbose = FALSE, variable.features.n = 10000)

GSE164897 <- RunPCA(GSE164897, assay = "SCT", verbose = FALSE,
              rev.pca = TRUE, reduction.name = "pca.rev",
              reduction.key="PCR_", npcs = 50)


E <- GSE164897@reductions$pca.rev@feature.loadings

pathwaysDF <- msigdbr(species = "Homo sapiens", collection = "C6")
pathways   <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)

set.seed(1)
gesecaRes <- geseca(pathways, E, minSize = 5, maxSize = 500, center = FALSE, eps=1e-100)

head(gesecaRes, 10)


library(data.table)
## Ensure data.table
setDT(gesecaRes)

## 1) Filter pathways starting with AKT_, RAF_, KRAS_, MTOR_ and take top 4
topPathways <- gesecaRes[grepl("^(AKT_|RAF_|KRAS_|MTOR_)", pathway)][1:10]

## 2) Character vector of selected pathway IDs (for subsetting the named list `pathways`)
sel <- topPathways[["pathway"]]

## Optional: guard against missing names in `pathways`
missing <- setdiff(sel, names(pathways))
if (length(missing)) warning("Missing in `pathways`: ", paste(missing, collapse = ", "))


GSE164897 <- RunUMAP(GSE164897, dims = 1:30, reduction = "integrated.cca")



ps <- plotCoregulationProfileReduction(pathways[sel], GSE164897,
                                       title=sel,
                                       reduction="umap")

cowplot::plot_grid(plotlist=ps[1:4], ncol=2)

cowplot::plot_grid(plotlist=ps[5:8], ncol=2)

cowplot::plot_grid(plotlist=ps[9:10], ncol=2)