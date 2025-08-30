suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(purrr)
  library(ggplot2)
})

dir.create("results/tonsil", recursive = TRUE, showWarnings = FALSE)

## Collect tonsil QC objects
qc_files <- list.files("data/qc", pattern = "(?i)(tonsil|palatine|pharyngeal).*\\.qc\\.rds$", full.names = TRUE)
stopifnot(length(qc_files) > 0)
tonsil.list <- setNames(lapply(qc_files, readRDS),
                        gsub("\\.qc\\.rds$", "", basename(qc_files)))

## LogNormalize
tonsil.list <- lapply(tonsil.list, function(x){
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  x
})

## CCA anchors and integration
features <- SelectIntegrationFeatures(tonsil.list, nfeatures = 3000)
tonsil.list <- lapply(tonsil.list, function(x){
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, npcs = 20, verbose = FALSE)
  x
})
anchors <- FindIntegrationAnchors(object.list = tonsil.list,
                                  anchor.features = features,
                                  reduction = "cca",
                                  dims = 1:20)
tonsil <- IntegrateData(anchorset = anchors,
                        dims = 1:20,
                        normalization.method = "LogNormalize")

DefaultAssay(tonsil) <- "integrated"
tonsil <- ScaleData(tonsil, verbose = FALSE)
tonsil <- RunPCA(tonsil, npcs = 20, verbose = FALSE)
tonsil <- RunUMAP(tonsil, reduction = "pca", dims = 1:20)
tonsil <- FindNeighbors(tonsil, reduction = "pca", dims = 1:20)
tonsil <- FindClusters(tonsil, resolution = 0.6)

p <- DimPlot(tonsil, reduction = "umap", group.by = group_col) +
  theme(aspect.ratio = 1)

ggsave("results/tonsil/UMAP_Tonsil_B.pdf", p, width = 5.5, height = 5.0)
saveRDS(tonsil, "results/tonsil/tonsil.rds")
