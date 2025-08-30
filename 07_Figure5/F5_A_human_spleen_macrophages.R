suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

dir.create("results/cross_species", recursive = TRUE, showWarnings = FALSE)
human_file <- "data/cross_species/human_spleen_myeloid.rds"   ## <- edit path if needed
out_pdf    <- "results/cross_species/panel_A_human_myeloid_NR1H3.pdf"

choose_label_col <- function(obj, prefs = c("Manually_curated_celltype",
                                            "Majority_voting_CellTypist",
                                            "Predicted_labels_CellTypist",
                                            "final_fine","celltype","seurat_clusters")) {
  ok <- prefs[prefs %in% colnames(obj@meta.data)]
  if (length(ok) == 0) "seurat_clusters" else ok[1]
}
find_gene <- function(obj, symbol){
  if (symbol %in% rownames(obj)) return(symbol)
  hits <- grep(paste0("^", symbol, "$"), rownames(obj), ignore.case = TRUE, value = TRUE)
  if (length(hits)) hits[1] else stop(symbol, " not found.")
}

human <- readRDS(human_file)
DefaultAssay(human) <- "RNA"

organ_col <- intersect(c("Organ","source","tissue","Tissue","organ"), colnames(human@meta.data))
if (length(organ_col)) {
  organ_col <- organ_col[1]
  is_spleen <- human@meta.data[[organ_col]] %in% c("SPL","Spleen","spleen")
  if (any(is_spleen)) human <- subset(human, cells = rownames(human@meta.data)[is_spleen])
}

if (is.null(human@reductions$umap)) {
  human <- NormalizeData(human, verbose = FALSE) |>
           FindVariableFeatures(nfeatures = 3000, verbose = FALSE) |>
           ScaleData(verbose = FALSE) |>
           RunPCA(npcs = 30, verbose = FALSE) |>
           FindNeighbors(dims = 1:30) |>
           RunUMAP(dims = 1:30)
}

## Left: cluster UMAP
label_col <- choose_label_col(human)
p_left <- DimPlot(human, reduction = "umap", group.by = label_col,
                  label = TRUE, repel = TRUE) +
  labs(title = "Human spleen myeloid") +
  theme(aspect.ratio = 1)

## Right: NR1H3 feature plot
gene_to_plot <- tryCatch(find_gene(human, "NR1H3"),
                         error = function(e) find_gene(human, "CSF1R"))
p_right <- FeaturePlot(human, reduction = "umap", features = gene_to_plot,
                       min.cutoff = "q5", max.cutoff = "q95", order = TRUE,
                       cols = c("grey90", "#5B2A86")) +
  labs(title = "NR1H3") + theme(aspect.ratio = 1)

## Save
(p_left | p_right)
ggsave(out_pdf, width = 9.5, height = 4.6)
