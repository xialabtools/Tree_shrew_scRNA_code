suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(harmony)
  library(patchwork)
})

dir.create("results/macrophage", recursive = TRUE, showWarnings = FALSE)

## Load the integrated atlas##
atlas <- readRDS("results/integration/integrated_seurat.rds")

## Panel A: CSF1R feature UMAP on the full atlas ##
DefaultAssay(atlas) <- if ("integrated" %in% Assays(atlas)) "integrated" else DefaultAssay(atlas)
pA <- FeaturePlot(
  atlas,
  features = "CSF1R",
  reduction = "umap",
  order = TRUE,
  min.cutoff = "q5",
  max.cutoff = "q95",
  cols = c("grey90", "#5b2a86")
) + theme(aspect.ratio = 1) + labs(title = "CSF1R")

ggsave("results/macrophage/panel_A_CSF1R_feature_umap.pdf", pA, width = 5.5, height = 5)


label_col <- dplyr::case_when(
  "final_fine" %in% colnames(atlas@meta.data) ~ "final_fine",
  "ct_fine"    %in% colnames(atlas@meta.data) ~ "ct_fine",
  TRUE ~ "seurat_clusters"
)

macro_labels <- c("Kupffer","Type_1 Macrophage","Type_2 Macrophage")
mac <- subset(atlas, subset = get(label_col) %in% macro_labels)

## Tissue column (source) is required
stopifnot("source" %in% colnames(mac@meta.data))


## Filter tissues with macrophage cell count < 100 ##
keep_src <- names(which(table(mac$source) >= 100))
mac <- subset(mac, subset = source %in% keep_src)

src_levels <- c("Heart","Liver","Lung","Lymph","Marrow","PBMC","Spleen","Thymus","Tonsil")
mac$source <- factor(as.character(mac$source),
                     levels = intersect(src_levels, unique(as.character(mac$source))))

grp_var <- dplyr::case_when(
  "sample_id" %in% colnames(mac@meta.data) ~ "sample_id",
  "batch"     %in% colnames(mac@meta.data) ~ "batch",
  "source"    %in% colnames(mac@meta.data) ~ "source",
  TRUE ~ "orig.ident"
)

mac <- SCTransform(mac, verbose = FALSE)
mac <- RunPCA(mac, npcs = 50, verbose = FALSE)
mac <- RunHarmony(mac, group.by.vars = grp_var, assay.use = "SCT", max.iter.harmony = 20)

dims_use <- 1:15
mac <- RunUMAP(mac, reduction = "harmony", dims = dims_use, verbose = FALSE)
mac <- FindNeighbors(mac, reduction = "harmony", dims = dims_use, verbose = FALSE)
mac <- FindClusters(mac, resolution = 0.1, verbose = FALSE)

## Panel C: UMAP of macrophages colored by tissue (source) ##
pC <- DimPlot(
  mac,
  reduction = "umap",
  group.by = "source",
  label = FALSE
) + theme(aspect.ratio = 1) + labs(title = "Macrophages")

ggsave("results/macrophage/panel_C_macrophage_umap_by_tissue.pdf", pC, width = 7.2, height = 4.8)

saveRDS(mac, file = "results/macrophage/macrophage_SCT_harmony.rds")