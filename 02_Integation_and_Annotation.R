suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(readr)
  library(pheatmap)
})


dir.create("results/integration", recursive = TRUE, showWarnings = FALSE)
dir.create("results/annotation",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/markers",     recursive = TRUE, showWarnings = FALSE)

qc_files <- list.files("data/qc", pattern = "\\.qc\\.rds$", full.names = TRUE)
stopifnot(length(qc_files) > 0)

read_qc <- function(fp){
  x <- readRDS(fp)
  if (is.null(x$sample_id)) x$sample_id <- sub("\\.qc\\.rds$", "", basename(fp))
  x
}
sc.list <- setNames(lapply(qc_files, read_qc),
                    gsub("\\.qc\\.rds$", "", basename(qc_files)))

## --- 1) Integration: RPCA ---
sc.list <- lapply(sc.list, function(o){
  o <- NormalizeData(o, verbose = FALSE)
  o <- FindVariableFeatures(o, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  o
})
features <- SelectIntegrationFeatures(object.list = sc.list, nfeatures = 3000)
sc.list <- lapply(sc.list, function(o){
  o <- ScaleData(o, features = features, verbose = FALSE)
  o <- RunPCA(o, features = features, npcs = 50, verbose = FALSE)
  o
})

anchors <- FindIntegrationAnchors(
  object.list     = sc.list,
  anchor.features = features,
  reduction       = "rpca",
  k.anchor        = 20,
  dims            = 1:50
)
integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:50)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:50)
integrated <- FindClusters(integrated, resolution = 1.6)

ggsave("results/integration/UMAP_clusters.pdf",
       DimPlot(integrated, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend(),
       width = 8, height = 6, limitsize = FALSE)

## --- 2) CellTypist  annotation---
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated, verbose = FALSE)
SaveH5Seurat(integrated, filename = "results/integration/tmp.h5seurat", overwrite = TRUE)
Convert("results/integration/tmp.h5seurat", "results/integration/tmp.h5ad", overwrite = TRUE)

py <- "
import scanpy as sc, celltypist, pandas as pd
ad = sc.read_h5ad('results/integration/tmp.h5ad')
try:
    celltypist.models.download_models('Immune_All_Low')
except Exception:
    pass
res = celltypist.annotate(ad, model='Immune_All_Low.pkl', majority_voting=True)
obs = res.to_adata().obs.reset_index().rename(columns={'index':'cell'})
obs.to_csv('results/annotation/celltypist_labels.csv', index=False)
"
tf <- tempfile(fileext = ".py"); writeLines(py, tf); system2("python", tf)

ct <- suppressMessages(read_csv("results/annotation/celltypist_labels.csv", show_col_types = FALSE))
integrated$ct_fine <- ct$predicted_labels[match(colnames(integrated), ct$cell)]
if ("majority_voting" %in% names(ct)) {
  integrated$ct_conf <- ct$majority_voting[match(colnames(integrated), ct$cell)]
}

## --- 3) FindAllMarkers ---
DefaultAssay(integrated) <- "RNA"
markers <- FindAllMarkers(
  integrated,
  only.pos        = TRUE,
  test.use        = "wilcox",
  min.pct         = 0.25,
  logfc.threshold = 0.25
) |>
  arrange(cluster, desc(avg_log2FC))
write.csv(markers, "results/markers/FindAllMarkers_all.csv", row.names = FALSE)
top5 <- markers |>
  group_by(cluster) |>
  slice_max(order_by = avg_log2FC, n = 5) |>
  ungroup()
write.csv(top5, "results/markers/top5_per_cluster.csv", row.names = FALSE)

## Create a suggestion table for manual review (CellTypist majority + FindAllMarkers genes)
ct_major <- integrated@meta.data |>
  mutate(cluster = as.character(seurat_clusters)) |>
  count(cluster, ct_fine, name = "n") |>
  group_by(cluster) |>
  mutate(frac = n/sum(n)) |>
  slice_max(frac, n = 1, with_ties = FALSE) |>
  ungroup() |>
  transmute(cluster, ct_majority = ct_fine, ct_frac = frac)

markers_per_cluster <- markers |>
  group_by(cluster) |>
  summarise(markers = paste(gene, collapse = "; "), .groups = "drop")

suggested <- left_join(ct_major, markers_per_cluster, by = "cluster")
write.csv(suggested, "results/annotation/suggested_labels.csv", row.names = FALSE)

## --- 4) Finalize fine/major labels ---
## 读取可选人工修订文件 results/annotation/manual_override.csv（两列：cluster,final_fine）
final_map <- suggested |> transmute(cluster, final_fine = ct_majority)
if (file.exists("results/annotation/manual_override.csv")) {
  over <- suppressMessages(read_csv("results/annotation/manual_override.csv", show_col_types = FALSE))
  stopifnot(all(c("cluster","final_fine") %in% names(over)))
  over$cluster <- as.character(over$cluster)
  final_map <- final_map |> select(cluster) |> left_join(over, by = "cluster")
}

integrated$final_fine <- final_map$final_fine[match(as.character(integrated$seurat_clusters), final_map$cluster)]

fine_to_major <- c(
  "CD8 NKT"="T","Cytotoxic CD8 T"="T","DN T"="T","DP T"="T","DP T cycling"="T",
  "Exhausted CD4 T"="T","Memory CD4 T"="T","Naïve CD4 T"="T","Naïve CD8 T"="T",
  "Pre T"="T","Proliferating T"="T",
  "Erythroid"="Myeloid","Kupffer"="Myeloid","Monocyte-Macrophage"="Myeloid","Mast"="Myeloid",
  "mDC"="Myeloid","Monocytes"="Myeloid","Myelocyte"="Myeloid","NK"="Myeloid","pDC"="Myeloid",
  "Pre-Monocytes"="Myeloid","Type_1 Macrophage"="Myeloid","Type_2 Macrophage"="Myeloid",
  "Activated Follicular B"="B","Follicular B"="B","Germinal center B"="B",
  "Marginal zone B"="B","Memory B"="B","Naïve B"="B","Plasma B"="B",
  "Ciliated"="Endothelial","Endothelial"="Endothelial",
  "Fibroblasts"="Fibroblasts",
  "Platelets"="Platelets",
  "Granulocyte–Monocyte Progenitor"="HSC/MPPs",
  "Hepatocytes"="Parenchyma","SLC5A10+ Proximal tubule"="Parenchyma",
  "Surface mucous"="Parenchyma","Ventricular cardiomyocytes"="Parenchyma"
)
integrated$final_major <- fine_to_major[integrated$final_fine]

major_panel <- list(
  "B"           = c("MS4A1","CD79A","CD79B","CD22","MZB1","XBP1"),
  "T"           = c("CD3D","CD3E","TRAC","IL7R","CCR7","CD4","CD8A","NKG7","GNLY","PRF1","GZMB"),
  "Myeloid"     = c("LYZ","LST1","S100A8","S100A9","FCGR3A","MS4A7","MRC1","MARCO","IL1B"),
  "Endothelial" = c("PECAM1","VWF","KDR","CLDN5"),
  "Fibroblasts" = c("COL1A1","COL1A2","COL3A1","DCN"),
  "Platelets"   = c("PPBP","PF4","ITGA2B"),
  "HSC/MPPs"    = c("GATA2","KIT","CD34"),
  "Parenchyma"  = c("ALB","EPCAM","KRT8","KRT18","SLC5A10")
)
major_panel <- lapply(major_panel, function(gs) intersect(gs, rownames(integrated)))
major_panel <- major_panel[sapply(major_panel, length) > 0]

tmp <- AddModuleScore(integrated, features = unname(major_panel), name = "major_", assay = "RNA")
score_cols <- tail(colnames(tmp@meta.data), length(major_panel))
colnames(tmp@meta.data)[match(score_cols, colnames(tmp@meta.data))] <- paste0("score_major_", names(major_panel))
integrated <- tmp

missing_major <- is.na(integrated$final_major)
if (any(missing_major)) {
  score_mat <- integrated@meta.data[missing_major, paste0("score_major_", names(major_panel)), drop = FALSE]
  best <- names(major_panel)[max.col(as.matrix(score_mat), ties.method = "first")]
  integrated$final_major[missing_major] <- best
}
integrated$final_major[is.na(integrated$final_major)] <- "Other"

## --- 5) Plot Figures ---
pal_major <- c(
  "B"="#1F77B4","Endothelial"="#2CA02C","Fibroblasts"="#17BECF",
  "HSC/MPPs"="#E377C2","Myeloid"="#FF7F0E","Parenchyma"="#8C564B",
  "Platelets"="#D62728","T"="#9467BD","Other"="#7F7F7F"
)
integrated$final_major <- factor(integrated$final_major,
                                 levels = c("T","B","Myeloid","Endothelial","Fibroblasts","Parenchyma","Platelets","HSC/MPPs","Other"))

ggsave("results/integration/UMAP_major_with_legend.pdf",
       DimPlot(integrated, reduction = "umap", group.by = "final_major") +
         scale_color_manual(values = pal_major, drop = FALSE),
       width = 8, height = 6, limitsize = FALSE)

fine_levels_plot <- c(
  "CD8 NKT","Cytotoxic CD8 T","DN T","DP T","DP T cycling",
  "Exhausted CD4 T","Memory CD4 T","Naïve CD4 T","Naïve CD8 T","Pre T","Proliferating T",
  "Erythroid","Kupffer","Monocyte-Macrophage","Mast","mDC","Monocytes","Myelocyte","NK","pDC","Pre-Monocytes","Type_1 Macrophage","Type_2 Macrophage",
  "Activated Follicular B","Follicular B","Germinal center B","Marginal zone B","Memory B","Naïve B","Plasma B",
  "Ciliated","Endothelial","Fibroblasts","Platelets","Granulocyte–Monocyte Progenitor",
  "Hepatocytes","SLC5A10+ Proximal tubule","Surface mucous","Ventricular cardiomyocytes",
  "Other"
)
integrated$final_fine <- factor(integrated$final_fine, levels = fine_levels_plot)

ggsave("results/integration/UMAP_fine_with_legend.pdf",
       DimPlot(integrated, reduction = "umap", group.by = "final_fine"),
       width = 10, height = 7, limitsize = FALSE)

if (!"source" %in% colnames(integrated@meta.data)) integrated$source <- integrated$sample_id
comp <- integrated@meta.data |>
  count(source, final_major, name = "n") |>
  group_by(source) |>
  mutate(pct = 100 * n / sum(n)) |>
  ungroup()
ggsave("results/integration/composition_by_tissue.pdf",
       ggplot(comp, aes(x = source, y = pct, fill = final_major)) +
         geom_col(width = 0.8) + coord_flip() +
         scale_fill_manual(values = pal_major, drop = FALSE) +
         labs(x = NULL, y = "Percent", fill = "Cluster") +
         theme_classic(base_size = 11),
       width = 7, height = 5, limitsize = FALSE)

genes_panel <- unique(unlist(major_panel))
ggsave("results/integration/dotplot_major_panel.pdf",
       DotPlot(integrated, features = genes_panel, group.by = "final_major",
               cols = c("lightgrey","#D55E00")) + RotatedAxis() + theme_classic(base_size = 11),
       width = 9, height = 6, limitsize = FALSE)

Idents(integrated) <- integrated$seurat_clusters
topN <- markers |>
  group_by(cluster) |>
  slice_max(order_by = avg_log2FC, n = 8) |>
  pull(gene) |>
  unique()
topN <- intersect(topN, rownames(integrated))

DefaultAssay(integrated) <- "RNA"
integrated <- ScaleData(integrated, features = topN, verbose = FALSE)
avg_mat <- AverageExpression(integrated, features = topN, group.by = "seurat_clusters",
                             assays = "RNA", slot = "scale.data")$RNA
cl2maj <- integrated@meta.data |>
  group_by(seurat_clusters) |>
  summarise(major = names(sort(table(final_major), decreasing = TRUE))[1], .groups = "drop")
maj_vec <- cl2maj$major[match(colnames(avg_mat), cl2maj$seurat_clusters)]
ord <- order(match(maj_vec, levels(integrated$final_major)), as.numeric(colnames(avg_mat)))
avg_mat <- avg_mat[, ord, drop = FALSE]
ann_col <- data.frame(Major = maj_vec[ord]); rownames(ann_col) <- colnames(avg_mat)
pheatmap(avg_mat,
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = TRUE, show_colnames = TRUE,
         annotation_col = ann_col, annotation_colors = list(Major = pal_major),
         filename = "results/integration/marker_heatmap_topN.pdf",
         width = 9, height = 12)

saveRDS(integrated, "results/integration/integrated_seurat.rds")
