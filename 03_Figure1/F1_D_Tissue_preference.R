## Dependencies
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(epitools)
  library(ComplexHeatmap)
  library(grid)
  library(circlize)
})

## Load integrated object (created by integrate_annotate_markers_final.R)
if (!exists("integrated")) {
  integrated <- readRDS("results/integration/integrated_seurat.rds")
}

## Sanity checks
stopifnot(all(c("source","final_fine","final_major") %in% colnames(integrated@meta.data)))

## Optional: normalize tissue labels to match the figure
normalize_tissue <- function(x) {
  x <- sub("^Nasopharynx$", "Pharyngeal_tonsil", x)
  x <- sub("^Tonsil$",      "Palatine_tonsil",  x)
  x
}
integrated$source2 <- normalize_tissue(as.character(integrated$source))

## Observed × Expected and R_{o/e}
obs <- as.matrix(table(integrated$source2, integrated$final_fine))
exp <- epitools::expected(obs)
exp[exp == 0] <- NA
roe <- obs / exp

## Row (tissue) order
row_order <- c("Marrow","Thymus","PBMC","Lymph","Spleen","Pharyngeal_tonsil","Palatine_tonsil")
row_order <- row_order[row_order %in% rownames(roe)]

## Column order
col_order <- c(
  "Naïve B","Germinal center B","Marginal zone B","Memory B","Follicular B",
  "Activated Follicular B","Plasma B",
  "DN T","DP T","DP T cycling","Pre T","Naïve CD4 T","Memory CD4 T","Exhausted CD4 T",
  "Naïve CD8 T","Cytotoxic CD8 T","CD8 NKT","Proliferating T",
  "Erythroid","Myelocyte","Platelets","Fibroblasts","Granulocyte–Monocyte Progenitor",
  "Monocyte-Macrophage","Mast","Pre-Monocytes","Monocytes","mDC","pDC","NK",
  "Type_1 Macrophage","Type_2 Macrophage"
)
col_order <- col_order[col_order %in% colnames(roe)]

## Subset and cap values for display
roe <- roe[row_order, col_order, drop = FALSE]
roe_cap <- pmin(roe, 5)

## Discrete fill palette and thresholds
fill_cols <- c("1" = "#FEE6CE", "1.5" = "#FDC08C", "3" = "#F5904B", ">5" = "#E6550D")

## Plus sign text by thresholds
plus_txt <- matrix("", nrow = nrow(roe), ncol = ncol(roe), dimnames = dimnames(roe))
plus_txt[roe >= 1   & roe < 1.5] <- "+"
plus_txt[roe >= 1.5 & roe < 3  ] <- "++"
plus_txt[roe >= 3             ] <- "+++"

fine_to_major <- c(
  "CD8 NKT"="T","Cytotoxic CD8 T"="T","DN T"="T","DP T"="T","DP T cycling"="T",
  "Exhausted CD4 T"="T","Memory CD4 T"="T","Naïve CD4 T"="T","Naïve CD8 T"="T","Pre T"="T","Proliferating T"="T",
  "Erythroid"="Myeloid","Kupffer"="Myeloid","Monocyte-Macrophage"="Myeloid","Mast"="Myeloid","mDC"="Myeloid",
  "Monocytes"="Myeloid","Myelocyte"="Myeloid","NK"="Myeloid","pDC"="Myeloid","Pre-Monocytes"="Myeloid",
  "Type_1 Macrophage"="Myeloid","Type_2 Macrophage"="Myeloid",
  "Activated Follicular B"="B","Follicular B"="B","Germinal center B"="B","Marginal zone B"="B","Memory B"="B","Naïve B"="B","Plasma B"="B",
  "Ciliated"="Myeloid","Endothelial"="Myeloid","Fibroblasts"="Myeloid","Platelets"="Myeloid",
  "Granulocyte–Monocyte Progenitor"="Myeloid"
)
maj_for_cols <- unname(fine_to_major[colnames(roe)])
maj_for_cols[is.na(maj_for_cols)] <- "Myeloid"
maj_for_cols <- factor(maj_for_cols, levels = c("B","T","Myeloid"))

## Column split and block labels/colors (header stripes)
ha_top <- HeatmapAnnotation(
  blocks = anno_block(
    gp = gpar(fill = c("B" = "#1F77B4", "T" = "#C6CC65", "Myeloid" = "#8E44AD")),
    labels = c("B","T","Myeloid"),
    labels_gp = gpar(col = "white", fontsize = 10, fontface = "bold")
  ),
  which = "column"
)

## Custom cell drawer to reproduce your “+ / ++ / +++” overlay
cell_fun <- function(j, i, x, y, width, height, fill) {
  v <- roe[i, j]
  col <- if (is.na(v) || v < 1) {
    fill_cols[1]
  } else if (v < 1.5) {
    fill_cols[2]
  } else if (v < 3) {
    fill_cols[3]
  } else {
    fill_cols[4]
  }
  grid.rect(x = x, y = y, width = width, height = height,
            gp = gpar(col = "white", fill = col))
  lab <- plus_txt[i, j]
  if (nzchar(lab)) grid.text(lab, x = x, y = y, gp = gpar(fontsize = 8))
}

## Build heatmap
ht <- Heatmap(
  matrix = roe_cap,
  cell_fun = cell_fun,
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_heatmap_legend = FALSE,
  row_names_side = "left",
  column_split = maj_for_cols,
  top_annotation = ha_top,
  border = NA,
  name = "Ro/e"
)

## Legend matching your bins
lgd <- Legend(
  labels = names(fill_cols),
  title = expression(R[o/e]),
  legend_gp = gpar(fill = fill_cols)
)

## Save figure
dir.create("results/integration", showWarnings = FALSE, recursive = TRUE)
pdf("results/integration/cell_preference_heatmap.pdf", width = 10, height = 4)
draw(ht, annotation_legend_list = list(lgd))
dev.off()
