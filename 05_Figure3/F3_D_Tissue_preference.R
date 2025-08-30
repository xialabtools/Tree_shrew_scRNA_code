suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(epitools)
  library(ComplexHeatmap)
  library(grid)
  library(circlize)
})

mac <- readRDS("results/macrophage/macrophage_SCT_harmony.rds")
stopifnot(all(c("source", "seurat_clusters") %in% colnames(mac@meta.data)))

normalize_tissue <- function(x) {
  x <- as.character(x)
  x[x == "Tonsil"] <- "Palatine tonsil"
  x
}
mac$source2 <- normalize_tissue(mac$source)

obs <- as.matrix(table(mac$source2, mac$seurat_clusters))
exp <- epitools::expected(obs)
exp[exp == 0] <- NA
roe <- obs / exp

roe_cap <- pmin(roe, 5)
plus_txt <- matrix("", nrow = nrow(roe), ncol = ncol(roe), dimnames = dimnames(roe))
plus_txt[roe >= 1   & roe < 1.5] <- "+"
plus_txt[roe >= 1.5 & roe < 3  ] <- "++"
plus_txt[roe >= 3             ] <- "+++"

row_order <- c("Heart","Liver","Lung","Lymph","Marrow","PBMC","Spleen","Thymus","Palatine tonsil")
row_order <- row_order[row_order %in% rownames(roe)]
roe_cap <- roe_cap[row_order, , drop = FALSE]
plus_txt <- plus_txt[row_order, , drop = FALSE]

col_order <- sort(as.numeric(colnames(roe_cap)))
col_order <- as.character(col_order)
roe_cap <- roe_cap[, col_order, drop = FALSE]
plus_txt <- plus_txt[, col_order, drop = FALSE]

row_totals <- rowSums(obs[row_order, , drop = FALSE], na.rm = TRUE)
rownames(roe_cap) <- paste0(row_order, " (", row_totals, ")")
rownames(plus_txt) <- rownames(roe_cap)

bin_cols <- c("1" = "#FEE6CE", "1.5" = "#FDC08C", "3" = "#F5904B", ">5" = "#E6550D")

cell_fun <- function(j, i, x, y, width, height, fill) {
  v <- roe_cap[i, j]
  col <- if (is.na(v) || v < 1) {
    bin_cols[1]
  } else if (v < 1.5) {
    bin_cols[2]
  } else if (v < 3) {
    bin_cols[3]
  } else {
    bin_cols[4]
  }
  grid.rect(x = x, y = y, width = width, height = height,
            gp = gpar(col = "white", fill = col))
  lab <- plus_txt[i, j]
  if (nzchar(lab)) grid.text(lab, x = x, y = y, gp = gpar(fontsize = 8))
}

lgd <- Legend(
  labels = names(bin_cols),
  title = expression(R[o/e]),
  legend_gp = gpar(fill = bin_cols)
)

top_anno <- HeatmapAnnotation(
  cluster = anno_text(colnames(roe_cap), gp = gpar(fontsize = 8)),
  annotation_height = unit(6, "mm")
)

dir.create("results/macrophage", showWarnings = FALSE, recursive = TRUE)

ht <- Heatmap(
  matrix = roe_cap,
  name = "Ro/e",
  cell_fun = cell_fun,
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE, show_column_names = FALSE,
  top_annotation = top_anno,
  border = NA
)

pdf("results/macrophage/panel_D_tissue_preference_heatmap.pdf", width = 6.5, height = 5.2)
draw(ht, annotation_legend_list = list(lgd))
dev.off()
