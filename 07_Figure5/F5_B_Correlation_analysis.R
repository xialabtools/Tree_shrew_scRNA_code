suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
})

dir.create("results/cross_species", recursive = TRUE, showWarnings = FALSE)
human_file     <- "data/cross_species/human_spleen_myeloid.rds"       ## <- edit path if needed
treeshrew_file <- "results/macrophage/macrophage_SCT_harmony.rds"     ## <- tree shrew macrophage object
out_pdf        <- "results/cross_species/panel_B_NR1H3pos_macrophage_correlation.pdf"

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
subset_NR1H3_macrophage_spleen <- function(obj, label_col,
                                           spleen_keys = c("SPL","Spleen","spleen"),
                                           macro_keys  = c("Macrophage","Macrophages","Monocyte-Macrophage",
                                                           "Type_1 Macrophage","Type_2 Macrophage","Kupffer")) {
  md <- obj@meta.data
  src_col <- intersect(c("Organ","source","tissue","Tissue","organ","orig.ident"), colnames(md))
  if (length(src_col)) {
    src_col <- src_col[1]
    keep <- md[[src_col]] %in% spleen_keys
    if (any(keep)) obj <- subset(obj, cells = rownames(md)[keep])
  }
  
  labs <- as.character(obj@meta.data[[label_col]])
  if (!is.null(labs)) {
    is_macro <- grepl(paste(macro_keys, collapse = "|"), labs, ignore.case = TRUE)
    if (any(is_macro)) obj <- subset(obj, cells = colnames(obj)[is_macro])
  }
  DefaultAssay(obj) <- "RNA"
  if (!nrow(GetAssayData(obj, slot = "data"))) obj <- NormalizeData(obj, verbose = FALSE)
  g <- find_gene(obj, "NR1H3")
  val <- FetchData(obj, vars = g)[,1]
  subset(obj, cells = colnames(obj)[val > 0])
}
avg_log_expr <- function(obj, genes = NULL) {
  DefaultAssay(obj) <- "RNA"
  if (!nrow(GetAssayData(obj, slot = "data"))) obj <- NormalizeData(obj, verbose = FALSE)
  m <- GetAssayData(obj, slot = "data")
  if (!is.null(genes)) m <- m[genes, , drop = FALSE]
  Matrix::rowMeans(m)
}

## Load objects
human <- readRDS(human_file)
ts    <- readRDS(treeshrew_file)
DefaultAssay(human) <- "RNA"; DefaultAssay(ts) <- "RNA"

## Subset NR1H3+ spleen macrophages for both species
label_h  <- choose_label_col(human)
label_ts <- choose_label_col(ts)
human_nr1h3 <- subset_NR1H3_macrophage_spleen(human, label_col = label_h)
ts_nr1h3    <- subset_NR1H3_macrophage_spleen(ts,    label_col = label_ts)

## Downsample up to 500 cells per species for balance
set.seed(123)
if (ncol(human_nr1h3) > 500) human_nr1h3 <- subset(human_nr1h3, cells = sample(colnames(human_nr1h3), 500))
if (ncol(ts_nr1h3)    > 500) ts_nr1h3    <- subset(ts_nr1h3,    cells = sample(colnames(ts_nr1h3),    500))

common_genes <- intersect(rownames(human_nr1h3), rownames(ts_nr1h3))
if (length(common_genes) < 50) warning("Few shared gene symbols; consider providing an ortholog map.")

human_mean <- avg_log_expr(human_nr1h3, genes = common_genes)
ts_mean    <- avg_log_expr(ts_nr1h3,    genes = common_genes)

## Pearson correlation
ct <- suppressWarnings(cor.test(human_mean, ts_mean, method = "pearson"))
r_val <- as.numeric(ct$estimate); p_val <- ct$p.value

## Plot
df <- data.frame(Human = human_mean, TreeShrew = ts_mean)
p <- ggplot(df, aes(Human, TreeShrew)) +
  geom_point(alpha = 0.6, color = "#1B9E77") +
  geom_smooth(method = "lm", se = FALSE, color = "#FF7F0E") +
  labs(
    title = "NR1H3+ Macrophages",
    subtitle = paste0("r = ", sprintf("%.3f", r_val), ",  p = ", format.pval(p_val, digits = 2)),
    x = "Human (log expression)",
    y = expression(italic("T. belangeri")*" (log expression)")
  ) +
  theme_classic(base_size = 11)

ggsave(out_pdf, p, width = 4.8, height = 4.2)
