suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ComplexHeatmap)
  library(grid)
})

tonsil <- readRDS("results/tonsil/tonsil.rds")

group_col <- dplyr::case_when(
  "final_fine" %in% colnames(tonsil@meta.data) ~ "final_fine",
  "ct_fine"    %in% colnames(tonsil@meta.data) ~ "ct_fine",
  TRUE ~ "seurat_clusters"
)
Idents(tonsil) <- tonsil@meta.data[[group_col]]

## Find markers with default settings (positive markers only)
all_markers <- FindAllMarkers(
  tonsil,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0.25
) %>% arrange(cluster, desc(avg_log2FC))

## Top 50 markers per cell type
top50 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50) %>%
  ungroup()

## Map to Entrez IDs (human tonsil assumed)
gene_map <- bitr(unique(top50$gene), fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db) %>% distinct(SYMBOL, .keep_all = TRUE)

gene_lists <- split(top50$gene, top50$cluster)

## Enrich KEGG for each cell type
ekegg_list <- lapply(names(gene_lists), function(cl){
  sym <- gene_lists[[cl]]
  entrez <- gene_map$ENTREZID[match(sym, gene_map$SYMBOL)]
  entrez <- unique(entrez[!is.na(entrez)])
  if (length(entrez) < 5) return(NULL)
  res <- enrichKEGG(gene = entrez, organism = "hsa",
                    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
  df <- as.data.frame(res)
  df$cluster <- cl
  df
})
ekegg_df <- bind_rows(ekegg_list)

## Keep significant (FDR < 0.05)
ekegg_sig <- ekegg_df %>%
  filter(p.adjust < 0.05) %>%
  group_by(cluster) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  ungroup()

## Build matrix
mat_df <- ekegg_sig %>%
  mutate(term = Description, value = -log10(p.adjust)) %>%
  select(term, cluster, value)

## Wide matrix
mat <- mat_df %>%
  tidyr::pivot_wider(names_from = cluster, values_from = value, values_fill = 0) %>%
  as.data.frame()
rownames(mat) <- mat$term; mat$term <- NULL
mat <- as.matrix(mat)

## Heatmap
col_fun <- colorRamp2(c(0, max(mat, na.rm = TRUE)),
                      c("#EFF3FF", "#08519C"))

ht <- Heatmap(
  mat,
  name = "-log10(FDR)",
  col = col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE, show_column_names = TRUE,
  heatmap_legend_param = list(title = "-log10(FDR)")
)

dir.create("results/tonsil", showWarnings = FALSE, recursive = TRUE)
pdf("results/tonsil/KEGG_Heatmap_Tonsil_D.pdf", width = 6.5, height = 6.5)
draw(ht)
dev.off()
