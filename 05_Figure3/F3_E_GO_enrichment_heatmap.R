suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(cowplot)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(Matrix)
})

dir.create("results/macrophage", recursive = TRUE, showWarnings = FALSE)

mac <- readRDS("results/macrophage/macrophage_SCT_harmony.rds")
stopifnot("seurat_clusters" %in% colnames(mac@meta.data))

cl_num <- as.integer(as.character(mac$seurat_clusters))
cl_lvls <- sort(unique(cl_num))
mac$cluster_short <- factor(paste0("C", cl_num), levels = paste0("C", cl_lvls))
Idents(mac) <- "cluster_short"

rna_counts <- GetAssayData(mac, assay = "RNA", slot = "counts")
lib_size <- Matrix::colSums(rna_counts)
cpm <- t(t(rna_counts) / lib_size * 1e5)
logcpm <- log2(cpm + 1)
mac[["RNA"]]@data <- logcpm

degs <- FindAllMarkers(
  mac, only.pos = TRUE, test.use = "roc",
  logfc.threshold = 0.5, min.pct = 0.30, return.thresh = 0.25
)

degs_sig <- degs %>%
  filter(pct.1 > 0.30, power > 0.25) %>%
  arrange(cluster, desc(power))

write.csv(degs_sig, "results/macrophage/E_degs_sig.csv", row.names = FALSE)

## Top-50 per subcluster (unique genes; keep first hit if duplicated)
top50 <- degs_sig %>%
  group_by(cluster) %>% slice_max(order_by = power, n = 50, with_ties = FALSE) %>%
  ungroup() %>% distinct(gene, .keep_all = TRUE)

write.csv(top50, "results/macrophage/E_top50_markers.csv", row.names = FALSE)

## Heatmap ##
avg_mat <- AverageExpression(
  mac, features = top50$gene, assays = "RNA",
  group.by = "cluster_short", slot = "data"
)$RNA

avg_mat <- avg_mat[top50$gene, levels(mac$cluster_short), drop = FALSE]

z <- t(scale(t(as.matrix(avg_mat))))
z[is.na(z)] <- 0

ann_row <- data.frame(cluster = top50$cluster)
rownames(ann_row) <- top50$gene
cl_cols <- setNames(brewer.pal(max(3, min(9, length(levels(mac$cluster_short)))), "Set1"),
                    levels(mac$cluster_short))

ph <- pheatmap(
  z,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(99),
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = TRUE,
  annotation_row = ann_row,
  annotation_colors = list(cluster = cl_cols),
  border_color = NA,
  main = "Row z-score"
)

ggsave("results/macrophage/E_left_heatmap.pdf", ph, width = 5.2, height = 5)

## GO-BP enrichment (top-50 per subcluster, FDR<0.05; keep top-3 terms/cluster)
gene_map <- bitr(unique(top50$gene), fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db) %>% distinct(SYMBOL, .keep_all = TRUE)

gene_lists <- split(top50$gene, top50$cluster)

enrich_per_cluster <- function(sym_vec) {
  eg <- gene_map$ENTREZID[match(sym_vec, gene_map$SYMBOL)]
  eg <- unique(eg[!is.na(eg)])
  if (length(eg) < 5) return(NULL)
  res <- enrichGO(
    gene = eg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
    ont = "BP", pAdjustMethod = "BH",
    pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE
  )
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
  as.data.frame(res)
}

ego_list <- lapply(names(gene_lists), function(cl){
  df <- enrich_per_cluster(gene_lists[[cl]])
  if (is.null(df)) return(NULL)
  df$cluster <- cl
  df$n_deg <- length(gene_lists[[cl]])
  df
})
ego <- bind_rows(ego_list)

ego_top <- ego %>%
  filter(p.adjust < 0.05) %>%
  group_by(cluster) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  ungroup()

ego_top$set_size <- as.numeric(sub("/.*", "", ego_top$BgRatio))
ego_top$GeneRatio2 <- ego_top$Count / ego_top$set_size

ego_top <- ego_top %>%
  mutate(term_lab = Description) %>%
  arrange(cluster, p.adjust) %>%
  mutate(Order = row_number())

xmax <- max(-log10(ego_top$p.adjust), na.rm = TRUE)
xbreaks <- pretty(c(0, xmax), n = 5)

## Bubble plot
p_dot <- ggplot(ego_top,
  aes(x = -log10(p.adjust), y = reorder(term_lab, Order))) +
  geom_point(aes(size = Count, fill = GeneRatio2), shape = 21, color = "grey30") +
  scale_fill_gradientn(colours = c("grey80", "gold", "red")) +
  scale_size(range = c(2, 8)) +
  scale_x_continuous(breaks = xbreaks, limits = c(0, max(xbreaks))) +
  labs(x = expression(-log[10](P-value)), y = NULL,
       size = "Gene count", fill = "GeneRatio") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank())

ggsave("results/macrophage/E_right_GO_bubble.pdf", p_dot, width = 5.3, height = 5)

## Combined panel: heatmap+ GO bubble
gh <- ph$gtable
gd <- ggplotGrob(p_dot)

gd$heights[1] <- gh$heights[1]
combined <- cowplot::plot_grid(gh, gd, rel_widths = c(1, 1.35), nrow = 1, align = "v")

ggsave("results/macrophage/E_combined_heatmap_GO.pdf", combined, width = 10.5, height = 5)