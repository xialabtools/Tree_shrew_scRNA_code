suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(ReactomePA)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)

counts_file <- "data/bulk/NR1H3_OE_WT_counts.tsv"  # rows: SYMBOL; cols: samples (raw counts)
meta_file   <- "data/bulk/NR1H3_OE_WT_meta.tsv"    # columns: sample, group (WT|NR1H3_OE)

cts  <- read_tsv(counts_file)
meta <- read_tsv(meta_file)

gene_col <- colnames(cts)[1]
cts_mat  <- as.matrix(cts[,-1]); rownames(cts_mat) <- cts[[gene_col]]
cts_mat  <- cts_mat[, meta$sample, drop = FALSE]
keep     <- rowSums(cts_mat >= 10) >= 2
cts_mat  <- cts_mat[keep, ]

## -------- DESeq2 --------
meta$group <- factor(meta$group, levels = c("WT","NR1H3_OE"))
dds <- DESeqDataSetFromMatrix(round(cts_mat), as.data.frame(meta), design = ~ group)
dds <- DESeq(dds)

coef_name <- grep("group_", resultsNames(dds), value = TRUE)[1]
res_shr   <- lfcShrink(dds, coef = coef_name, type = "apeglm")
res_raw   <- results(dds, contrast = c("group","NR1H3_OE","WT"))

de_tbl <- as.data.frame(res_shr) |>
  rownames_to_column("gene") |>
  left_join(as.data.frame(res_raw)[,c("stat")] |> rownames_to_column("gene"),
            by = "gene") |>
  arrange(padj, desc(abs(log2FoldChange)))

write.csv(de_tbl, "results/deseq2/deseq2_all_results.csv", row.names = FALSE)

## -------- Panel E: volcano --------
padj_thr <- 0.05; lfc_thr <- 1
volcano_df <- de_tbl |>
  mutate(cat = case_when(
    is.na(padj) ~ "Not Sig",
    padj < padj_thr & abs(log2FoldChange) > lfc_thr ~ "Pvalue and Log2FC",
    padj < padj_thr ~ "Pvalue",
    abs(log2FoldChange) > lfc_thr ~ "Log2FC",
    TRUE ~ "Not Sig"
  ),
  neg_log10_p = -log10(pvalue))

cols <- c("Not Sig" = "grey70", "Log2FC" = "#66bd63", "Pvalue" = "#2c7fb8",
          "Pvalue and Log2FC" = "#f03b20")

genes_to_label <- c("TNF","TNFSF14","C1QA","C1QC","NFKBIZ","NFKB2","MAP3K14","LTBR","NR1H3")
volcano_df$label <- ifelse(volcano_df$gene %in% genes_to_label, volcano_df$gene, NA)

p_volcano <- ggplot(volcano_df, aes(log2FoldChange, neg_log10_p, color = cat)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = cols, name = NULL) +
  geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = 2, color = "grey50") +
  geom_hline(yintercept = -log10(padj_thr), linetype = 2, color = "grey50") +
  geom_text_repel(aes(label = label), size = 3, min.segment.length = 0,
                  box.padding = 0.3, point.padding = 0.2, max.overlaps = 50) +
  labs(x = "Log2 fold change", y = expression(-log[10](pvalue))) +
  theme_classic(base_size = 11)

ggsave("results/deseq2/panel_E_volcano.pdf", p_volcano, width = 4.6, height = 3.8)

## -------- Panel F: Reactome GSEA  --------
sym <- de_tbl$gene
stat <- de_tbl$stat
conv <- bitr(sym, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) |>
  distinct(SYMBOL, .keep_all = TRUE)

geneList <- stat[match(conv$SYMBOL, sym)]
names(geneList) <- conv$ENTREZID
geneList <- geneList[!is.na(geneList)]
geneList <- sort(geneList, decreasing = TRUE)

## Run Reactome GSEA
gsea_res <- gsePathway(geneList, organism = "human",
                       pvalueCutoff = 0.2, pAdjustMethod = "BH",
                       minGSSize = 10, maxGSSize = 500, verbose = FALSE)
saveRDS(gsea_res, "results/deseq2/reactome_gsea_result.rds")


target <- "TNF receptor superfamily members mediating non canonical NF-kB pathway"
hit_idx <- which(grepl("TNF.*non.?canonical.*NF.?kB", gsea_res@result$Description, ignore.case = TRUE))
if (length(hit_idx) == 0) hit_idx <- 1  ## fallback to top term

## ----- Custom GSEA plot with gene labels -----
## Compute running ES and add ticks + symbols for core-enrichment genes
plot_gsea_labeled <- function(gsea_res, geneList, geneSetID, label_genes = NULL,
                              label_n = 12, title_prefix = NULL) {
  res <- gsea_res@result[geneSetID, ]
  ord <- sort(geneList, decreasing = TRUE)
  N   <- length(ord)

  core_eg  <- unlist(strsplit(res$core_enrichment, "/"))
  core_pos <- match(core_eg, names(ord))
  core_pos <- core_pos[!is.na(core_pos)]
  map_tab  <- bitr(core_eg, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  map_tab  <- map_tab[match(core_eg, map_tab$ENTREZID), , drop = FALSE]
  core_sym <- map_tab$SYMBOL

  p <- 1
  hits     <- rep(FALSE, N); hits[core_pos] <- TRUE
  sum_hit  <- sum(abs(ord[hits])^p)
  P_hit    <- ifelse(hits, abs(ord)^p / sum_hit, 0)
  P_miss   <- ifelse(!hits, 1/(N - length(core_pos)), 0)
  running  <- cumsum(P_hit - P_miss)
  df_es    <- data.frame(rank = 1:N, ES = running)

  ticks <- data.frame(rank = core_pos)

  label_df <- data.frame(rank = core_pos, SYMBOL = core_sym)
  if (!is.null(label_genes)) {
    label_df <- subset(label_df, SYMBOL %in% label_genes)
  } else {
    abs_w <- abs(ord[core_pos])
    ord_ix <- order(abs_w, decreasing = TRUE)
    label_df <- label_df[ord_ix[seq_len(min(label_n, nrow(label_df)))], , drop = FALSE]
  }
  label_df$y <- max(df_es$ES, na.rm = TRUE) * 0.15

  ttl <- if (is.null(title_prefix)) res$Description else paste0(title_prefix, "\n", res$Description)
  stat_lab <- paste0("NES: ", sprintf("%.2f", res$NES),
                     "   P: ", format(res$pvalue, digits = 2),
                     "   FDR: ", format(res$qvalues, digits = 2))

  ggplot(df_es, aes(rank, ES)) +
    geom_line(color = "#B2182B", linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_segment(data = ticks, aes(x = rank, xend = rank, y = 0, yend = 0.08 * max(ES)),
                 color = "grey40", linewidth = 0.2) +
    geom_text_repel(data = label_df, aes(x = rank, y = y, label = SYMBOL),
                    size = 3, min.segment.length = 0, box.padding = 0.25,
                    segment.color = "grey50", max.overlaps = Inf) +
    labs(title = ttl, subtitle = stat_lab,
         x = "Rank in ordered dataset", y = "Running Enrichment Score") +
    theme_classic(base_size = 11)
}

label_genes_highlight <- c("MAP3K14","LTBR","TNFSF14","CD40LG","TNFRSF12A","TNFRSF13B","LTA",
                           "TNFSF12","NFKB2","NFKBIZ")

p_gsea <- plot_gsea_labeled(gsea_res, geneList, hit_idx,
                            label_genes = label_genes_highlight,
                            label_n = 12,
                            title_prefix = "TNF Receptor Superfamily Members Mediating\nNon‑Canonical NF‑κB Pathway")

ggsave("results/deseq2/panel_F_gsea_reactome_labeled.pdf", p_gsea, width = 5.2, height = 3.8)
