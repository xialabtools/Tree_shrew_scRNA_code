suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(purrr)
  library(Seurat)
  library(patchwork)
})

samples <- c(
  "Heart_1","Kidney_1","Liver_1","Lung_1",
  "Lymph_1","Lymph_2","Lymph_3",
  "Marrow_1","Marrow_2","Marrow_3",
  "Nasopharynx_1",
  "PBMC_1","PBMC_2","PBMC_3",
  "Spleen_1","Spleen_2","Spleen_3",
  "Stomach_1",
  "Thymus_1","Thymus_2","Thymus_3",
  "Tonsil_1","Tonsil_2","Tonsil_3"
)

cfg <- list(qc_dir = "QC", rds_dir = "data/qc", out = "QC/qc_overview.pdf")

read_meta_one <- function(s) {
  csv <- file.path(cfg$qc_dir, paste0(s, "_meta.csv"))
  if (file.exists(csv)) {
    x <- suppressMessages(read_csv(csv, show_col_types = FALSE))
    x$sample_id <- s
    return(x)
  }
  rds <- file.path(cfg$rds_dir, paste0(s, ".qc.rds"))
  if (file.exists(rds)) {
    obj <- readRDS(rds)
    md <- obj@meta.data
    md$sample_id <- s
    return(as_tibble(md[, c("sample_id","nFeature_RNA","nCount_RNA")]))
  }
  tibble(sample_id = character(0), nFeature_RNA = numeric(0), nCount_RNA = numeric(0))
}

meta <- map_dfr(samples, read_meta_one) |>
  mutate(sample_id = factor(sample_id, levels = rev(samples))) |>
  select(sample_id, nFeature_RNA, nCount_RNA)

pal <- setNames(scales::hue_pal()(length(samples)), samples)

counts_df <- meta |>
  count(sample_id, name = "cells")

p_left <- ggplot(counts_df, aes(x = cells, y = sample_id, fill = sample_id)) +
  geom_col(width = 0.7, color = NA) +
  scale_fill_manual(values = pal, guide = "none") +
  labs(x = "Number", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.title.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.major.x = element_blank()
  )

p_mid <- ggplot(meta, aes(x = nFeature_RNA, y = sample_id, fill = sample_id)) +
  geom_boxplot(width = 0.6, outlier.shape = 21, outlier.size = 0.8, outlier.stroke = 0.2) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_continuous(trans = "log2", breaks = c(512,1024,2048,4096,8192)) +
  labs(x = "Genes", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.major.x = element_blank()
  )

p_right <- ggplot(meta, aes(x = nCount_RNA, y = sample_id, fill = sample_id)) +
  geom_boxplot(width = 0.6, outlier.shape = 21, outlier.size = 0.8, outlier.stroke = 0.2) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_continuous(trans = "log2", breaks = c(4096,8192,16384,32768)) +
  labs(x = "UMI", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.major.x = element_blank()
  )

g <- p_left + p_mid + p_right + plot_layout(widths = c(1.25, 1, 1))
ggsave(cfg$out, g, width = 10, height = 7, limitsize = FALSE)