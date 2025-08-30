suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(muscat)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
})

dir.create("results/cross_species", recursive = TRUE, showWarnings = FALSE)

## Load integrated object and basic checks
atlas <- readRDS("results/cross_species/atlas_integrated_SCT_CCA.rds")
stopifnot(all(c("celltype_coarse","orig.ident") %in% colnames(atlas@meta.data)))

keep_types <- c("B_cells","Macrophages","Plasma_cells","T_cells")
atlas <- subset(atlas, subset = celltype_coarse %in% keep_types)

if ("unknown" %in% atlas$celltype_coarse) {
  atlas <- subset(atlas, subset = celltype_coarse != "unknown")
}

sample_max <- 1000
md <- atlas@meta.data |>
  tibble::rownames_to_column("cell_id") |>
  group_by(orig.ident, celltype_coarse) |>
  slice_sample(n = pmin(n(), sample_max)) |>
  ungroup()
atlas <- subset(atlas, cells = md$cell_id)

## Convert to SCE and set required columns for muscat
sce <- as.SingleCellExperiment(atlas)
colData(sce)$cluster_id <- atlas$celltype_coarse
colData(sce)$group_id   <- atlas$orig.ident            
colData(sce)$sample_id  <- atlas$orig.ident

sce <- muscat::prepSCE(sce, kid = "cluster_id", gid = "group_id", sid = "sample_id", drop = TRUE)

## Aggregate counts to pseudobulk 
pb <- muscat::aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id","sample_id"))
pb_counts <- assay(pb, "counts")

libsz   <- colSums(pb_counts)
logcpm  <- log1p(1e6 * t(t(pb_counts) / libsz))

var_genes <- head(order(matrixStats::rowVars(logcpm), decreasing = TRUE), 1000)
logcpm_top <- logcpm[var_genes, , drop = FALSE]

d   <- dist(t(logcpm_top), method = "euclidean")
mds <- cmdscale(d, k = 2)
df  <- as.data.frame(mds)
colnames(df) <- c("MDS1","MDS2")

cd <- as.data.frame(colData(pb))
df$Species  <- as.character(cd$sample_id)
df$Celltype <- factor(as.character(cd$cluster_id), levels = keep_types)

species_labs <- c(human = "human", monkey = "monkey", mouse = "mouse",
                  pig = "pig", rat = "rat", treeshrew = "T. belangeri",
                  zebrafish = "zebrafish")
df$Species_lab <- species_labs[df$Species]

pal_celltype <- c(
  "B_cells"       = "#DC050C",
  "Macrophages"   = "#FB8072",
  "Plasma_cells"  = "#1965B0",
  "T_cells"       = "#7BAFDE"
)

shape_map <- c(
  "human"        = 17,  # triangle
  "monkey"       = 4,   # cross
  "mouse"        = 5,   # diamond open
  "pig"          = 6,   # triangle down open
  "rat"          = 7,   # square cross
  "T. belangeri" = 8,   # star
  "zebrafish"    = 9    # diamond plus
)

p_mds <- ggplot(df, aes(x = MDS1, y = MDS2,
                        color = Celltype, shape = Species_lab)) +
  geom_hline(yintercept = 0, linetype = 3, color = "grey80") +
  geom_vline(xintercept = 0, linetype = 3, color = "grey80") +
  geom_point(size = 4, alpha = 0.8, stroke = 1) +
  scale_color_manual(values = pal_celltype, name = "Celltype") +
  scale_shape_manual(values = shape_map, name = "Species") +
  labs(x = "MDS1", y = "MDS2") +
  coord_equal() +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")

ggsave("results/cross_species/FigC_pseudobulk_MDS.pdf", p_mds, width = 4.2, height = 3.6)
