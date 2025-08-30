suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

dir.create("results/cross_species", recursive = TRUE, showWarnings = FALSE)

species_files <- list(
  human     = "data/cross_species/human.rds",
  monkey    = "data/cross_species/monkey.rds",
  treeshrew = "data/cross_species/treeshrew.rds",
  mouse     = "data/cross_species/mouse.rds",
  pig       = "data/cross_species/pig.rds",
  rat       = "data/cross_species/rat.rds",
  zebrafish = "data/cross_species/zebrafish.rds"
)

objs <- lapply(names(species_files), function(sp) {
  obj <- readRDS(species_files[[sp]])
  DefaultAssay(obj) <- "RNA"
  if (!"orig.ident" %in% colnames(obj@meta.data)) obj$orig.ident <- sp
  obj
})
names(objs) <- names(species_files)

atlas <- merge(
  x = objs[[1]],
  y = objs[-1],
  add.cell.ids = names(objs),
  project = "CrossSpeciesSpleen"
)

## QC (same logic you used previously)
atlas$log10GenesPerUMI <- log10(atlas$nFeature_RNA) / log10(atlas$nCount_RNA)
before_n <- ncol(atlas)
atlas <- subset(
  atlas,
  subset = nCount_RNA > 500 & nCount_RNA < 20000 &
           nFeature_RNA < 3500 &
           log10GenesPerUMI > 0.80
)
message("Cells kept after QC: ", ncol(atlas), " / ", before_n)

## SCT + CCA integration
atlas_list <- SplitObject(atlas, split.by = "orig.ident")
atlas_list <- lapply(atlas_list, function(x) SCTransform(x, verbose = FALSE))

features <- SelectIntegrationFeatures(atlas_list, nfeatures = 3000)
atlas_list <- PrepSCTIntegration(atlas_list, anchor.features = features)

anchors <- FindIntegrationAnchors(
  object.list = atlas_list,
  normalization.method = "SCT",
  anchor.features = features
)
atlas_int <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

DefaultAssay(atlas_int) <- "integrated"
atlas_int <- ScaleData(atlas_int, verbose = FALSE)
atlas_int <- RunPCA(atlas_int, npcs = 50, verbose = FALSE)
dims_use <- 1:10
atlas_int <- FindNeighbors(atlas_int, reduction = "pca", dims = dims_use)
atlas_int <- FindClusters(atlas_int, resolution = 0.2)
atlas_int <- RunUMAP(atlas_int, reduction = "pca", dims = dims_use)

## Annotation
DefaultAssay(atlas_int) <- "RNA"
atlas_int <- NormalizeData(atlas_int, verbose = FALSE)

panel <- list(
  B_cells        = c("MS4A1","CD79A","CD79B"),
  Plasma_cells   = c("JCHAIN","MZB1","XBP1"),
  T_cells        = c("CD3D","CD3E","TRAC"),
  Macrophages    = c("CSF1R","C1QA","LST1","MS4A7")
)
panel <- lapply(panel, function(gs) intersect(gs, rownames(atlas_int)))

scores <- AddModuleScore(atlas_int, features = unname(panel), name = "coarse_", assay = "RNA")
sc_cols <- tail(colnames(scores@meta.data), length(panel))
colnames(scores@meta.data)[match(sc_cols, colnames(scores@meta.data))] <- paste0("score_", names(panel))
atlas_int <- scores

score_mat <- atlas_int@meta.data[, grep("^score_", colnames(atlas_int@meta.data)), drop = FALSE]
best_lab  <- names(panel)[max.col(as.matrix(score_mat), ties.method = "first")]
atlas_int$celltype_coarse <- factor(best_lab, levels = c("B_cells","Plasma_cells","T_cells","Macrophages"))

## Color palette
pal <- c(B_cells = "#DC050C", Plasma_cells = "#FB8072",
         T_cells = "#1965B0", Macrophages = "#7BAFDE")

## Figure A: UMAP by cell type
Idents(atlas_int) <- "celltype_coarse"

lab_counts <- table(atlas_int$celltype_coarse)
levels(atlas_int$celltype_coarse) <- paste0(levels(atlas_int$celltype_coarse),
                                            " (", as.integer(lab_counts[levels(atlas_int$celltype_coarse)]), ")")

pA <- DimPlot(
  atlas_int, reduction = "umap", group.by = "celltype_coarse",
  label = TRUE, repel = TRUE, cols = pal
) + theme(aspect.ratio = 1) + labs(title = "Integrated spleen across seven species")

ggsave("results/cross_species/FigA_UMAP_by_celltype.pdf", pA, width = 6.0, height = 5.0)

## Figure B: Feature plots of canonical markers
feat <- c("CD79A","JCHAIN","CD3E","CSF1R")
feat <- feat[feat %in% rownames(atlas_int)]
pB <- FeaturePlot(
  atlas_int, features = feat, reduction = "umap",
  min.cutoff = "q10", max.cutoff = "q90",
  order = TRUE, keep.scale = "feature", ncol = 2, pt.size = 0.1
) & theme(aspect.ratio = 1)

ggsave("results/cross_species/FigB_FeaturePlots_markers.pdf", pB, width = 6.5, height = 5.6)

saveRDS(atlas_int, "results/cross_species/atlas_integrated_SCT_CCA.rds")