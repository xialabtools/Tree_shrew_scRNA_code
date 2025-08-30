suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(WGCNA)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(patchwork)
  library(enrichR)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

dir.create("results/wgcna", recursive = TRUE, showWarnings = FALSE)

## Load macrophage object
obj_file <- "results/cross_species/macrophage_reclustered_SCT_CCA.rds"
if (file.exists(obj_file)) {
  sce <- readRDS(obj_file)
} else if (exists("macro.final")) {
  sce <- macro.final
} else {
  stop("No macrophage object found. Supply 'macro.final' or the RDS from Fig D.")
}
DefaultAssay(sce) <- "RNA"

if (!"celltype_coarse" %in% colnames(sce@meta.data) &&
    !"celltype" %in% colnames(sce@meta.data)) {
  sce$celltype <- "Macrophages"
}
label_col <- if ("celltype_coarse" %in% colnames(sce@meta.data)) "celltype_coarse" else "celltype"
sce <- subset(sce, subset = get(label_col) == "Macrophages")

sp_to_cap <- c("pig","rat","treeshrew","zebrafish")
keep_cells <- unlist(lapply(split(colnames(sce), sce$orig.ident), function(ids) {
  sp <- unique(sce$orig.ident[ids])[1]
  if (sp %in% sp_to_cap && length(ids) > 200) sample(ids, 200) else ids
}), use.names = FALSE)
sce <- subset(sce, cells = keep_cells)

DefaultAssay(sce) <- "RNA"
mt_genes   <- grep("^MT",  rownames(sce), ignore.case = TRUE, value = TRUE)
ribo_genes <- grep("^(RPS|RPL)", rownames(sce), ignore.case = TRUE, value = TRUE)
genes_keep <- setdiff(rownames(sce), unique(c(mt_genes, ribo_genes)))
sce <- subset(sce, features = genes_keep)

if (!"celltype" %in% colnames(sce@meta.data)) sce$celltype <- "Macrophages"

sce <- SetupForWGCNA(
  sce, wgcna_name = "crossSpec",
  gene_select = "fraction", fraction = 0.05,
  group.by = "orig.ident"
)

sce <- MetacellsByGroups(
  seurat_obj = sce, group.by = "orig.ident",
  k = 25, max_shared = 10, reduction = "pca",
  min_cells = 50, ident.group = "orig.ident"
)
sce <- NormalizeMetacells(sce)
sce <- SetDatExpr(sce, assay = "RNA", slot = "data")

sce <- TestSoftPowers(sce, networkType = "signed")
sce <- ConstructNetwork(
  sce, soft_power = 10,
  tom_name = "crossSpec_TOM", overwrite_tom = TRUE,
  setDatExpr = FALSE
)

sce <- ModuleEigengenes(sce, group.by.vars = "orig.ident")
sce <- ModuleConnectivity(sce)
sce <- ResetModuleNames(sce, new_name = "M")

## Panel F: module–species correlation heatmap with significance stars
sp_levels <- c("human","monkey","treeshrew","mouse","pig","rat","zebrafish")
sce$orig.ident <- factor(as.character(sce$orig.ident), levels = sp_levels)
for (sp in sp_levels) sce[[sp]] <- as.factor(sce$orig.ident == sp)

traits <- sp_levels               
sce <- ModuleTraitCorrelation(
  sce, traits = traits,
  group.by = "celltype",        
  features = "hMEs"
)

pdf("results/wgcna/panel_F_module_species_correlation.pdf", width = 12, height = 7)
PlotModuleTraitCorrelation(
  sce,
  label = "fdr",                   
  label_symbol = "stars",        
  text_size = 4,
  high_color = "#b2182b",
  mid_color  = "#f7fbff",
  low_color  = "#2166ac",
  combine = TRUE
)
dev.off()

## Panel G: enrichment analysis
dbs <- c('GO_Biological_Process_2023')

sce <- RunEnrichr(
  sce,
  dbs=dbs, 
  max_genes = 100
)

enrich_df <- GetEnrichrTable(sce)
write.csv(enrich_df,"enrich_GO.csv")

pdf("10_module_GO.species.pdf",width = 10,height =10)
EnrichrDotPlot(
  sce,
  mods = c("M1","M8","M5","M10","M17","M2","M4"),
  database = "GO_Biological_Process_2023", 
  n_terms=3 )

dev.off()