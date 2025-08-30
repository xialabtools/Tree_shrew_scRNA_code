suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

dir.create("results/cross_species", recursive = TRUE, showWarnings = FALSE)

atlas_file <- "results/cross_species/atlas_integrated_SCT_CCA.rds"
atlas <- readRDS(atlas_file)

label_col <- if ("celltype_coarse" %in% colnames(atlas@meta.data)) "celltype_coarse" else "celltype"
stopifnot(label_col %in% colnames(atlas@meta.data), "orig.ident" %in% colnames(atlas@meta.data))
mac <- subset(atlas, subset = get(label_col) == "Macrophages")

mac_list <- SplitObject(mac, split.by = "orig.ident")
mac_list <- lapply(mac_list, function(x) SCTransform(x, verbose = FALSE))
features <- SelectIntegrationFeatures(mac_list, nfeatures = 3000)
mac_list <- PrepSCTIntegration(mac_list, anchor.features = features)
anchors  <- FindIntegrationAnchors(object.list = mac_list,
                                   normalization.method = "SCT",
                                   anchor.features = features)
mac_int  <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

DefaultAssay(mac_int) <- "integrated"
mac_int <- ScaleData(mac_int, verbose = FALSE)
mac_int <- RunPCA(mac_int, npcs = 50, verbose = FALSE)
dims_use <- 1:30
mac_int <- FindNeighbors(mac_int, dims = dims_use)
mac_int <- FindClusters(mac_int, resolution = 0.2)
mac_int <- RunUMAP(mac_int, dims = dims_use)

species_order <- c("human","monkey","treeshrew","mouse","pig","rat","zebrafish")
mac_int$orig.ident <- factor(as.character(mac_int$orig.ident),
                             levels = intersect(species_order, unique(as.character(mac_int$orig.ident))))

pD <- DimPlot(mac_int, reduction = "umap", group.by = "orig.ident") +
      theme(aspect.ratio = 1) + labs(title = "Macrophages")

ggsave("results/cross_species/FigD_macrophage_umap_by_species.pdf", pD, width = 5.6, height = 4.4)

saveRDS(mac_int, "results/cross_species/macrophage_reclustered_SCT_CCA.rds")
