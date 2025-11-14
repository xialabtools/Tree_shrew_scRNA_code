suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(cowplot); library(Matrix); library(dplyr)
  library(purrr);  library(tibble);  library(patchwork); library(magrittr)
  library(SingleR); library(SummarizedExperiment)
  library(future)
})

theme_set(theme_cowplot())
set.seed(123)
options(future.globals.maxSize = 200 * 1024^3)

OUTDIR <- "./data/infection.spleen"
F4DIR  <- file.path(OUTDIR, "Fig4")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(F4DIR,  showWarnings = FALSE, recursive = TRUE)

if (!exists("mycolor_6")) {
  mycolor_6 <- c("#1f78b4","#33a02c","#e31a1c","#ff7f00","#6a3d9a","#b15928",
                 "#a6cee3","#b2df8a","#fb9a99","#fdbf6f","#cab2d6","#ffff99")
}

read_or_get <- function(obj_name, rds_path) {
  if (exists(obj_name, envir = .GlobalEnv)) get(obj_name, envir = .GlobalEnv) else readRDS(rds_path)
}
pth <- list(
  ctl = c("data/spleen/control/Spleen1.qc.rds","data/spleen/control/Spleen2.qc.rds","data/spleen/control/Spleen3.qc.rds"),
  d3  = c("data/spleen/3dpi/Day3_0.qc.rds","data/spleen/3dpi/Day3_1.qc.rds","data/spleen/3dpi/Day3_2.qc.rds","data/spleen/3dpi/Day3_3.qc.rds"),
  d7  = c("data/spleen/7dpi/Day7_0.qc.rds","data/spleen/7dpi/Day7_1.qc.rds","data/spleen/7dpi/Day7_2.qc.rds","data/spleen/7dpi/Day7_3.qc.rds")
)

add_inf <- function(obj, lab){ obj$infection <- lab; obj }

## Control sample
Spleen1.qc <- read_or_get("Spleen1.qc", pth$ctl[1])
Spleen2.qc <- read_or_get("Spleen2.qc", pth$ctl[2])
Spleen3.qc <- read_or_get("Spleen3.qc", pth$ctl[3])
Spleen <- merge(Spleen1.qc, c(Spleen2.qc,Spleen3.qc), add.cell.ids=c("Spleen1","Spleen2","Spleen3")) %>% add_inf("0dpi")
Spleen@meta.data <- Spleen@meta.data[, c("orig.ident","nCount_RNA","nFeature_RNA","infection"), drop=FALSE]

## 3dpi sample
Day3_0 <- read_or_get("Spleen_EBV_Day3_0.qc", pth$d3[1]) %>% add_inf("3dpi")
Day3_1 <- read_or_get("Spleen_EBV_Day3_1.qc", pth$d3[2]) %>% add_inf("3dpi")
Day3_2 <- read_or_get("Spleen_EBV_Day3_2.qc", pth$d3[3]) %>% add_inf("3dpi")
Day3_3 <- read_or_get("Spleen_EBV_Day3_3.qc", pth$d3[4]) %>% add_inf("3dpi")
Spleen_day3 <- merge(Day3_0, c(Day3_1,Day3_2,Day3_3), add.cell.ids=c("Day3_0","Day3_1","Day3_2","Day3_3"))

## 7dpi sample
Day7_0 <- read_or_get("Spleen_EBV_Day7_0.qc", pth$d7[1]) %>% add_inf("7dpi")
Day7_1 <- read_or_get("Spleen_EBV_Day7_1.qc", pth$d7[2]) %>% add_inf("7dpi")
Day7_2 <- read_or_get("Spleen_EBV_Day7_2.qc", pth$d7[3]) %>% add_inf("7dpi")
Day7_3 <- read_or_get("Spleen_EBV_Day7_3.qc", pth$d7[4]) %>% add_inf("7dpi")
Spleen_day7 <- merge(Day7_0, c(Day7_1,Day7_2,Day7_3), add.cell.ids=c("Day7_0","Day7_1","Day7_2","Day7_3"))

##SingleR 
singleR_annotate <- function(obj, ref_rdata = "data/ref/ref_BE_259s.RData"){
  load(ref_rdata)  
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj, features = rownames(obj))
  obj <- RunPCA(obj, features = VariableFeatures(obj))
  obj <- RunUMAP(obj, dims = 1:20) %>% FindNeighbors(dims = 1:20) %>% FindClusters(resolution = 0.5)
  pred <- SingleR(test=as.matrix(GetAssayData(obj, "RNA", "data")), ref=ref.BE, labels=ref.BE$label.main)
  obj$singleR <- pred$labels[match(rownames(obj@meta.data), rownames(pred))]
  obj@meta.data <- obj@meta.data[, c("orig.ident","nCount_RNA","nFeature_RNA","infection","singleR")]
  obj
}

Spleen       <- singleR_annotate(Spleen)
Spleen_day3  <- singleR_annotate(Spleen_day3)
Spleen_day7  <- singleR_annotate(Spleen_day7)

scRNA <- merge(Spleen, c(Spleen_day3,Spleen_day7))

## SCT/CCA  intergrated
split_seurat <- SplitObject(scRNA, split.by = "orig.ident")
split_seurat <- lapply(split_seurat, SCTransform)
integ_features <- SelectIntegrationFeatures(split_seurat, nfeatures = 2000)
split_seurat  <- PrepSCTIntegration(split_seurat, anchor.features = integ_features)
plan("multisession", workers = 6)
integ_anchors <- FindIntegrationAnchors(split_seurat, normalization.method = "SCT", anchor.features = integ_features)
seurat_integrated <- IntegrateData(integ_anchors, normalization.method = "SCT")
plan("sequential"); gc()

scRNA <- RunPCA(seurat_integrated, npcs = 50, verbose = FALSE)
pc.use <- 1:10
scRNA <- RunUMAP(scRNA, dims = pc.use) %>% FindNeighbors(dims = pc.use) %>% FindClusters(resolution = 0.8)
saveRDS(scRNA, file = file.path(OUTDIR, "infection.spleen.SCT.CCA.rds"))

cluster_to_type <- c("0"="B","1"="B","4"="B","2"="CD4T","8"="CD4T","5"="CD8T","9"="CD8T",
                     "6"="IFNG_NKT","3"="NKT","10"="Plasma_B","7"="Neutrophil",
                     "13"="Cycling_B","14"="Cycling_B","12"="NR1H3_Macro","11"="IRF8_Macro")
scRNA$celltype <- unname(cluster_to_type[as.character(scRNA$seurat_clusters)])
Idents(scRNA) <- "celltype"
scRNA$infection <- factor(scRNA$infection, levels = c("0dpi","3dpi","7dpi"))
saveRDS(scRNA, file = file.path(OUTDIR, "scRNA.finall.rds"))

## Fig.4A UMAP
pA1 <- DimPlot(scRNA, group.by="celltype", label=TRUE, cols = mycolor_6)
pA2 <- DimPlot(scRNA, group.by="celltype", split.by="infection", label=TRUE, cols = mycolor_6, ncol=3)
ggsave(file.path(F4DIR, "F4A_umap_celltype.pdf"), pA1, width=6.5, height=5)
ggsave(file.path(F4DIR, "F4A_umap_celltype_byTime.pdf"), pA2, width=10, height=4.6)