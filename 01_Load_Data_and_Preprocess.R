suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(Matrix)
  library(dplyr)
  library(purrr)
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

cfg <- list(
  data_dir = "./data",
  suffix = "_NCBI_filtered_feature_bc_matrix",
  out_dir = "QC",
  obj_dir = "data/qc",
  dims = 1:15,
  doublet_rate = 0.076,
  lower_pct = 0.01,
  upper_pct = 0.99,
  min_cells_per_gene = 10,
  complex_min = 0.8,
  mt_max = 30,
  hb_max = 1
)

dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$obj_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("CaseMatch")) {
  CaseMatch <- function(search, match) {
    match[match(tolower(search), tolower(match))]
  }
}

build_sample_info <- function(x, batch_map = NULL, species = "Tupaia belangeri", platform = "10x") {
  data.frame(
    sample_id = x,
    tissue = sub("_\\d+$", "", x),
    replicate = sub("^.*_", "", x),
    batch = if (is.null(batch_map)) NA_character_ else unname(batch_map[x]),
    species = species,
    platform = platform,
    stringsAsFactors = FALSE
  )
}

load_sample <- function(id) {
  m <- Read10X(file.path(cfg$data_dir, paste0(id, cfg$suffix)))
  CreateSeuratObject(counts = m, project = id, min.cells = 3, min.features = 200)
}

add_sample_meta <- function(obj, info) {
  obj$sample_id <- info$sample_id
  obj$source <- info$tissue
  obj$batch <- info$batch
  obj$replicate <- info$replicate
  obj
}

call_doublets <- function(obj) {
  obj <- SCTransform(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- RunUMAP(obj, dims = cfg$dims, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = cfg$dims, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  sweeps <- paramSweep_v3(obj, PCs = cfg$dims, sct = TRUE)
  stats <- summarizeSweep(sweeps, GT = FALSE)
  bcmvn <- find.pK(stats)
  pk <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  hom <- modelHomotypic(obj$seurat_clusters)
  nExp <- round(cfg$doublet_rate * ncol(obj))
  nExpAdj <- round(nExp * (1 - hom))
  obj <- doubletFinder_v3(
    obj, PCs = cfg$dims, pN = 0.25, pK = pk, nExp = nExpAdj,
    reuse.pANN = FALSE, sct = TRUE
  )
  cls_col <- grep("^DF\\.classifications", colnames(obj@meta.data), value = TRUE)
  obj$doublet_class <- obj@meta.data[[cls_col[length(cls_col)]]]
  obj <- subset(obj, subset = doublet_class == "Singlet")
  DefaultAssay(obj) <- "RNA"
  obj@assays$SCT <- NULL
  attr(obj, "pK") <- pk
  obj
}

compute_qc_metrics <- function(obj) {
  obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)
  if (any(grepl("^MT-", rownames(obj)))) {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  } else {
    obj[["percent.mt"]] <- 0
  }
  obj[["percent.rb"]] <- PercentageFeatureSet(obj, pattern = "^RP[LS]")
  hb <- CaseMatch(
    c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ"),
    rownames(obj)
  )
  hb <- hb[!is.na(hb)]
  if (length(hb)) {
    obj[["percent.hb"]] <- PercentageFeatureSet(obj, features = hb)
  } else {
    obj[["percent.hb"]] <- 0
  }
  obj
}

apply_cell_qc <- function(obj) {
  nc <- quantile(obj$nCount_RNA, probs = c(cfg$lower_pct, cfg$upper_pct), na.rm = TRUE, names = FALSE)
  nf <- quantile(obj$nFeature_RNA, probs = c(cfg$lower_pct, cfg$upper_pct), na.rm = TRUE, names = FALSE)
  obj <- subset(
    obj,
    subset =
      log10GenesPerUMI > cfg$complex_min &
      percent.mt < cfg$mt_max &
      percent.hb < cfg$hb_max &
      nCount_RNA > nc[1] & nCount_RNA < nc[2] &
      nFeature_RNA > nf[1] & nFeature_RNA < nf[2]
  )
  attr(obj, "ncount_cut") <- nc
  attr(obj, "nfeature_cut") <- nf
  obj
}

score_cell_cycle <- function(obj) {
  obj <- NormalizeData(obj, verbose = FALSE)
   {s <- CaseMatch(cc.genes$s.genes, rownames(obj)); s <- s[!is.na(s)]
    g <- CaseMatch(cc.genes$g2m.genes, rownames(obj)); g <- g[!is.na(g)]
  } else {
    s <- CaseMatch(cc.genes$s.genes, rownames(obj)); s <- s[!is.na(s)]
    g <- CaseMatch(cc.genes$g2m.genes, rownames(obj)); g <- g[!is.na(g)]
  }
  CellCycleScoring(obj, s.features = s, g2m.features = g)
}

filter_genes_by_cells <- function(obj, min_cells = cfg$min_cells_per_gene) {
  cts <- GetAssayData(obj, slot = "counts")
  keep <- Matrix::rowSums(cts > 0) >= min_cells
  CreateSeuratObject(cts[keep, , drop = FALSE], meta.data = obj@meta.data, project = unique(obj$orig.ident))
}

process_one <- function(id, info_tbl) {
  info <- info_tbl[info_tbl$sample_id == id, , drop = FALSE]
  o0 <- load_sample(id)
  o0 <- add_sample_meta(o0, info)
  n_raw <- ncol(o0)
  o1 <- call_doublets(o0)
  o1 <- compute_qc_metrics(o1)
  n_post_doublet <- ncol(o1)
  o2 <- apply_cell_qc(o1)
  o2 <- score_cell_cycle(o2)
  o3 <- filter_genes_by_cells(o2)
  n_post_qc <- ncol(o3)
  pk <- attr(o1, "pK")
  nc <- attr(o2, "ncount_cut")
  nf <- attr(o2, "nfeature_cut")
  saveRDS(o3, file.path(cfg$obj_dir, paste0(id, ".qc.rds")))
  }

## Batch processing
if (sys.nframe() == 0) {
  dir.create(cfg$out_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(cfg$obj_dir,  recursive = TRUE, showWarnings = FALSE)

  if (!exists("process_one")) {
    stop("process_one() must be defined above before batch processing.")
  }

  .msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

  infer_input_dir <- function(sample_id) {
    cand <- c(
      file.path(cfg$data_dir, sample_id, paste0(sample_id, cfg$suffix)),
      file.path(cfg$data_dir, sample_id, "outs", "filtered_feature_bc_matrix"),
      file.path(cfg$data_dir, sample_id)
    )
    cand <- cand[dir.exists(cand)]
    if (length(cand) == 0) NA_character_ else cand[[1]]
  }

  fml <- names(formals(process_one))

  .msg("Found ", length(samples), " sample(s): ", paste(samples, collapse = ", "))

  for (sid in samples) {
    input_dir <- infer_input_dir(sid)
    out_file  <- file.path(cfg$obj_dir, paste0(sid, ".qc.rds"))

    if (file.exists(out_file)) { .msg("[Skip]", sid, "already has", out_file); next }

    if (is.na(input_dir)) {
      .msg("[Warn]", sid, ": input dir not found under ", cfg$data_dir, "; skip")
      next
    }

    .msg("[Start]", sid, "->", input_dir)
    t0 <- Sys.time()

    args <- list()
    if ("sample_id" %in% fml)      args$sample_id <- sid else if ("sid" %in% fml) args$sid <- sid
    if ("input_dir" %in% fml)      args$input_dir <- input_dir
    if ("cfg" %in% fml)            args$cfg       <- cfg

    res <- try(do.call(process_one, args), silent = TRUE)
    if (inherits(res, "try-error")) {
      .msg("[Error]", sid, ":", conditionMessage(attr(res, "condition"))); next
    }

    if (!file.exists(out_file) && !is.null(res)) {
      try(saveRDS(res, out_file), silent = TRUE)
    }

    .msg("[Done]", sid, " (",
         round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
         " sec ) -> ", out_file)
  }
}
