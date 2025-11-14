suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(CellChat); library(future)
})

OUTDIR <- "./data/infection.spleen"
F4DIR  <- file.path(OUTDIR, "Fig4")
dir.create(F4DIR, showWarnings = FALSE, recursive = TRUE)

scRNA <- readRDS(file.path(OUTDIR, "scRNA.finall.rds"))

Idents(scRNA) <- "infection"
data_list <- list(Control=subset(scRNA, idents="0dpi"),
                  EBV_Day3=subset(scRNA, idents="3dpi"),
                  EBV_Day7=subset(scRNA, idents="7dpi"))

TS_cellchat <- function(
    Seuratobject,
    select_meta_col,
    organism,
    min.cells
){
  
  datainput <- Seuratobject[["RNA"]]@data
  data_meta <- Seuratobject@meta.data
  data_meta <- data_meta[,select_meta_col]
  colnames(data_meta) <- c("group","labels")
  
  data_cellchat <- createCellChat(object = datainput, 
                                  meta = data_meta, 
                                  group.by = "labels")
  
  data_cellchat <- addMeta(data_cellchat, meta = data_meta)
  data_cellchat <- setIdent(data_cellchat, ident.use = "labels")
  groupSize <- as.numeric(table(data_cellchat@idents)) 
  
  if(organism=='human'){
    CellChatDB = CellChatDB.human
    data_cellchat@DB <- CellChatDB
    
    data_cellchat <- subsetData(data_cellchat, features = NULL)
    data_cellchat <- identifyOverExpressedGenes(data_cellchat)
    data_cellchat <- identifyOverExpressedInteractions(data_cellchat)
    
    data_cellchat <- projectData(data_cellchat, PPI.human)
    data_cellchat <- computeCommunProb(data_cellchat, raw.use=T)
    data_cellchat <- filterCommunication(data_cellchat, min.cells = min.cells)
    data_cellchat <- computeCommunProbPathway(data_cellchat)
    data_cellchat <- aggregateNet(data_cellchat)
    data_cellchat <- netAnalysis_computeCentrality(data_cellchat, slot.name = "netP")
    
  }else{
    
    CellChatDB = CellChatDB.mouse
    data_cellchat@DB <- CellChatDB
    
    data_cellchat <- subsetData(data_cellchat, features = NULL)
    data_cellchat <- identifyOverExpressedGenes(data_cellchat)
    data_cellchat <- identifyOverExpressedInteractions(data_cellchat)
    
    data_cellchat <- projectData(data_cellchat, PPI.mouse)
    data_cellchat <- computeCommunProb(data_cellchat, raw.use=T)
    data_cellchat <- filterCommunication(data_cellchat, min.cells = min.cells)
    data_cellchat <- computeCommunProbPathway(data_cellchat)
    data_cellchat <- aggregateNet(data_cellchat)
  }
  
  return(data_cellchat)
  
}

plan("multisession", workers = 24); options(future.globals.maxSize = 200 * 1024^3)
cellchat_list <- lapply(data_list, TS_cellchat, select_meta_col=c("infection","celltype"), min.cells=10)
plan("sequential"); names(cellchat_list) <- names(data_list)

saveRDS(cellchat_list, file = file.path(OUTDIR, "cellchat_list.EBV.infection.spleen.rds"))

## Fig.4C
cellchat <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))

pC1 <- compareInteractions(cellchat, show.legend=FALSE, group=c(1,2,3))
pC2 <- compareInteractions(cellchat, show.legend=FALSE, group=c(1,2,3), measure="weight")
ggsave(file.path(F4DIR, "F4C_interactions_count.pdf"),  pC1, width=5.2, height=3.6)
ggsave(file.path(F4DIR, "F4C_interactions_weight.pdf"), pC2, width=5.2, height=3.6)

## Fig.4D
num.link <- sapply(object_list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) 
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], 
                                               weight.MinMax = weight.MinMax,
                                               dot.size = c(3, 6),
                                               label.size = 6)
}

patchwork::wrap_plots(plots = gg)
