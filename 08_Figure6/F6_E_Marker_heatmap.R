suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(pheatmap)
  library(RColorBrewer)
})

dir.create("results/cross_species", recursive = TRUE, showWarnings = FALSE)

mac_int <- readRDS("results/cross_species/macrophage_reclustered_SCT_CCA.rds")

DefaultAssay(mac_int) <- "RNA"
mac_int <- NormalizeData(mac_int, verbose = FALSE)

genes_conserved <- c("CSF1R","CFP","CTSB","PSAP","LGMN","ASAH1","NR1H3")
genes_antigen   <- c("HLA.A","HLA.B","HLA.C","HLA.DMA","HLA.DRA","HLA.E","HLA.DPB1")
genes_cytokine  <- c("IL1B","IL4R","IL10RA","TGFB1","IFNAR1","IFNGR1","TNFAIP8","TNFRSF14","TNFRSF1B")
genes_chemokine <- c("CCL4","CCL5","CXCR4","CXCL10","CXCL16")
genes_complement<- c("C1QA","C1QB","C1QBP","C1QC")

heatmap_genes  <- c(genes_conserved, genes_antigen, genes_cytokine, genes_chemokine, genes_complement)

Idents(mac_int) <- "orig.ident"
avg <- AverageExpression(mac_int, assays = "RNA", features = heatmap_genes)$RNA
mat <- t(avg)

species_order <- c("human","monkey","treeshrew","mouse","pig","rat","zebrafish")
mat <- mat[intersect(species_order, rownames(mat)), , drop = FALSE]

group_labels <- c(rep("Conserved genes", length(genes_conserved)),
                  rep("Antigen Presentation", length(genes_antigen)),
                  rep("Cytokine", length(genes_cytokine)),
                  rep("Chemokine", length(genes_chemokine)),
                  rep("Complement", length(genes_complement)))
anno_col <- data.frame(Function = factor(group_labels,
                                         levels = c("Conserved genes","Antigen Presentation",
                                                    "Cytokine","Chemokine","Complement")))
rownames(anno_col) <- colnames(mat)

sp_counts <- table(mac_int$orig.ident)
rownames(mat) <- paste0(rownames(mat), " (", as.integer(sp_counts[rownames(mat)]), ")")

pal <- colorRampPalette(c("#2166AC","#F7FBFF","#B2182B"))(100)
gaps_col <- cumsum(c(length(genes_conserved),
                     length(genes_antigen),
                     length(genes_cytokine),
                     length(genes_chemokine),
                     length(genes_complement)))

ph <- pheatmap(mat,
               scale = "row", 
               color = pal,
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               show_rownames = TRUE,
               show_colnames = TRUE,
               annotation_col = anno_col,
               gaps_col = gaps_col,
               border_color = NA,
               legend_breaks = c(-2,0,2),
               legend_labels = c("-2","0","2"))
ggsave("results/cross_species/FigE_macrophage_conserved_function_heatmap.pdf",
       ph, width = 9.5, height = 3.5)
