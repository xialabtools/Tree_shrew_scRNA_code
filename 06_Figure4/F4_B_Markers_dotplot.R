suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(cowplot); library(scales)
})

OUTDIR <- "./data/infection.spleen"
F4DIR  <- file.path(OUTDIR, "Fig4")
dir.create(F4DIR, showWarnings = FALSE, recursive = TRUE)

grad4 <- c("#252CBD","#F8FAEC","#FFD097","#AA2526")

scRNA <- readRDS(file.path(OUTDIR, "scRNA.finall.rds"))
DefaultAssay(scRNA) <- "RNA"; scRNA <- NormalizeData(scRNA, verbose = FALSE)

## marker
markers_f4B <- c("CD3D","CD4","CD8A","GZMK","IFNG","CD79A","MKI67","JCHAIN","CSF3R","FCGR3B",
                 "CSF1R","PPARG","NR1H3","FLT3","IRF8")

## DotPlot
dp <- DotPlot(scRNA, features = markers_f4B, group.by = "celltype") + coord_flip()
dd <- dp[["data"]]
dd$id <- factor(dd$id, levels = rev(unique(as.character(dd$id))))

pB <- ggplot(dd, aes(x = features.plot, y = id, color = avg.exp.scaled, size = pct.exp)) +
  geom_point() +
  scale_size(breaks = c(0,25,75), labels = c(0,25,75), range = c(0,6), name = '%exp.cells') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, colour="black", size=10),
        axis.text.y = element_text(colour="black", size=10),
        panel.grid.major=element_line(colour="gray60"),
        panel.grid.minor=element_line(colour="gray40"),
        panel.border = element_rect(color="gray40", size=1, fill=NA),
        axis.ticks = element_blank()) +
  xlab(NULL) + ylab(NULL) +
  scale_color_gradientn(colours = grad4, breaks = c(-1,0,1,2), limits = c(-1,2),
                        oob = squish, name = "expression")

ggsave(file.path(F4DIR, "F4B_markers_dotplot.pdf"), pB, width=6.8, height=5.0)
