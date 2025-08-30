suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

tonsil <- readRDS("results/tonsil/tonsil.rds")

DefaultAssay(tonsil) <- "RNA"

heatmap_gene <- c("CD3D","CD4","LEF1","SELL",
                  "S100A4","TNFRSF25",
                  "FOXP3","TNFRSF18","LAG3",
                  "CTLA4","ICOS","TOX",
                  "CD8A","GZMK","NKG7","IL23R",
                  "CD79A","CD19","TNFRSF13C","CD40","CXCR4",
                  "CR2","TNFRSF13B","CD83","PAX5",
                  "MKI67","TOP2A",
                  "PCLAF","DUT",
                  "BCL6","POLD4","LY86",
                  "JCHAIN","FKBP11",
                  "KRT16","KRT4",
                  "MUC5B","TFF3","NUPR1",
                  "CD14","S100A12","S100A9",
                  "CD68","CD86","S100A4",
                  "CD163","C1QA")


heatmap_gene <- unique(as.character(heatmap_gene))

celltype <- c("Naïve CD4 T",
              "Memory CD4 T",
              "Treg CD4 T",
              "Exhausted CD4 T",
              "CD8 NKT",
              "NK",
              "Follicular B",
              "Marginal zone B",
              "Activated Follicular B",
              "Germinal center B",
              "Memory B",
              "Plasma B",
              "Epithelial",
              "Secretory cells",
              "DCs",
              "Type_1 Macrophage",
              "Type_2 Macrophage")

celltype <-rev(unique(as.character(celltype)))  
tonsil@meta.data$customclassif <- factor(tonsil@meta.data$customclassif, levels = celltype)

dot.data=p[["data"]]
colnames(dot.data)

dot.data$id <- factor(dot.data$id, levels = celltype)

p4=ggplot(dot.data,aes(x=features.plot, y = id, color = avg.exp.scaled, size = pct.exp)) + 
  geom_point() +
  scale_size(breaks = c(0,25,75),labels = c(0,25,75),range = c(0, 6), name = '%exp.cells')+
  cowplot::theme_cowplot() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust =1, hjust=1.0,colour="black",size = 10), 
        axis.text.y = element_text(colour="black",size = 10),
        legend.margin=margin(50, 0, 0, 0),
        panel.grid.major=element_line(colour="gray60"),
        panel.grid.minor=element_line(colour="gray40"),
        panel.border = element_rect(color = "gray40", size = 1, fill = NA)) +  
  ylab('') +
  xlab('') +  
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = c("#252CBD","#F8FAEC","#FFD097","#AA2526"),breaks = c(-1,0,1,2),limits = c(-1,2), oob = scales::squish, name = 'expression') 
p4

ggsave("Tonsil.maker.tiff", p4,width = 14, height = 5)
ggsave("Tonsil.maker.pdf", p4,width = 14, height = 5)
dev.off() 
