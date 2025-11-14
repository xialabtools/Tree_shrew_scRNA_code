suppressPackageStartupMessages({
  library(Seurat); library(UCell); library(ggplot2); library(ggpubr)
})

OUTDIR <- "./data/infection.spleen"
F4DIR  <- file.path(OUTDIR, "Fig4")
dir.create(F4DIR, showWarnings = FALSE, recursive = TRUE)

scRNA <- readRDS(file.path(OUTDIR, "scRNA.finall.rds"))
Idents(scRNA) <- "celltype"
macro_NR1H3 <- subset(scRNA, idents = "NR1H3_Macro")
macro_IRF8  <- subset(scRNA, idents = "IRF8_Macro")

DefaultAssay(macro_NR1H3) <- "RNA"; macro_NR1H3 <- NormalizeData(macro_NR1H3, verbose=FALSE)
DefaultAssay(macro_IRF8)  <- "RNA"; macro_IRF8  <- NormalizeData(macro_IRF8,  verbose=FALSE)

markers <- list()
markers$negtive <- c("NR1H3","PPARG","SYT11")
markers$postive <- c("IRF8","HAVCR2","CARD11","FLT3","SULF2","CADM1","CD83","FCER1A","CLNK","CD226","LRRK2")

macro_NR1H3 <- AddModuleScore_UCell(macro_NR1H3, features = markers)
macro_IRF8  <- AddModuleScore_UCell(macro_IRF8,  features = markers)

df_neg <- FetchData(macro_NR1H3, vars = c("infection","negtive_UCell")); df_neg$infection <- factor(df_neg$infection, levels=c("0dpi","3dpi","7dpi"))
df_pos <- FetchData(macro_IRF8,  vars = c("infection","postive_UCell")); df_pos$infection <- factor(df_pos$infection, levels=c("0dpi","3dpi","7dpi"))

cmp <- list(c("0dpi","3dpi"), c("0dpi","7dpi"), c("3dpi","7dpi"))

p_neg <- ggplot(df_neg, aes(x=infection, y=negtive_UCell, fill=infection)) +
  geom_violin(color='black', size=0.5, trim=FALSE) +
  geom_boxplot(width=0.12, fill="white", outlier.shape=NA) +
  theme_classic() + theme(legend.position="none") +
  stat_compare_means(method="wilcox.test", hide.ns=TRUE, comparisons=cmp, label="p.format", bracket.size=0.5, size=4) +
  ylim(0.25, 1.25) + ylab("anti-inflammatory (UCell)") + xlab(NULL)

p_pos <- ggplot(df_pos, aes(x=infection, y=postive_UCell, fill=infection)) +
  geom_violin(color='black', size=0.5, trim=FALSE) +
  geom_boxplot(width=0.12, fill="white", outlier.shape=NA) +
  theme_classic() + theme(legend.position="none") +
  stat_compare_means(method="wilcox.test", hide.ns=FALSE, comparisons=cmp, label="p.format", bracket.size=0.5, size=4) +
  ylim(0.25, 1.50) + ylab("pro-inflammatory (UCell)") + xlab(NULL)

ggsave(file.path(F4DIR, "F4G_UCell_violin.pdf"), p_neg / p_pos, width=4.8, height=5.6)
