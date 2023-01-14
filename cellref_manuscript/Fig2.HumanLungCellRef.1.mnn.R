library(Seurat)
library(monocle3)
library(SingleR)
library(SingleCellExperiment)
library(harmony)
library(RobustRankAggreg)
library(dplyr)
library(reshape2)
library(ggplot2)
library(pheatmap)

options(future.globals.maxSize = 1280000 * 1024^2)

source("./CellRef_functions.R")

obj = readRDS(file=".Data/full.object.rds")

obj = subset(obj, Cohort != "GSE136831")
obj@meta.data = droplevels(obj@meta.data)

obj = NormalizeData(obj)
obj = CellCycleScoring(obj, s.features = Seurat::cc.genes.updated.2019$s.genes, g2m.features = Seurat::cc.genes.updated.2019$g2m.genes)


print(obj)

# batch correction of data from different donors
obj = doDataIntegration(obj, integration.batch ="DonorID",
                        method="Monocle3-mnn",
                        npcs=200, umap.min_dist = 0.1,  do.clustering = T)

g1 = DimPlot(obj, reduction="umap",
             group.by="integrated_clusters", label=T)
g2 = DimPlot(obj, reduction="umap",
             group.by="DonorID", label=F)
g2 = g2 + theme(legend.text = element_text(size=8))
g = g1+g2

ggsave(file="mnn.tiff", plot=g, width=16, height=8, dpi=300, units="in", compression="lzw")

saveRDS(obj, file="mnn.rds")
