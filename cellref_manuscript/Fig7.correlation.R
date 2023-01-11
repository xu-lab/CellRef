library(Seurat)
library(pheatmap)
library(ggplot2)

ihlca = readRDS(file=".Data/HLCA_v1.rds")
hlcellref = readRDS(file=".Data/LungMAP_HumanLung_CellRef.v1.rds")

# Create pseudobulk profiles

DefaultAssay(ihlca) = "RNA"
DefaultAssay(hlcellref) = "RNA"

ihlca = SetIdent(ihlca, value=ihlca@meta.data$ann_finest_level)
ihlca.avg = AverageExpression(ihlca, assay="RNA", return.seurat = T)
ihlca.avg = FindVariableFeatures(ihlca.avg, nfeatures = 2000)

hlcellref = SetIdent(hlcellref, value=hlcellref@meta.data$celltype_level3)
cellref.avg = AverageExpression(hlcellref, assay="RNA", return.seurat = T)
cellref.avg = FindVariableFeatures(cellref.avg, nfeatures = 2000)

# identify highly variable genes (HVG)

hvg.common = union(ihlca.avg@assays$RNA@var.features,
                   cellref.avg@assays$RNA@var.features)

hvg.common = hvg.common[which(hvg.common %in% rownames(cellref.avg@assays$RNA@data))]
hvg.common = hvg.common[which(hvg.common %in% rownames(ihlca.avg@assays$RNA@data))]

# scale the expression of HVG

ihlca.avg = ScaleData(ihlca.avg, features = hvg.common)
cellref.avg = ScaleData(cellref.avg, features = hvg.common)

ihlca.data = ihlca.avg@assays$RNA@scale.data
cellref.data = cellref.avg@assays$RNA@scale.data

colnames(ihlca.data) = paste0(colnames(ihlca.data), ".iHLCA")
colnames(cellref.data) = paste0(colnames(cellref.data), ".CellRef")

tmp = cbind(ihlca.data, cellref.data[rownames(ihlca.data), ])
any(is.na(tmp))

# calcualte correlations

viz = cor(tmp)

# visualize correlations and hierarchical clustering of pseudobulk profiles

pheatmap::pheatmap(viz, clustering_method = "average",
                   filename = "Figure_7A.tiff", width=20, height=19)

# save source data

sourcedata = list(hvg.common=hvg.common, correlations=viz)
save(sourcedata, file="Figure_7A.sourcedata.rda")

viz = sourcedata$correlations
write.table(viz, file="Figure_7A.correlations.txt", sep=",", col.names = T, row.names = T, quote = F)


#sessionInfo
sink("sessionInfo.txt")
sessionInfo()
sink()
