# transfer the cell type annotations of LungMAP_HumanLung_CellRef_Seed to all other scRNA and snRNA-seq data using SingleR

library(Seurat)
library(ggplot2)
library(patchwork)
library(SingleR)

scRNA = readRDS(file=".Data/scRNA.rds")

seed = readRDS(file=".Data/LungMAP_HumanLung_CellRef_Seed.v1.rds")
DefaultAssay(seed) = "RNA"
seed = NormalizeData(seed)
seed = FindVariableFeatures(seed)
seed = ScaleData(seed, vars.to.regress = c("S.Score","G2M.Score","pMT"))

seed.sce = as.SingleCellExperiment(seed)


predictions_singleR = NULL


for (dataset in c("GSE161382","GSE171524", "GSE136831", "CCHMC_LungMAP", "EGAS00001004082",
                  "GSE135893", "GSE122960", "syn21041850", "GSE134174", "UPenn_LungMAP")) {

  if (dataset %in% c("GSE161382","GSE171524")) {
    query = readRDS(file=paste0("../.Data/snRNA.", dataset, ".rds"))

  } else {
    query = subset(scRNA, Dataset==dataset)
    query@meta.data = droplevels(query@meta.data)
  }

  cat("In total:", dim(query@meta.data)[1], "cells in the Dataset", dataset, "\n")

 # exclude cells that are in the seed
  cells_in_seed = query@meta.data[which( rownames(query@meta.data) %in% rownames(seed@meta.data) ), ]
  cat(dim(cells_in_seed)[1], "cells in the CellRef seed\n")

  cells_new = query@meta.data[which(!(rownames(query@meta.data) %in% rownames(seed@meta.data))), ]

  cat("Mapping the remaining", dim(cells_new)[1], "cells\n")

  query = subset(query, cells=rownames(cells_new))

  query@meta.data = droplevels(query@meta.data)

  query = NormalizeData(query)
  query = FindVariableFeatures(query)
  query = ScaleData(query, vars.to.regress = c("S.Score","G2M.Score","pMT"))

  query.sce = as.SingleCellExperiment(query)

  pred.singleR.celltype <- SingleR(test=query.sce, ref=seed.sce,
                                   labels=as.character(seed.sce@colData$celltype_level3), de.method="wilcox")
  pred.singleR.celltype$Cell = rownames(pred.singleR.celltype)
  pred.singleR.celltype$Dataset = dataset

  if (is.null(predictions_singleR)) {
    predictions_singleR = pred.singleR.celltype
  } else {
    predictions_singleR = rbind(predictions_singleR, pred.singleR.celltype)
  }
}

saveRDS(predictions_singleR, file="prediction.SingleR.rds")
