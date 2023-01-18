# transfer the cell type annotations of LungMAP_HumanLung_CellRef_Seed to all other scRNA and snRNA-seq data using Seurat 4 reference mapping

library(Seurat)
library(ggplot2)
library(patchwork)

scRNA = readRDS(file=".Data/scRNA.rds")

seed = readRDS(file=".Data/LungMAP_HumanLung_CellRef_Seed.v1.rds")

predictions_refmap = NULL

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
  cat(dim(cells_in_seed)[1], "cells in the reference core\n")

  cells_new = query@meta.data[which(!(rownames(query@meta.data) %in% rownames(seed@meta.data))), ]

  cat("Mapping the remaining", dim(cells_new)[1], "cells\n")

  query = subset(query, cells=rownames(cells_new))

  query@meta.data = droplevels(query@meta.data)

  query = SCTransform(query, vars.to.regress = c("S.Score","G2M.Score","pMT"))

  anchors <- FindTransferAnchors(
    reference = seed,
    query = query,
    normalization.method = "SCT",
    reference.reduction = "pca",
    dims = 1:dim(seed@reductions$pca@cell.embeddings)[2]
  )

  mapped <- MapQuery(anchorset = anchors, reference = seed, query = query,
                     refdata = list(celltype = "celltype_level3"), reference.reduction = "pca", reduction.model = "umap")

  saveRDS(query, file=paste0(dataset, ".mapped.rds"))


  if (is.null(predictions_refmap)) {
    predictions_refmap = mapped@meta.data
  } else {
    predictions_refmap = rbind(predictions_refmap, mapped@meta.data)
  }
}

saveRDS(predictions_refmap, file = "predictions_refmap.rds")
