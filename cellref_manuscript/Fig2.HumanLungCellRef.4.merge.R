# merge the label transfer results from Seurat refmap and SingleR, prune the predictions with low confidence, prune the predictions with inconsistent annotations in the two methods

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)

options(future.globals.maxSize = 128000 * 1024^2) # for 50 Gb RAM

# load objects from Seurat refmap
objlist = list()
datasets = c("CCHMC_LungMAP", "EGAS00001004082", "GSE122960", "GSE134174", "GSE135893", "GSE136831", "UPenn_LungMAP", "syn21041850", "GSE161382", "GSE171524")
for (i in 1:length(datasets)) {
  i.dataset = datasets[i]
  i.obj = readRDS(file=paste0(i.dataset, ".mapped.rds"))

  i.obj@assays = i.obj@assays['RNA']
  i.obj$id = i.dataset
  DefaultAssay(i.obj) <- 'RNA'
  i.obj@meta.data = i.obj@meta.data[, meta.use]

  objlist[[i]] = i.obj
}
names(objlist) = datasets

# load LungMAP CellRef Seed
Seed = readRDS(file=".Data/LungMAP_HumanLung_CellRef_Seed.v1.rds")
Seed$predicted.celltype = Seed@meta.data$celltype_level3
Seed$predicted.celltype.score = 1
Seed@assays = Seed@assays['RNA']
Seed$id = "Seed"
DefaultAssay(Seed) <- 'RNA'

# merge all objects
obj = Seed
for (i in 1:length(objlist)) {
  i.obj = objlist[[i]]
  obj = merge(obj, i.obj)
}
DefaultAssay(obj) = "RNA"


# merge refmap and singleR predictions
cells = subset(obj@meta.data, id != "Seed")
cells.selected = NULL

# for the Seurat 4 refmap predictions, for each dataset, top 90% of cells with highest prediction scores were kept
datasets = names(table(cells$Dataset))
for (i in 1:length(datasets)) {
  i.cells = subset(cells, Dataset==datasets[i])
  i.cells = i.cells[order(-i.cells$predicted.celltype.score), ]

  i.value = quantile(i.cells$predicted.celltype.score, probs = 0.1)
  i.cells.selected = subset(i.cells, predicted.celltype.score>=i.value)

  if (is.null(cells.selected)) {
    cells.selected = i.cells.selected
  } else {
    cells.selected = rbind(cells.selected, i.cells.selected)
  }
}
cells.seed = subset(obj@meta.data, id == "Seed")
cells.selected = rbind(cells.selected, cells.seed)

obj@meta.data$pred_RefMap = NA
obj@meta.data[rownames(cells.selected), "pred_RefMap"] = as.character(cells.selected$predicted.celltype)

# load singleR predictions
pred.singleR = readRDS(file="prediction.SingleR.rds")
pred.singleR = data.frame(pred.singleR)
obj@meta.data$pred_SingleR = as.character(obj@meta.data$predicted.celltype)
obj@meta.data[rownames(pred.singleR), "pred_SingleR"] = as.character(pred.singleR$pruned.labels)

# Since the Seurat refmap didn't perform well in distinguishing Basal and Suprabasal annotations, we temporarily combined Basal and Suprabasal predictions
obj@meta.data$pred_RefMap[which(obj@meta.data$pred_RefMap %in% c("Basal","Suprabasal"))] = "Basal_Suprabasal"
obj@meta.data$pred_SingleR[which(obj@meta.data$pred_SingleR %in% c("Basal","Suprabasal"))] = "Basal_Suprabasal"

# Remove label transfer results that are not consistent in the two methods
obj@meta.data$cellref_type = ifelse(as.character(obj@meta.data$pred_RefMap) == as.character(obj@meta.data$pred_SingleR),
                                    obj@meta.data$pred_RefMap, "Unclassified")
obj@meta.data$cellref_type[which(is.na(obj@meta.data$cellref_type))] = "Unclassified"
obj = subset(obj, cellref_type != "Unclassified")
obj@meta.data = droplevels(obj@meta.data)

# The predictions of Basal and Suprabasal annotations were based on the results from SingleR
cells.use = intersect(rownames(obj@meta.data), rownames(pred.singleR))
pred.singleR = pred.singleR[cells.use, ]
obj@meta.data[rownames(pred.singleR), "cellref_type"] = as.character(pred.singleR$pruned.labels)
obj@meta.data[rownames(cells.seed), "cellref_type"] = as.character(cells.seed$predicted.celltype)

# Save the merged results
saveRDS(obj, file="merged.rds")
