# With the consistenly predicted cell types, we integrated the data from different donors using RPCA, identified 20 nearest neighbors for each cell in the integration,
# pruned cells with less than 60% kNN purity, and obtained the final set of cells with cell type annotations for the LungMAP Human Lung CellRef.
# The final data were integrated using RPCA for visualization and reference mapping applications.

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SingleCellExperiment)

options(future.globals.maxSize = 480000 * 1024^2)


#################################################################################################################################################
#################################################################################################################################################

batch_var = "DonorID"
ident_var = "cellref_type"

# integrated merged data using RPCA for finding nearest neighbors of each cell

obj = readRDS(file="merged.rds")

obj = CreateSeuratObject(counts=obj@assays$RNA@counts, meta.data=obj@meta.data)
obj = NormalizeData(obj)


# As the number of batches (104 donors) is large for the RPCA integration, we selected a set of donors as "reference" in the RPCA integration

# Sort donors based on the number of predicted cell types. Cell types with at least 9 cells were included in the calculation.
donor_celltype = table(obj@meta.data$DonorID, obj@meta.data[, ident_var])
donor_celltype = sort(rowSums(donor_celltype>9))
donor_celltype = donor_celltype[which(donor_celltype>9)]    # donors with at least 10 cell types

# Get the Dataset information of each donor
donor_dataset = obj@meta.data[, c("DonorID","Dataset")]
donor_dataset = unique(donor_dataset)
rownames(donor_dataset) = donor_dataset$DonorID

# For each dataset, we selected the donors with top 2 highest number of cell types predicted
donor_celltype = data.frame(donor_celltype)
colnames(donor_celltype) = "Types"
donor_celltype$Dataset = donor_dataset[rownames(donor_celltype),"Dataset"]
donor_celltype = donor_celltype[order(donor_celltype$Dataset, -donor_celltype$Types), ]
donor_celltype$Donor = rownames(donor_celltype)

refid = donor_celltype %>% group_by(Dataset) %>% slice_max(order_by = Types, n = 2)
reference.donors = as.character(refid$Donor)


# Perform SCT based integration usign RPCA

DefaultAssay(obj) = "RNA"
# split the dataset into a list of seurat objects, one for the data from a donor
objlist <- SplitObject(obj, split.by = batch_var)
# the positions of selected reference donors in the list
ref.ids = which(names(objlist) %in% reference.donors)
objlist <- lapply(X = objlist, FUN = function(x) {
    x <- SCTransform(x, vars.to.regress = c("S.Score","G2M.Score","pMT"), verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = objlist, nfeatures = 2000)
objlist <- PrepSCTIntegration(object.list = objlist, anchor.features = features)
objlist <- lapply(X = objlist, FUN = function(x) {
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Find integration anchors
obj.anchors <- FindIntegrationAnchors(object.list = objlist, anchor.features = features, reduction = 'rpca',
                                      normalization.method = "SCT",
                                      k.anchor = 10, l2.norm = TRUE,
                                      reference = ref.ids)

# prune anchors connecting two cells with different cell type predictions
anchors.df = obj.anchors@anchors
anchors.df$type1 = NA
anchors.df$type2 = NA
for (i in 1:length(obj.anchors@object.list)) {
  i.idx = i.type = NULL
  i.idx = which(anchors.df$dataset1 == i)
  if (length(i.idx)>0){
    i.type = as.character(objlist[[i]]@meta.data[as.numeric(anchors.df$cell1[i.idx]), ident_var])
    anchors.df$type1[i.idx] = i.type
  }
}
for (i in 1:length(obj.anchors@object.list)) {
  i.idx = i.type = NULL
  i.idx = which(anchors.df$dataset2 == i)
  if (length(i.idx)>0){
    i.type = as.character(objlist[[i]]@meta.data[as.numeric(anchors.df$cell2[i.idx]), ident_var])
    anchors.df$type2[i.idx] = i.type
  }
}
anchors.df$type_matched = (anchors.df$type1 ==  anchors.df$type2)
anchors.df = droplevels(subset(anchors.df, type_matched==TRUE))

obj.anchors1 = obj.anchors
obj.anchors1@anchors$type1 = anchors.df$type1
obj.anchors1@anchors$type2 = anchors.df$type2
obj.anchors1@anchors$type_matched = anchors.df$type_matched
obj.anchors1@anchors = anchors.df

# Integrate donor from different donors
combined = NULL
combined <- IntegrateData(anchorset = obj.anchors1, normalization.method = "SCT", k.weight=100)

#################################################################################################################################################
#################################################################################################################################################

DefaultAssay(combined) <- "integrated"
combined <- RunPCA(combined, npcs = npcs, verbose = FALSE)

# Find 20 nearest neighbors for each celll
combined = FindNeighbors(combined, reduction = "pca", dims=1:200, k.param = 20, nn.method="annoy", annoy.metric="cosine")

# for each cell, calculate the purity of cell type prediction in the neighbors of the cell
nn = combined@graphs$integrated_nn
purity = data.frame(cell=rownames(nn), ident=obj@meta.data[rownames(nn), ident.var])
purity$score = 0
rownames(purity) = as.character(purity$cell)

purity_helper <- function(x) {
  x.cell = as.character(x["cell"])
  x.ident = as.character(x["ident"])
  x.nn = names(which(nn[x.cell, ]==1))
  x.nn.types = purity[x.nn, "ident"]
  return(length(which(x.nn.types==x.ident)))
}
score = apply(purity, 1, FUN=function(x) purity_helper(x))
purity$score = score

obj@meta.data$nn.purity = nn.purity[rownames(obj@meta.data), "score"]/20


# removed cells with purity less than 60%. Seed cells were not removed.
cells.selected = obj@meta.data
cells.selected$selected = FALSE
cells.selected$selected[which(cells.selected$nn.purity>=0.6)] = TRUE
cells.selected$selected[which(cells.selected$id=="Seed")] = TRUE

ells.selected = subset(cells.selected, selected==TRUE)

# use the pruned data to construct the final Seurat object for the LungMAP Human Lung CellRef
counts.selected = obj@assays$RNA@counts[, rownames(cells.selected)]

obj = CreateSeuratObject(counts=counts.selected, meta.data=cells.selected)
obj = NormalizeData(obj)

#################################################################################################################################################
#################################################################################################################################################

# For visualization and automated cell type annotation, we used RPCA to integrate the final set of cells from different donors
DefaultAssay(obj) = "RNA"
# split the dataset into a list of seurat objects, one for the data from a donor
objlist <- SplitObject(obj, split.by = batch_var)
# the positions of selected reference donors in the list
ref.ids = which(names(objlist) %in% reference.donors)
objlist <- lapply(X = objlist, FUN = function(x) {
  x <- SCTransform(x, vars.to.regress = c("S.Score","G2M.Score","pMT"), verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = objlist, nfeatures = 2000)
objlist <- PrepSCTIntegration(object.list = objlist, anchor.features = features)
objlist <- lapply(X = objlist, FUN = function(x) {
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Find integration anchors
obj.anchors <- FindIntegrationAnchors(object.list = objlist, anchor.features = features, reduction = 'rpca',
                                      normalization.method = "SCT",
                                      k.anchor = 10, l2.norm = TRUE,
                                      reference = ref.ids)

# prune anchors linking cells of different cell type predictions
anchors.df = obj.anchors@anchors
anchors.df$type1 = NA
anchors.df$type2 = NA
for (i in 1:length(obj.anchors@object.list)) {
  i.idx = i.type = NULL
  i.idx = which(anchors.df$dataset1 == i)
  if (length(i.idx)>0){
    i.type = as.character(objlist[[i]]@meta.data[as.numeric(anchors.df$cell1[i.idx]), ident_var])
    anchors.df$type1[i.idx] = i.type
  }
}
for (i in 1:length(obj.anchors@object.list)) {
  i.idx = i.type = NULL
  i.idx = which(anchors.df$dataset2 == i)
  if (length(i.idx)>0){
    i.type = as.character(objlist[[i]]@meta.data[as.numeric(anchors.df$cell2[i.idx]), ident_var])
    anchors.df$type2[i.idx] = i.type
  }
}
anchors.df$type_matched = (anchors.df$type1 ==  anchors.df$type2)
anchors.df = droplevels(subset(anchors.df, type_matched==TRUE))

obj.anchors1 = obj.anchors
obj.anchors1@anchors$type1 = anchors.df$type1
obj.anchors1@anchors$type2 = anchors.df$type2
obj.anchors1@anchors$type_matched = anchors.df$type_matched
obj.anchors1@anchors = anchors.df

# Integrate donor from different donors
combined = NULL
combined <- IntegrateData(anchorset = obj.anchors1, normalization.method = "SCT", k.weight=100)

DefaultAssay(combined) <- "integrated"
combined <- RunPCA(combined, npcs = npcs, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:200, min.dist=0.5, return.model = T, n.neighbors = 50, seed.use=123L)

p1 <- DimPlot(combined, reduction = "umap", group.by = batch_var, pt.size=0.001) + NoLegend()
p2 <- DimPlot(combined, reduction = "umap", group.by = ident_var,label = TRUE, label.size=3, repel = TRUE, pt.size=0.001) + NoLegend()
g = p1 + p2
ggsave(file=paste0("LungMAP_HumanLung_CellRef_v1_umap.tiff"), plot=g, width=15, height=7.25, dpi=300, units="in", compression="lzw")


saveRDS(combined, file="LungMAP_HumanLung_CellRef.v1.rds")

