library(ggplot2)
library(Seurat)
library(patchwork)
library(mclust)

# load the CellRef Seed
ref = readRDS(file=".Data/LungMAP_HumanLung_CellRef_Seed.v1.rds")

# annotation using CellRef seed
query = readRDS(file=".Data/LAM.scRNAseq.rds")
DefaultAssay(query) = "RNA"
query = SCTransform(query, vars.to.regress=c("G2M.Score","S.Score","pMT"))

anchors <- FindTransferAnchors(
  reference = ref,
  query = query,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:dim(ref@reductions$pca@cell.embeddings)[2]
)

obj <- MapQuery(anchorset = anchors, reference = ref, query = query,
                   refdata = list(celltype = "celltype_level2"),
                   reference.reduction = "pca", reduction.model = "umap")

score.mean = mean(obj@meta.data$predicted.celltype.score)
score.sd = sd(obj@meta.data$predicted.celltype.score)
score.cutoff = score.mean-score.sd

obj@meta.data$predicted.celltype.pruned = as.character(obj@meta.data$predicted.celltype)
obj@meta.data$predicted.celltype.pruned[which(obj@meta.data$predicted.celltype.score < score.cutoff)] = "Pruned"

# saveRDS(obj, file="./Fig6.LAM.A.object.rds")

# exclude cells with prediction score >= the cutoff
obj.sub = subset(obj, predicted.celltype.pruned != "Pruned")

# exclude singleton cell type predictions
obj.sub = subset(obj.sub, predicted.celltype.pruned != "AF2")
obj.sub = subset(obj.sub, predicted.celltype.pruned != "CD4_T")
obj.sub = subset(obj.sub, predicted.celltype.pruned != "Pericyte")

obj.sub@meta.data = droplevels(obj.sub@meta.data)

# saveRDS(obj.sub, file="./Fig6.LAM.B.object.rds")

g1 = DimPlot(obj, reduction = "umap", group.by="celltype_original",
             label=F, label.size = 2, repel = T,
             pt.size = 0.001)  + NoAxes() + theme(plot.title = element_blank())
g2 = DimPlot(obj.sub, reduction = "umap",
             group.by="predicted.celltype.pruned", 
             label=F, label.size = 2, repel = T, pt.size = 0.001)  + NoAxes() + theme(plot.title = element_blank())
g = g1+g2
ggsave(file="Fig6.LAM.AB.tiff", width=16, height=5, dpi=300, units="in", compression="lzw")

g = ggplot(obj@meta.data, aes(x=celltype_original, y=predicted.celltype.score, fill=celltype_original))
g = g + geom_boxplot(outlier.shape = NA, alpha=0.5, position=position_dodge2())
g = g + geom_hline(yintercept = score.mean, color="black")
g = g + geom_hline(yintercept = score.cutoff, color="red")
g = g + ylab("Prediction score") + xlab("Original cell type annotation")  + guides(fill="none")
g = g + theme_bw() + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, color="black", size=14),
                           axis.title = element_text(size=12, color="black"))

ggsave(file="Fig6.LAM.D.tiff", width=9.5, height=4, dpi=300, units="in", compression="lzw")
