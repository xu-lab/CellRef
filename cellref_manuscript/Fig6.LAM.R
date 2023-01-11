library(ggplot2)
library(Seurat)
library(patchwork)
library(mclust)


options(future.globals.maxSize = 198000 * 1024^2) # for 50 Gb RAM

createColorPanel <- function(num.color){
  colPanel = c(
    "#E31A1C", "#FFD700", "#771122", "#777711", "#1F78B4", "#68228B", "#AAAA44",
    "#60CC52", "#771155", "#DDDD77", "#774411", "#AA7744", "#AA4455", "#117744",
    "#000080", "#44AA77", "#AA4488", "#DDAA77", "#D9D9D9", "#BC80BD", "#FFED6F",
    "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17",
    "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
    "#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
    "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928", "#FBB4AE", "#B3CDE3",
    "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2",
    "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC",
    "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FFFF33", "#A65628",
    "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
    "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
    "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5"
  )
  if(num.color > length(colPanel)){
    colPanel = c(colPanel, colVector(num.color - length(colPanel)));
  }else{
    colPanel = colPanel[1:num.color];
  }
  return(colPanel)
}

color_panel <- createColorPanel(50)

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


# color legend
types_original = names(table(obj@meta.data$celltype_original))
types_predicted = names(table(obj@meta.data$predicted.celltype))

types_original_colors = color_panel[1:length(types_original)]
types_predicted_colors = color_panel[1:length(types_predicted)]

names(types_original_colors) = types_original
names(types_predicted_colors) = types_predicted

types_original_colors["LAM+"] = rgb(255,0,255,maxColorValue = 255)


g1 = DimPlot(obj, reduction = "umap", group.by="celltype_original",
             label=F, label.size = 2, repel = T,
             cols=types_original_colors,
             pt.size = 0.001)  + NoAxes() + theme(plot.title = element_blank())

g2 = DimPlot(obj.sub, reduction = "umap",
             group.by="predicted.celltype.pruned", cols=types_predicted_colors,
             label=F, label.size = 2, repel = T, pt.size = 0.001)  + NoAxes() + theme(plot.title = element_blank())
g = g1+g2
ggsave(file="Fig6.LAM.AB.tiff", width=16, height=5, dpi=300, units="in", compression="lzw")


g = ggplot(obj@meta.data, aes(x=celltype_original, y=predicted.celltype.score, fill=celltype_original))
g = g + scale_fill_manual(values=types_original_colors)
g = g + geom_boxplot(outlier.shape = NA, alpha=0.5, position=position_dodge2())
g = g + geom_hline(yintercept = score.mean, color="black")
g = g + geom_hline(yintercept = score.cutoff, color="red")
g = g + ylab("Prediction score") + xlab("Original cell type annotation")  + guides(fill="none")
g = g + theme_bw() + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, color="black", size=14),
                           axis.title = element_text(size=12, color="black"))
g

ggsave(file="Fig6.LAM.D.tiff", width=9.5, height=4, dpi=300, units="in", compression="lzw")
