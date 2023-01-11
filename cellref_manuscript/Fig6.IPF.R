library(Seurat)
library(ggplot2)
library(readxl)
library(dplyr)
options(future.globals.maxSize= 1024*1024^2)

cellref.seed = readRDS(file = '.Data/LungMAP_HumanLung_CellRef_Seed.v1.rds')
ILD = readRDS(file = './.Data/GSE135893_ILD_annotated_fullsize.rds')

ipf = subset(ILD, Diagnosis == 'IPF')
# Dataset Info
table(ipf$Sample_Name)
length(table(ipf$Sample_Name)) #12

# Map IPF to CellRef ------------------------------------------------------
#Combine CD4_T, CD8_T, Treg as T Cells
#Combine iMon, pMon as Monocytes
cellref.seed$celltype_level2[which(cellref.seed$celltype_level2 %in% c('CD4_T','CD8_T', 'Treg'))] = 'T Cells'
cellref.seed$celltype_level2[which(cellref.seed$celltype_level2 %in% c('pMON','iMON'))] = 'Monocytes'

DefaultAssay(ipf) = 'RNA'
ipf = NormalizeData(ipf)
DefaultAssay(ipf) = 'SCT'

anchors = FindTransferAnchors(
  reference = cellref.seed,
  query = ipf,
  normalization.method = 'SCT',
  reference.reduction = 'pca',
  dims = 1:dim(cellref.seed@reductions$pca@cell.embeddings)[2])

mapped = MapQuery(
  anchorset = anchors,
  query = ipf,
  reference = cellref.seed,
  refdata = list(celltype_level = 'celltype_level2'),
  reference.reduction = 'pca',
  reduction.model = 'umap'
)

ipf@meta.data$predicted.celltype_level = mapped@meta.data[rownames(ipf@meta.data),'predicted.celltype_level']
ipf@meta.data$predicted.celltype_level.score = mapped@meta.data[rownames(ipf@meta.data), 'predicted.celltype_level.score']


# Final Plots -------------------------------------------------------------
# Color Scale
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
sub.colors <- c("magenta", createColorPanel(50))

g1 = DimPlot(ipf, reduction = 'umap', group.by = 'predicted.celltype_level', label = F, repel = T, pt.size = 0.001, label.size = 3,cols = sub.colors) + ggtitle('CellRef Seed Annotation')+ theme(legend.text = element_text(size=15))
ggsave('./Figure_6E.tiff', dpi = 600, units = 'in', compression = 'lzw', height = 8, width = 12)

# Prediction Score Featureplots
FeaturePlot(ipf, reduction = "umap", features = "predicted.celltype_level.score",pt.size = 0.1, order = F) + 
  scale_colour_gradient2(low = 'red',mid = 'grey85', high = "grey85", midpoint = mean(ipf$predicted.celltype_level.score)- sd(ipf$predicted.celltype_level.score))+
  ggtitle('Prediction Score')
ggsave('./Figure_6F.tiff', dpi = 600, units = 'in', compression = 'lzw', height = 8, width = 8.5)


g = ipf@meta.data %>% ggplot(aes(x=celltype, y=predicted.celltype_level.score, fill=celltype)) + geom_boxplot(outlier.shape = NA) 
g = g + theme_bw() 
g = g + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size = 16, colour = 'black'), axis.text.y = element_text(size = 16, colour = 'black'), axis.title = element_text(size=16, color = 'black',family = 'TT Arial'))
g = g + labs(x = '', y = 'Prediction Score')
g = g + geom_hline(aes(yintercept = mean(ipf$predicted.celltype_level.score)), color='black')
g = g + geom_hline(aes(yintercept = c(mean(ipf$predicted.celltype_level.score) - sd(ipf$predicted.celltype_level.score))),color='red',size=1)
g= g + ylim(0, NA)
g = g + NoLegend()
ggsave(filename = './Figure_6G_right_score_boxplot.tiff', dpi = 600, units = 'in', compression = 'lzw', height = 7, width = 12)

#UMAP and FeaturePlots

g1 = DimPlot(ipf, reduction = 'umap', group.by = 'celltype', label = T, repel = T, pt.size = 0.001, label.size = 3) + ggtitle('Original Annotation') + theme(legend.text = element_text(size=15))
ggsave('./Figure_6G_left_umap.tiff', dpi = 600, units = 'in', compression = 'lzw', height = 8, width = 14)



# SessionInfo -------------------------------------------------------------
sink('sessionInfo.txt')
sessionInfo()
sink()