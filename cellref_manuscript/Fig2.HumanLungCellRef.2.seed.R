# code to reproduce the LungMAP Human Lung CellRef Seed identification

library(CellRef)
library(Seurat)
library(ggplot2)
library(RobustRankAggreg)
library(SingleR)
library(mclust)
library(BiocParallel)
library(BiocNeighbors)

options(future.globals.maxSize = 1280000 * 1024^2)

# load the MNN based data integration,
# the "Cluster" column in the metadata contains Leiden based cell clustering
obj = readRDS(file=".Data/mnn.rds")

# as shown in the vigenette, we used findCandidateClusters function followed by manual curation
# to identify candidate cell clusters in the integrated data for each of the 48 cell types in the dictionary (Supplementary Data 2 of the manuscript).
# here, we load the data frame that contains the cell cluster and cell type mapping
candidate_clusters = readRDS(file=".Data/mnn_cluster_celltype_mapping.rds")

# load the cell type dictionary for the 48 cell types (Supplementary Data 2 of the manuscript)
# For cell types that shared the same candidate cluster in the above mapping, including ASMC and SCMF, AF2 and Chondrocyte, cDC1 and maDC, 
# we added the positive markers of one cell type (e.g., ASMC) as the negative markers of the cell type (i.e., SCMF) that shared the cell cluster, and vice versa, for the seed cell identification.
ctd = readRDS(file=".Data/HumanLung_CellRef_dictionary.rds")

# Now we identified the seed cells based on the data integration, the candidate cell clusters, and the cell type markers
obj@meta.data$Cluster = paste0("C", obj@meta.data$clusters)
obj_seed = findSeedCells(obj, ctd, candidate_clusters,
                         cluster.var="Cluster",
                         score.thresh=Inf,
                         seed.n.min=-Inf, seed.n.max=200,
                         verbose = F)

saveRDS(obj_seed, file="./CellRef_seed.rds")
