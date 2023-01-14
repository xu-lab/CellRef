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

source("../CellRef_functions.R")

# load the mnn based data integration,
# the "Cluster" column in the metadata contains Leiden based cell clustering
obj = readRDS(file=".Data/mnn.rds")

# curated candidate cell clusters for the 48 cell types
candidate_clusters = readRDS(file=".Data/mnn_cluster_celltype_mapping.rds")

# cell type dictionary for the 48 cell types
# In our curated candidate cell cluster and cell type mapping above,
# there are some cell types shared the same candidate cluster.
# To identify better seed cells for those cell types, we added the positive markers of one cell type
# as negative markers of the other cell type that share the same candidate cluster, and vice visa.
# We did this dictionary extension for ASMC and SCMF, AF2 and Chondrocyte, and cDC1 and maDC.
ctd = readRDS(file=".Data/HumanLung_CellRef_dictionary.rds")

# Now we identified the seed cells based on the data integration, the candidate cell clusters, and the cell type markers
obj@meta.data$Cluster = paste0("C", obj@meta.data$clusters)
obj_seed = findSeedCells(obj, ctd, candidate_clusters,
                         cluster.var="Cluster",
                         score.thresh=Inf,
                         seed.n.min=-Inf, seed.n.max=200,
                         verbose = F)

saveRDS(obj_seed, file="./CellRef_seed.rds")
