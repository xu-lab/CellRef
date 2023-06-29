#' Perform batch correction of single cell RNA-seq data
#'
#' This function loads a Seurat object that contains all data for batch correction and integration.
#' It will use the "RNA" assay to perform the data integration and assume that the "data" slot contains the normalized data.
#' Currently, the function supports three single cell batch correction methods
#' Monocle3-mnn: the mutual nearest neighbor matching based algorithm (Haghverdi et al., 2018) implemented in Monocle 3 (Cao et al., 2019)
#' Seurat4-rpca: the reciprocal PCA based integration in Seurat 4 (Hao and Hao, et al., 2021)
#' Harmony: the Harmony based integration (Korsunsky et al., 2019)
#'
#' For large data integration tasks, we suggest running the integration in a cluster environment.
#' For example, the data integration for LungMAP Human Lung CellRef was performed in a linux cluster.
#'
#' @param object A Seurat object
#' @param integration.batch a character string indicating the meta.data column that contains the data batch for correction
#' @param ident.var a character string indicating the meta.data column that contains the identities of each cell if available
#' @param method the method that will be used for single cell batch correction
#' @param npcs the number of reduced dimensions used for batch correction and umap projection
#' @param vars.to.regression Variables to regress out. For example, cell cycle (S.Score, G2M.Score) and mitochondrial percentage (pMT).
#' @param umap.min_dist The effective minimum distance between embedded points, controlling how tightly the embedding is allowed to compress points together.
#' @param umap.n_neighbors The size of local neighborhood (in terms of number of neighboring sample points) used for manifold approximation
#' @param do.clustering Whether to perform clustering analysis after batch correction
#' @param clustering.method The clustering method. Default is leiden algorithm.
#' @param clustering.resolution The resolution of graph-based clustering
#' @param rpca.normalization The normalization method used in Seurat's RPCA pipeline
#' @param rpca.k.anchor The k.anchor parameter used in Seurat's RPCA pipeline
#' @param rpca.k.weight The k.weight parameter used in Seurat's RPCA pipeline
#' @param rpca.prune_anchors When performing RPCA based integration, whether to prune anchors connecting cells of different identities. Need to set the ident.var parameter.
#' @param verbose Whether to show verbose output
#' @return A Seurat object
#' @export
doDataIntegration <- function(object, integration.batch, ident.var=NULL, method=c("Monocle3-mnn","Seurat4-rpca", "Harmony"),
                          npcs=200, vars.to.regress=c("S.Score","G2M.Score","pMT"),
                          umap.min_dist = 0.3, umap.n_neighbors = 50,
                          do.clustering=TRUE, clustering.method="leiden", clustering.resolution=NULL,
                          rpca.normalization="SCT", rpca.k.anchor=10, rpca.k.weight=50, rpca.prune_anchors=T,
                          verbose=F, ...) {

  #TODO: add input validation
  if(!(integration.batch %in% colnames(object@meta.data))) {
    stop("The integration.batch was not found in the meta.data")
  }

  if (!is.null(ident.var)) {
    if(!(ident.var %in% colnames(object@meta.data))) {
      stop("The ident.var was not found in the meta.data")
    }
  }

  if (method == "Monocle3-mnn") {

    object = doDataIntegration_mnn(object=object, integration.batch=integration.batch, ident.var=ident.var,
                                   npcs=npcs, vars.to.regress=vars.to.regress,
                                   umap.min_dist = umap.min_dist, umap.n_neighbors = umap.n_neighbors,
                                   do.clustering=do.clustering, clustering.method=clustering.method, clustering.resolution=clustering.resolution,
                                   rpca.normalization=rpca.normalization, rpca.k.anchor=rpca.k.anchor,
                                   rpca.k.weight=rpca.k.weight, rpca.prune_anchors=rpca.prune_anchors,
                                   verbose=verbose, ...)

  } else if (method=="Seurat4-rpca") {

    object = doDataIntegration_rpca(object=object, integration.batch=integration.batch, ident.var=ident.var,
                                   npcs=npcs, vars.to.regress=vars.to.regress,
                                   umap.min_dist = umap.min_dist, umap.n_neighbors = umap.n_neighbors,
                                   do.clustering=do.clustering, clustering.method=clustering.method, clustering.resolution=clustering.resolution,
                                   rpca.normalization=rpca.normalization, rpca.k.anchor=rpca.k.anchor,
                                   rpca.k.weight=rpca.k.weight, rpca.prune_anchors=rpca.prune_anchors,
                                   verbose=verbose, ...)

  } else if (method=="Harmony") {

    object = doDataIntegration_harmony(object=object, integration.batch=integration.batch, ident.var=ident.var,
                                    npcs=npcs, vars.to.regress=vars.to.regress,
                                    umap.min_dist = umap.min_dist, umap.n_neighbors = umap.n_neighbors,
                                    do.clustering=do.clustering, clustering.method=clustering.method, clustering.resolution=clustering.resolution,
                                    rpca.normalization=rpca.normalization, rpca.k.anchor=rpca.k.anchor,
                                    rpca.k.weight=rpca.k.weight, rpca.prune_anchors=rpca.prune_anchors,
                                    verbose=verbose, ...)
  }

  return(object)
}


doDataIntegration_mnn <- function(object, integration.batch, ident.var=NULL,
                                  npcs=200, vars.to.regress=c("S.Score","G2M.Score","pMT"),
                                  umap.min_dist = 0.3, umap.n_neighbors = 50,
                                  do.clustering=TRUE, clustering.method="leiden", clustering.resolution=NULL,
                                  rpca.normalization="SCT", rpca.k.anchor=10, rpca.k.weight=100, rpca.prune_anchors=T,
                                  verbose=F, ...) {
  if(!(integration.batch %in% colnames(object@meta.data))) {
    stop("The integration.batch was not found in the meta.data")
  }

  if (!is.null(ident.var)) {
    if(!(ident.var %in% colnames(object@meta.data))) {
      stop("The ident.var was not found in the meta.data")
    }
  }

  gene_annotation = data.frame(gene_short_name = rownames(object@assays$RNA@counts),
                               gene_id=rownames(object@assays$RNA@counts))
  rownames(gene_annotation) =  as.character(gene_annotation$gene_id)

  cds <- new_cell_data_set(object@assays$RNA@counts,
                           cell_metadata = object@meta.data,
                           gene_metadata = gene_annotation)

  cds <- preprocess_cds(cds, method="PCA", num_dim = npcs)
  if (length(vars.to.regress)>0) {
    formula_str = paste0("~", paste0(vars.to.regress, collapse = "+"))
    cds <- align_cds(cds, preprocess_method = "PCA", alignment_group = integration.batch, residual_model_formula_str=formula_str)
  } else {
    cds <- align_cds(cds, preprocess_method = "PCA", alignment_group = integration.batch)
  }

  cds <- reduce_dimension(cds, preprocess_method = "Aligned",
                          umap.min_dist = umap.min_dist,
                          umap.n_neighbors = umap.n_neighbors)

  #dims.pca = reducedDims(cds)$PCA
  #object@reductions$pca = CreateDimReducObject(embeddings = dims.pca, assay = "RNA", key="pca_")

  dims.aligned = reducedDims(cds)$Aligned
  object@reductions$pca = CreateDimReducObject(embeddings = dims.aligned, assay = "RNA", key="pca_")

  dims.umap = reducedDims(cds)$UMAP
  colnames(dims.umap) = c("umap_1", "umap_2")

  object@reductions$umap = CreateDimReducObject(embeddings = dims.umap, assay = "RNA", key="umap_")


  if (do.clustering==TRUE) {
    cds = cluster_cells(cds, cluster_method = clustering.method, clustering.resolution=clustering.resolution)

    clusters = data.frame(cluster=as.numeric(cds@clusters$UMAP$clusters), cell=names(cds@clusters$UMAP$clusters), stringsAsFactors = F)
    rownames(clusters) = clusters$cell

    object@meta.data$integrated_clusters = clusters[rownames(object@meta.data), "cluster"]
  }

  return(object)
}

doDataIntegration_rpca <- function(object, integration.batch, ident.var=NULL,
                                   npcs=200, vars.to.regress=c("S.Score","G2M.Score","pMT"),
                                   umap.min_dist = 0.3, umap.n_neighbors = 50,
                                   do.clustering=TRUE, clustering.method="leiden", clustering.resolution=NULL,
                                   rpca.normalization="SCT", rpca.k.anchor=10, rpca.k.weight=100, rpca.prune_anchors=T,
                                   verbose=F, ...) {

  if(!(integration.batch %in% colnames(object@meta.data))) {
    stop("The integration.batch was not found in the meta.data")
  }

  if (!is.null(ident.var)) {
    if(!(ident.var %in% colnames(object@meta.data))) {
      stop("The ident.var was not found in the meta.data")
    }
  }

  DefaultAssay(object) = "RNA"

  # split the dataset into a list of two seurat objects (stim and CTRL)
  objlist <- SplitObject(object, split.by = integration.batch)

  if (rpca.normalization=="LogNormalize") {
    # normalize and identify variable features for each dataset independently
    objlist <- lapply(X = objlist, FUN = function(x) {
      x <- NormalizeData(x, verbose=F)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    })

    # select features that are repeatedly variable across datasets for integration
    # run PCA on each dataset using these features
    features <- SelectIntegrationFeatures(object.list = objlist, verbose = F)

    objlist <- lapply(X = objlist, FUN = function(x) {

      if (length(vars.to.regress)>0) {
        x <- ScaleData(x, features = features, vars.to.regress = vars.to.regress, verbose = FALSE)
      } else {
        x <- ScaleData(x, features = features, verbose = FALSE)
      }

      x <- RunPCA(x, features = features, verbose = FALSE)
    })
  } else if (rpca.normalization=="SCT") {

    # normalize and identify variable features for each dataset independently
    objlist <- lapply(X = objlist, FUN = function(x) {

      if (length(vars.to.regress)>0) {
        x <- SCTransform(x, vars.to.regress = vars.to.regress, verbose = FALSE)
      } else {
        x <- SCTransform(x, verbose = FALSE)
      }

    })
    features <- SelectIntegrationFeatures(object.list = objlist, nfeatures = 2000, verbose=F)
    objlist <- PrepSCTIntegration(object.list = objlist, anchor.features = features, verbose=F)

    objlist <- lapply(X = objlist, FUN = function(x) {
      x <- RunPCA(x, features = features, verbose = FALSE)
    })
  }

  # find anchors
  obj.anchors <- FindIntegrationAnchors(object.list = objlist, anchor.features = features, reduction = 'rpca',
                                        normalization.method = rpca.normalization,
                                        k.anchor = rpca.k.anchor, l2.norm = TRUE, verbose = verbose)

  if (rpca.prune_anchors & !is.null(ident.var)) {

    anchors.df = obj.anchors@anchors
    anchors.df$type1 = NA
    anchors.df$type2 = NA

    for (i in 1:length(obj.anchors@object.list)) {
      i.idx = i.type = NULL
      i.idx = which(anchors.df$dataset1 == i)
      if (length(i.idx)>0){
        i.type = as.character(objlist[[i]]@meta.data[as.numeric(anchors.df$cell1[i.idx]), ident.var])
        anchors.df$type1[i.idx] = i.type
      }
    }

    for (i in 1:length(obj.anchors@object.list)) {
      i.idx = i.type = NULL
      i.idx = which(anchors.df$dataset2 == i)
      if (length(i.idx)>0){
        i.type = as.character(objlist[[i]]@meta.data[as.numeric(anchors.df$cell2[i.idx]), ident.var])
        anchors.df$type2[i.idx] = i.type
      }
    }

    anchors.df$type_matched = (anchors.df$type1 ==  anchors.df$type2)

    obj.anchors1 = obj.anchors
    all(paste0(obj.anchors1@anchors$cell1,obj.anchors1@anchors$cell2) == paste0(anchors.df$cell1, anchors.df$cell2))

    obj.anchors1@anchors$type1 = anchors.df$type1
    obj.anchors1@anchors$type2 = anchors.df$type2
    obj.anchors1@anchors$type_matched = anchors.df$type_matched

    # 14142 -> 13174
    anchors.df = droplevels(subset(anchors.df, type_matched==TRUE))
    obj.anchors1@anchors = anchors.df

    obj.anchors = obj.anchors1
  }

  object <- IntegrateData(anchorset = obj.anchors, normalization.method = rpca.normalization, k.weight=rpca.k.weight, verbose = verbose)


  # specify that we will perform downstream analysis on the corrected data
  # note that the original unmodified data still resides in the 'RNA' assay
  DefaultAssay(object) <- "integrated"

  # Run the standard workflow for visualization and clustering
  if (rpca.normalization=="LogNormalize") {
    object <- ScaleData(object, verbose = FALSE)
  }

  object <- RunPCA(object, npcs = npcs, verbose = FALSE)
  object <- RunUMAP(object, reduction = "pca", dims = 1:npcs, return.model = T,
                    min.dist=umap.min_dist, n.neighbors = umap.n_neighbors, verbose = verbose)

  if (do.clustering==TRUE) {
    object = FindNeighbors(object, reduction = "pca", dims=1:npcs, verbose=verbose)
    if (is.null(clustering.resolution)) clustering.resolution=0.8
    if (clustering.method == "leiden") {
      object = FindClusters(object, graph.name = "integrated_nn", method = 4, resolution = clustering.resolution, verbose=verbose)
    } else {
      object = FindClusters(object, graph.name = "integrated_nn", method = 3, resolution = clustering.resolution, verbose=verbose)

    }
    object@meta.data$integrated_clusters = object@meta.data$seurat_clusters
  }

  return(object)

}

doDataIntegration_harmony <- function(object, integration.batch, ident.var=NULL,
                                      npcs=200, vars.to.regress=c("S.Score","G2M.Score","pMT"),
                                      umap.min_dist = 0.3, umap.n_neighbors = 50,
                                      do.clustering=TRUE, clustering.method="leiden", clustering.resolution=NULL,
                                      rpca.normalization="SCT", rpca.k.anchor=10, rpca.k.weight=50, rpca.prune_anchors=T,
                                      verbose=F, ...) {

  if(!(integration.batch %in% colnames(object@meta.data))) {
    stop("The integration.batch was not found in the meta.data")
  }

  if (!is.null(ident.var)) {
    if(!(ident.var %in% colnames(object@meta.data))) {
      stop("The ident.var was not found in the meta.data")
    }
  }

  DefaultAssay(object) = "RNA"
  object = NormalizeData(object)
  object = FindVariableFeatures(object)
  object = ScaleData(object, vars.to.regress=vars.to.regress)

  object = RunPCA(object, npcs=npcs)

  hmat = HarmonyMatrix(object@reductions$pca@cell.embeddings[, 1:npcs],
                       meta_data = object@meta.data, vars_use = integration.batch, do_pca = FALSE)
  object@reductions$harmony = CreateDimReducObject(embeddings = hmat, assay= DefaultAssay(object), key="harmony_")
  object = RunUMAP(object, dims=1:npcs, reduction = "harmony",
                   reduction.name = "umap", return.model = TRUE,
                   min.dist=umap.min_dist, n.neighbors = umap.n_neighbors)

  if (do.clustering==TRUE) {
    object = FindNeighbors(object, reduction = "pca", dims=1:npcs, verbose=verbose)
    if (is.null(clustering.resolution)) clustering.resolution=0.8
    if (clustering.method == "leiden") {
      object = FindClusters(object, method = 4, resolution = clustering.resolution, verbose=verbose)
    } else {
      object = FindClusters(object, method = 3, resolution = clustering.resolution, verbose=verbose)

    }
    object@meta.data$integrated_clusters = object@meta.data$seurat_clusters
  }

  return(object)
}

#' Identify candidate cell clusters for each cell type defined in the dictionary
#'
#' @param object A Seurat object
#' @param ctd A data.frame object that contains the cell type dictionary. The data frame should contain at least three columns: CellType, Marker, and MarkerType (p - positive markers; n - negative markers).
#' @param cluster.var The name of meta.data column that contains cluster assignment of each cell. By default, "integrated_clusters".
#' @param expr.thresh A marker gene is expressed in a cell if its expression value >expr.thresh. By default, 0.
#' @param use.scaled.expr If TRUE, use the scaled expression in the "scale.data" slot; otherwise, use the expression in the "data" slot. By default, TRUE.
#' @param score.thresh Cells with aggregated rank scores >=score.thresh will be excluded.
#' @param seed.n.min Minimum number of seed cells for a cell type
#' @param seed.n.max Maximum number of seed cells for a cell type
#' @param group.size.min Minimum size of a cell cluster to be included in the candidate cell cluster identification
#' @param precision.thresh Precision threshold for determining candidate cell clusters
#' @param recall.thresh Recall threshold for determining candidate cell clusters
#' @param verbose Whether to show verbose output
#' @return a data.frame containing the mapping of candidate cell clusters (Cluster) and cell types (CellType).
#' @export
findCandidateClusters <- function(object, ctd, cluster.var="integrated_clusters",
                                  expr.thresh = 0, use.scaled.expr=T,
                                  score.thresh=0.1, seed.n.min=5, seed.n.max=Inf,
                                  group.size.min=10, precision.thresh=0.05, recall.thresh=0.2, verbose=T) {
  # TODO input validation

  candidate_clusters=NULL
  seed_types = names(table(ctd$CellType))

  celldist = table(object@meta.data[, cluster.var])

  DefaultAssay(object) = "RNA"
  object = ScaleData(object, features = as.character(ctd$Marker), verbose=F)

  for (i in 1:length(seed_types)) {

    i.type = seed_types[i]

    if (verbose) {
      message(paste("\nFind candidate cell clusters for cell type:", i.type, "\n"))
    }

    # get positive and negative markers
    i.markers = subset(ctd, CellType==i.type)

    i.markers.p = which(i.markers$MarkerType=="p")
    i.markers.n = which(i.markers$MarkerType=="n")

    if (length(i.markers.p)>1) { # require at lest 2 positive markers

      i.candidates.n = c()

      # cells with negative marker expression
      if (length(i.markers.n)>0) {
        for (tt in 1:length(i.markers.n)) {
          tt.marker = i.markers$Marker[i.markers.n[tt]]
          i.candidates.n = union(i.candidates.n, which(object@assays$RNA@counts[as.character(tt.marker), ]>0))
        }
      }

      # remove negative cells from all cells
      i.candidates = 1:dim(object@meta.data)[1]
      i.candidates = setdiff(i.candidates,  i.candidates.n)

      if (length(i.candidates)>=seed.n.min) {

        # positive marker expression
        i.x = NULL
        if (use.scaled.expr==T) {
          i.x = as.matrix(object@assays$RNA@scale.data[i.markers[i.markers.p, "Marker"], i.candidates])
        } else {
          i.x = as.matrix(object@assays$RNA@data[i.markers[i.markers.p, "Marker"], i.candidates])
        }

        i.x = i.x[, which(colSums(i.x>expr.thresh)>1)]

        i.cells = c()
        i.list <- list()
        for (j in 1:dim(i.x)[1]) {
          i.j.cells = NULL
          i.j.cells <- i.x[j,]
          i.j.cells <- i.j.cells[i.j.cells > expr.thresh]
          if(length(i.j.cells) > 0) {
            i.j.cells <- sort(i.j.cells, decreasing=T)
            i.list[[as.character(rownames(i.x)[j])]] <- names(i.j.cells)
            i.cells = union(i.cells, names(i.j.cells))
          }
        }

        if (verbose) {
          message(paste(length(i.cells), "candidate cells\n"))
        }

        i.rank <- aggregateRanks(glist = i.list, N = length(i.cells))
        i.rank = i.rank[order(i.rank$Score, decreasing = F), ]

        i.cells = subset(i.rank, Score<score.thresh)

        if (!is.null(dim(i.cells)) & dim(i.cells)[1]>1) {
          if (verbose) {
            message(paste(dim(i.cells)[1], " cells with Score <", score.thresh, "\n"))
          }

          i.cells = object@meta.data[rownames(i.cells), ]

          i.groups = data.frame(table(i.cells[, cluster.var]))
          colnames(i.groups)= c("Cluster","Cell_count")

          i.groups$Recall = i.groups$Cell_count/sum(i.groups$Cell_count)
          i.groups$Cell_Total = celldist[as.character(i.groups$Cluster)]

          i.groups$Precision = i.groups$Cell_count/i.groups$Cell_Total

          i.groups = i.groups[order(-i.groups$Cell_count), ]

          if (verbose) {
            message(paste("\nKeep groups with at least", group.size.min, "cells"))
          }
          i.groups = subset(i.groups, Cell_count>group.size.min)

          if (verbose) {
            message(paste("keep groups with recall>", recall.thresh))
          }
          i.groups.1 = subset(i.groups, Recall>=recall.thresh)

          if (verbose) {
            message(paste("keep groups with precision>", precision.thresh))
          }
          i.groups.2 = subset(i.groups, Precision>=precision.thresh)

          i.groups.selected = unique(rbind(i.groups.1, i.groups.2))
          rownames(i.groups.selected) = NULL

          if (dim(i.groups.selected)[1]>0) {

            i.mapping=data.frame(i.groups.selected, CellType=i.type)

            if (is.null(candidate_clusters)) {
              candidate_clusters = i.mapping
            } else {
              candidate_clusters = rbind(candidate_clusters, i.mapping)
            }
          }
        }

      }

    } else {

      if (verbose) {
        message("less than 2 markers found in the data. Skip prediction.")
      }

    }
  }

  if (!is.null(candidate_clusters)) {
    candidate_clusters = candidate_clusters[, c("Cluster", "CellType")]
  }

  return(candidate_clusters)

}

#' Seed cell identification
#'
#' @param object A Seurat object
#' @param ctd A data.frame object that contains the cell type dictionary. The data frame should contain at least three columns: CellType, Marker, and MarkerType (p - positive markers; n - negative markers).
#' @param candidate.clusters A data.frame that contains the mapping between candidate cell clusters and cell types. Column names should be "Cluster and "CellType".
#' @param cluster.var The name of meta.data column that contains cluster assignment of each cell. By default, "integrated_clusters".
#' @param expr.thresh A marker gene is expressed in a cell if its expression value >expr.thresh. By default, 0.
#' @param use.scaled.expr If TRUE, use the scaled expression in the "scale.data" slot; otherwise, use the expression in the "data" slot. By default, TRUE.
#' @param score.thresh Cells with aggregated rank scores >=score.thresh will be excluded.
#' @param seed.n.min Minimum number of seed cells for a cell type
#' @param seed.n.max Maximum number of seed cells for a cell type
#' @param verbose Whether to show verbose output
#' @return A Seurat object containing the identified seed cells
#' @export
findSeedCells <- function(object, ctd, candidate.clusters,
                          cluster.var="integrated_clusters",
                          expr.thresh = 0, use.scaled.expr=T,  score.thresh=Inf,
                          seed.n.min=5, seed.n.max=200, verbose=T) {
  #TODO input validation

  seed_types = names(table(ctd$CellType))

  object = ScaleData(object, features = as.character(ctd$Marker), verbose = verbose)
  object@meta.data$Seed = NA

  ranklist = list()

  stats = data.frame(Type=seed_types, nCell_neg=0, nCell_cluster=0, nCell_cluster_noneg=0, nCell_pos=0, nCell_seed=0)

  for (i in 1:length(seed_types)) {

    i.type = seed_types[i]
    if (verbose) {
      message(paste0("\n", i.type, "\n"))
    }

    i.markers = subset(ctd, CellType==i.type)

    i.markers.p = which(i.markers$MarkerType=="p")
    i.markers.n = which(i.markers$MarkerType=="n")

    if (length(i.markers.p)>1) {

      i.candidates.n = c()
      if (length(i.markers.n)>0) {
        for (tt in 1:length(i.markers.n)) {
          tt.marker = i.markers$Marker[i.markers.n[tt]]
          i.candidates.n = union(i.candidates.n, which(object@assays$RNA@counts[as.character(tt.marker), ]>0))
        }

        if (verbose) {
          message(paste(length(i.candidates.n), "cells expressing at least 1 negative marker\n"))
        }

        stats$nCell_neg[i] = length(i.candidates.n)
      } else {

        if (verbose) {
          message("no negative markers used\n")
        }
      }

      i.clusters = subset(candidate.clusters, CellType==i.type)$Cluster

      i.candidates = c()
      i.candidates = which(object@meta.data[, cluster.var] %in% i.clusters)

      if (verbose) {
        message(paste(length(i.candidates), "candidate cells\n"))
      }


      stats$nCell_cluster[i] = length(i.candidates)

      i.candidates = setdiff(i.candidates,  i.candidates.n)

      stats$nCell_cluster_noneg[i] = length(i.candidates)

      if (length(i.candidates)>=seed.n.min) {

        i.x = NULL
        if (use.scaled.expr==T) {
          i.x = as.matrix(object@assays$RNA@scale.data[i.markers[i.markers.p, "Marker"], i.candidates])
        } else {
          i.x = as.matrix(object@assays$RNA@data[i.markers[i.markers.p, "Marker"], i.candidates])
        }

        # cells with at least two markers
        i.x = i.x[, which(colSums(i.x>expr.thresh)>1)]

        i.cells = c()
        i.list <- list()
        for (j in 1:dim(i.x)[1]) {
          i.j.cells = NULL
          i.j.cells <- i.x[j,]
          i.j.cells <- i.j.cells[i.j.cells > expr.thresh]
          if(length(i.j.cells) > 0) {
            i.j.cells <- sort(i.j.cells, decreasing=T)
            i.list[[as.character(rownames(i.x)[j])]] <- names(i.j.cells)
            i.cells = union(i.cells, names(i.j.cells))
          }
        }

        if (verbose) {
          message(paste(length(i.cells), "candidate cells with at least 2 positive marker expression\n"))
        }

        stats$nCell_pos[i] = length(i.cells)

        if (length(i.cells)> seed.n.min) {

          i.rank <- aggregateRanks(glist = i.list, N = length(i.cells))
          i.rank = i.rank[order(i.rank$Score, decreasing = F), ]

          ranklist[[i.type]] = i.rank

          i.rank.sig = subset(i.rank, Score<score.thresh)

          i.idx = NULL
          i.idx = which(i.rank$Score<score.thresh)

          if (verbose) {
            message(paste(length(i.idx), "cells have Score<", score.thresh, "\n"))
          }

          i.seed = NULL

          if (length(i.idx) >= seed.n.min) {

            i.seed = i.rank[i.idx, ]

            if (length(i.idx) > seed.n.max) {

              if (verbose) {
                message(paste("More than the maximum", seed.n.max, "threshold. Select top", seed.n.max, "cells as seed\n"))
              }

              i.seed = i.seed[1:seed.n.max, ]
            }

            object@meta.data[rownames(i.seed), "Seed"] = i.type

            stats$nCell_seed[i] = dim(i.seed)[1]

          } else {

            if (verbose) {
              message(paste(i.type,": skipped. Less than", seed.n.min, "cells\n"))
            }

            stats$nCell_seed[i] = dim(i.seed)[1]
          }

        } else {


          if (verbose) {
            message(paste(i.type,": skipped. Less than", seed.n.min, "cells\n"))
          }

          stats$nCell_seed[i] = dim(i.seed)[1]
        }

      }

    } else {

      if (verbose) {
        message("less than 2 markers found in the data. Skip prediction.")
      }

    }
  }

  object = subset(object, cells=which(!is.na(object@meta.data$Seed)))

  return(object)
}

####


annotateCells_RefMap <- function(ref, query, query.batch=NULL, reference.ident="Seed",
                                 normalization.method="SCT", reference.reduction="pca",
                                 reference.model = "umap",
                                 do.prune=T, prune.prob=0.1,
                                 verbose=T) {

  predictions = NULL

  if (is.null(query.batch)) {
    query.batch = paste0("batch", gsub("-", "", Sys.Date()))
    query@meta.data[, query.batch] = "All"
  }

  batches = names(table(query@meta.data[, query.batch]))

  DefaultAssay(query) = "RNA"

  for (batch in batches) {

    query1 = subset(query, cells=rownames(query@meta.data[which(query@meta.data[, query.batch]==batch), ]))

    if (verbose) {
      message(paste("\n\n", batch, ":", dim(query1@meta.data)[1], "cells"))
    }

    cells_in_ref = query1@meta.data[which( rownames(query1@meta.data) %in% rownames(ref@meta.data) ), ]

    if (verbose) {
      message(paste(dim(cells_in_ref)[1], "cells in the reference"))
    }


    cells_new = query1@meta.data[which(!(rownames(query1@meta.data) %in% rownames(ref@meta.data))), ]

    if (verbose) {
      message(paste("Mapping the remaining", dim(cells_new)[1], "cells\n"))
    }

    query1 = subset(query1, cells=rownames(cells_new))

    query1@meta.data = droplevels(query1@meta.data)

    anchors <- FindTransferAnchors(
      reference = ref,
      query = query1,
      normalization.method = normalization.method,
      reference.reduction = reference.reduction,
      dims = 1:dim(ref@reductions[[reference.reduction]]@cell.embeddings)[2],
      verbose=verbose
    )

    query1 <- MapQuery(anchorset = anchors, reference = ref, query = query1,
                       refdata = list(celltype = reference.ident),
                       reference.reduction = reference.reduction, reduction.model = reference.model,
                       verbose=verbose)

    i.annotation =
    if (is.null(predictions)) {
      predictions = query1@meta.data
    } else {
      predictions = rbind(predictions, query1@meta.data)
    }
  }


  predictions$pred_RefMap = as.character(predictions$predicted.celltype)
  if(do.prune==TRUE) {
    score.thresh = data.frame(predictions %>% group_by(Dataset) %>% summarize(thresh=quantile(predicted.celltype.score, probs=prune.prob)))
    for (i in 1:dim(score.thresh)[1]) {
      i.batch = score.thresh[i, query.batch]
      i.cells = which(predictions[, query.batch]==i.batch)
      i.cells = intersect(i.cells, which(predictions$predicted.celltype.score<score.thresh$thresh[i]))
      predictions$pred_RefMap[i.cells] = NA
    }
  }

  return(predictions)
}


annotateCells_SingleR <- function(ref, query, query.batch=NULL, reference.ident="Seed",
                                  vars.to.regress = c("S.Score","G2M.Score","pMT"),
                                  verbose=T,
                                  ...) {

  predictions = NULL

  DefaultAssay(ref) = "RNA"
  ref = NormalizeData(ref, verbose=F)
  ref = FindVariableFeatures(ref, verbose=F)

  if (length(vars.to.regress)>0) {
    ref = ScaleData(ref, vars.to.regress=vars.to.regress, verbose=F)
  } else {
    ref = ScaleData(ref, verbose=F)
  }


  ref.sce = as.SingleCellExperiment(ref)

  if (is.null(query.batch)) {
    query.batch = paste0("batch", gsub("-", "", Sys.Date()))
    query@meta.data[, query.batch] = "All"
  }

  batches = names(table(query@meta.data[, query.batch]))


  for (i in 1:length(batches)) {

    batch = batches[i]

    query1 = subset(query, cells=rownames(query@meta.data[which(query@meta.data[, query.batch]==batch), ]))

    if (verbose) {
      message(paste("\n\n", batch, ":", dim(query1@meta.data)[1], "cells"))
    }


    cells_in_ref = query1@meta.data[which( rownames(query1@meta.data) %in% rownames(ref@meta.data) ), ]

    if (verbose==T) {
      message(paste(dim(cells_in_ref)[1], "cells in the reference"))
    }

    cells_new = query1@meta.data[which(!(rownames(query1@meta.data) %in% rownames(ref@meta.data))), ]

    if (verbose==T) {
      message(paste("Mapping the remaining", dim(cells_new)[1], "cells\n"))
    }

    query1 = subset(query1, cells=rownames(cells_new))

    query1@meta.data = droplevels(query1@meta.data)

    DefaultAssay(query1) = "RNA"
    query1 = NormalizeData(query1, verbose=F)
    query1 = FindVariableFeatures(query1, verbose=F)

    if (length(vars.to.regress)>0) {
      query1 = ScaleData(query1, vars.to.regress=vars.to.regress, verbose=F)
    } else {
      query1 = ScaleData(query1, verbose=F)
    }

    query.sce = as.SingleCellExperiment(query1)

    pred.singleR.celltype <- SingleR(test=query.sce, ref=ref.sce,
                                     labels=as.character(ref.sce@colData[, reference.ident]),
                                     de.method="wilcox")
    pred.singleR.celltype$Cell = rownames(pred.singleR.celltype)
    pred.singleR.celltype[, query.batch] = batch

    if (is.null(predictions)) {
      predictions = pred.singleR.celltype
    } else {
      predictions = rbind(predictions, pred.singleR.celltype)
    }
  }

  predictions = data.frame(predictions)

  return(predictions)

}

#' Build a CellRef using the initial data integration and the identified seed cells
#'
#' @param object.all A Seurat object that contains all data for CellRef construction after batch correction
#' @param object.seed The Seurat object that contains the identified seed cells. The object should have run UMAP project and have stored the UMAP model.
#' @param mapping.batch a character string indicating the meta.data column that contains data batch information to perform mapping to the seed cells
#' @param integration.batch the name of the meta.data column that contains the batch information of CellRef cells for correction using Seurat's RPCA method.
#' @param integration.npcs the number of principal components
#' @param rpca.prune_anchors In Seurat's RPCA integration, whether to prune anchors connecting cells of different identities.
#' @param purity.pruning Whether to prune cells based on kNN purity
#' @param purity.nn.k The number of neighbors used in the kNN purity calculation
#' @param purity.thresh The threshold for kNN purity score
#' @param verbose Whether to show verbose output
#' @return A Seurat object
#' @export
buildCellRef <- function(object.all, object.seed,
                         mapping.batch="Dataset",
                         integration.batch="Dataset", integration.npcs=200,
                         rpca.prune_anchors=T,
                         purity.pruning=T, purity.nn.k=20, purity.thresh=0.6,
                         vars.to.regress = c("S.Score","G2M.Score","pMT"),
                         verbose=T) {
  # TODO: add input validation

  object =  NULL

  if (verbose) {
    message("Mapping all the cells to the seed cells using Seurat4's reference mapping algorithm")
  }
  preds = annotateCells_RefMap(ref=object.seed, query=object.all, query.batch=mapping.batch,
                               verbose=F)

  if (verbose) {
    message("Mapping all the cells to the seed cells using Seurat4's reference mapping algorithm")
  }
  pred_SingleR = annotateCells_SingleR(ref=object.seed, query=object.all, query.batch=mapping.batch,
                                       vars.to.regress = vars.to.regress,
                                       verbose=F)

  # add SingleR predictions to the Seurat predictions
  preds$pred_SingleR = pred_SingleR[rownames(preds), "pruned.labels"]

  # select cells with consistent predictions
  pred_consensus = subset(preds, pred_RefMap==pred_SingleR)

  if (verbose) {
    message(paste(dim(pred_consensus)[1], "have consistent cell type predictions"))
  }

  cellref.ann.name="cellref.celltype"

  pred_consensus$cellref.source = "consensus"
  pred_consensus[, cellref.ann.name] = as.character(pred_consensus$pred_RefMap)

  cellref_cells = object.seed@meta.data
  cellref_cells$cellref.source = "Seed"
  cellref_cells[, cellref.ann.name] = as.character(cellref_cells$Seed)
  cellref_cells = rbind(cellref_cells[, c("cellref.source",cellref.ann.name)], pred_consensus[, c("cellref.source",cellref.ann.name)])

  object = subset(object.all, cells=rownames(cellref_cells))

  object@meta.data$cellref.source = cellref_cells[rownames(object@meta.data), "cellref.source"]
  object@meta.data[, cellref.ann.name] = cellref_cells[rownames(object@meta.data), cellref.ann.name]

  if (verbose) {
    message(paste("Combined the predicted cells with", dim(object.seed@meta.data)[1], "seed cells"))
  }

  if (verbose) {
    message(paste("Integrate",dim(object@meta.data)[1], "cells using Seurat4's RPCA pipeline using SCT normalization"))
  }
  object = doDataIntegration(object, integration.batch=integration.batch, ident.var=cellref.ann.name, method="Seurat4-rpca",
                             npcs=integration.npcs, rpca.prune_anchors = T, vars.to.regress = vars.to.regress,
                             do.clustering=FALSE, verbose=F)

  # kNN
  if (purity.pruning) {

    if (verbose) {
      message(paste("Performing kNN purity pruning\nCalculating kNN purity using k=", purity.nn.k))
    }

    object = FindNeighbors(object, reduction = "pca", dims=1:integration.npcs, k.param = purity.nn.k, nn.method="annoy", annoy.metric="cosine",
                           graph.name="cellref.purity.knn",
                           verbose=verbose)

    nn=NULL
    nn = object@graphs[['cellref.purity.knn']]

    purity=NULL
    purity = data.frame(cell=rownames(nn), ident=object@meta.data[rownames(nn), cellref.ann.name])
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
    purity$score = purity$score/purity.nn.k

    if (verbose) {
      message(paste("Excluding", length(which(purity$score<purity.thresh)), "cells with kNN.purity<", purity.thresh))
    }
    object@meta.data$knn.purity = purity[rownames(object@meta.data), "score"]
    object = subset(object, knn.purity >= purity.thresh)

    if (verbose) {
      message(paste("Integrate",dim(object@meta.data)[1], "cells using Seurat4's RPCA pipeline using SCT normalization"))
    }
    object = doDataIntegration(object, integration.batch=integration.batch, ident.var=cellref.ann.name, method="Seurat4-rpca",
                               npcs=integration.npcs, rpca.prune_anchors = rpca.prune_anchors, do.clustering=FALSE,
                               vars.to.regress = vars.to.regress,
                               verbose=F)

  }

  return(object)
}


#' CellRef validation of automated cell type annotation
#'
#' This function loads a Seurat object with cell type annotations predicted by CellRef or CellRef seed.
#' It will produce a PDF, containing UMAP visualization of cells colored by the CellRef predicted cell types,
#' differential expression of CellRef markers in the predictions, dotplot visualization of expression of CellRef marker genes,
#' identification of differentially expressed genes (DEGs) in the predicted cell types,
#' heatmap visualization of DEGs, and most enriched KEGG pathways and Gene Ontology Biological Processes enriched by the DEGs.
#'
#' @param obj A Seurat object
#' @param markers A data frame that contains marker(column name: 'Marker') and cell type(column name: 'CellType') information,
#' @param ident.use A character string indicating the meta.data column that contains the identities of each cell
#' @param umap.width Width of output UMAP
#' @param umap.height Height of output UMAP
#' @param log.thresh The threshold of Log Fold Change when finding the expression of CellRef marker genes for given ident.use
#' @param min.pct The min frequency to find the expression of CellRef marker genes for given ident.use
#' @param min.cells.group The smallest number of cells for a cell type in your annotation
#' @param dotplot.height Height of output dotplot
#' @param dotplot.width Width of output dotplot
#' @param vlnplt.height Height of output violin plot
#' @param Vlnplt.width Width of output violin plot
#' @param fullsigs.log.thrsh The threshold of log.FC when finding full signature genes for given ident.use
#' @param fullsigs.min.pct The min frequency to find full signature genes for given ident
#' @param heatmap.height Height of output heatmap
#' @param heatmap.width Width of output heatmap
#' @param heatmap_n Number of DEGs per cell type in the heatmap
#' @param dotchart.height Height of dotchart
#' @param dotchart.width Width of dotchart
#' @param enrich.height Height of enrichment analysis plots
#' @param enrich.width Width of enrichment analysis plots
#' @return A PDF file and multiple data frames.
#' @export
#'
evaluatePrediction = function(obj,
                         markers,
                         ident.use,
                         umap.width = 10,
                         umap.height = 7,
                         log.thresh = log2(1.5),
                         min.pct = 0.2,
                         min.cells.group = 3,
                         dotplot.height = 13,
                         dotplot.width = 9,
                         vlnplt.height = 5,
                         vlnplt.width = 3,
                         fullsigs.log.thrsh = log2(1.5),
                         fullsigs.min.pct = 0.2,
                         heatmap.height = 10,
                         heatmap.width = 16,
                         heatmap_n = 50,
                         dotchart.height = 10,
                         dotchart.width = 16,
                         enrich.height = 6,
                         enrich.width = 8){

  #detect ident.use in meta.data
  if(!ident.use %in% colnames(obj@meta.data)){
    stop(paste0("Did not find <", ident.use, "> in object's meta data!"))
  }

  Idents(obj) = obj@meta.data[,ident.use]
  DefaultAssay(obj) = 'RNA'

  obj.celltypes = names(table(Idents(obj)))
  marker.celltypes = names(table(markers$CellType))

  #detect obj and marker consistency
  if (all(obj.celltypes %in% marker.celltypes) == F) {
    stop(paste0("The following cell types are inconsistent with marker table's cell type: ", obj.celltypes[!obj.celltypes %in% marker.celltypes]))
  }

  #UMAP: UMAP of the cell type predictions
  if(T){
    cat('Generating UMAP')
    umap_label = DimPlot(obj, reduction = "umap", group.by=ident.use, label=T, label.size=2, repel=T) + ggtitle('UMAP of Predicted Cell Types')
    ggsave(filename = "UMAP.tiff", width=umap.width, height=umap.height, dpi=300, units="in", compression="lzw")
  }
  # Dotplot : Get markers that remain DE in the predcited cell types
  if(T){
    cat('Generating dotplot')
    targets = names(table(Idents(obj)))
    tmp = markers

    df1 = data.frame(row.names = targets)
    df2 = data.frame(row.names = targets)

    for (x in 1:length(targets)) {
      tmp.markers = subset(tmp, CellType == targets[x])
      tmp.markers = tmp.markers$Marker
      tmp.markers = intersect(rownames(obj), tmp.markers)

      de = FindMarkers(obj,ident.1 = targets[x], assay = 'RNA',only.pos = T, logfc.threshold = log.thresh ,min.pct = min.pct,features = tmp.markers, min.cells.group = min.cells.group)
      de = subset(de, p_val<0.05)
      de$gene = rownames(de)
      de$cluster = targets[x]

      df1 = rbind(df1, de)
      write.table(df1, paste0('geneList.markers.remained.selective.txt'), sep = '\t', row.names = T)

      pct = length(de$gene) / length(tmp.markers)
      df2[targets[x], 'Selective.Percentage'] = pct
      write.table(df2, paste0('pct.markers.remained.selective.txt'), sep = '\t', row.names = T)
    }

    genes = df1$gene
    genes = genes[which(!duplicated(genes))]
    obj@meta.data[,ident.use] = factor(obj@meta.data[,ident.use], levels = targets)

    dotplot = DotPlot(obj, assay="RNA", features=rev(genes),  group.by= ident.use,
                      dot.min=0.01, cols=c("grey90",rgb(146,0,0,maxColorValue = 255)))  + coord_flip()
    dotplot = dotplot + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), axis.text.y=element_text(face="italic"), axis.title.y = element_blank(),
                              axis.title.x = element_blank(), panel.border = element_rect(color="black")) + NoLegend() + ggtitle('Markers that Remained Selective in Testing Data')
    ggsave(filename = paste0("dotplot.markers.remained.selective.tiff"), height = dotplot.height, width = dotplot.width, dpi = 300)
  }

  # Figure 5F: Percentage of genes remained selective in prediction
  if(T){
    cat('Generating violin plot')
    df = df2
    df$Selective.Percentage = round(df$Selective.Percentage*100,2)
    df$DataID = 'DataID'
    errbar_lims <- group_by(df, DataID) %>%
      summarize(mean=mean(Selective.Percentage), se=sd(Selective.Percentage)/sqrt(n()),
                upper=mean+se, lower=mean-se)

    vlnplt = ggplot() + geom_violin(data=df, aes(x=DataID, y=Selective.Percentage, fill=DataID),scale='width', width=0.8, trim = T) +
      geom_point(data=errbar_lims, aes(x=DataID, y=mean), size=3) +
      geom_errorbar(aes(x=errbar_lims$DataID, ymax=errbar_lims$upper,
                        ymin=errbar_lims$lower), stat='identity', width=.15)

    vlnplt = vlnplt + theme_bw()
    vlnplt = vlnplt + theme(panel.grid = element_blank())
    vlnplt = vlnplt + scale_fill_manual(values=c('grey'))
    vlnplt = vlnplt + labs(x = '', y= '% of marker genes of a cell type')
    vlnplt = vlnplt + theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 17, color =  'black'), axis.text.y.right = element_blank(), axis.title = element_text(size = 17, color = 'black'))
    vlnplt = vlnplt + coord_cartesian(ylim = c(0, 105))
    vlnplt = vlnplt + scale_y_continuous(breaks = seq(0, 100, 10)) + ggtitle('Summary of Percentage of Markers Remained Selective')
    ggsave('vlnplot.markers.remained.selective.tiff', width = vlnplt.width, height = vlnplt.height, units = 'in', compression='lzw', dpi =600)
  }

  # Figure 5G: Heatmap of DEGs for cell types in object
  if(T){
    cat('Generating heatmap...')
    sigs = FindAllMarkers(obj, logfc.threshold = fullsigs.log.thrsh, test.use = "wilcox", only.pos = T, min.pct = fullsigs.min.pct)
    sigs = sigs[order(sigs$cluster, -sigs$avg_log2FC), ]
    sigs1 = subset(sigs, p_val<0.05 & p_val_adj<0.1 & avg_log2FC >= log2(1.5) & pct.1>=0.2)
    write.table(sigs1, 'full.DEGs.txt', row.names = F)

    top50sigs <- sigs1 %>% group_by(cluster) %>% top_n(n = heatmap_n, wt = avg_log2FC)
    obj.avg <- ScaleData(obj,features = top50sigs$gene)
    obj.avg = AverageExpression(obj,return.seurat = T)
    hm = DoHeatmap(obj.avg, features = top50sigs$gene,draw.lines = F)+ theme(axis.text.y = element_blank()) + ggtitle(paste0('Heatmap of Top', heatmap_n, ' DEGs for Each Predicted Cell Type'))
    ggsave(filename = paste0("heatmap.top.DEGs.tiff"), height = heatmap.height ,width = heatmap.width,dpi = 300, compression = 'lzw')
  }

  # Figure 5H: Number of the DEGs for each cell type
  if(T){
    cat('Generating dotchart...')
    viz = data.frame(table(sigs1$cluster))
    viz = viz[order(viz$Var1, decreasing=T), ]
    viz$Var1 = factor(viz$Var1)
    dotchar = ggdotchart(viz, x = "Var1", y = "Freq",
                         color = "Var1",                                # Color by groups
                         sorting = "descending",                       # Sort value in descending order
                         add = "segments",                             # Add segments from y = 0 to dots
                         rotate = TRUE,                                # Rotate vertically
                         #group = "cyl",                                # Order by groups
                         dot.size = 3,                                 # Large dot size
                         label = round(viz$Freq,2),                        # Add mpg values as dot labels
                         font.label = list(color = "white", size = 0,
                                           vjust = 1),               # Adjust label parameters
                         ggtheme = theme_pubr()                        # ggplot2 theme
    )  + NoLegend()
    dotchar = dotchar + ylab("Number of genes") + xlab("celltype") + ggtitle('DEG Counts for Predicted Cell Type') # + geom_hline(yintercept = 50, color="blue")
    ggsave(filename = "count.full.DEGs.tiff", width=dotchart.width, height=dotchart.height, dpi=100, units="in", compression="lzw")
  }

  # Figure 5I: Functional enrichment analysis of cell type signature genes
  if (T) {
    cat('Performing enrichment analysis')
    sigs1 = read.table('full.DEGs.txt',header = T)
    top500sig <- sigs1 %>% group_by(cluster) %>% top_n(n = 500, wt = avg_log2FC) %>% data.frame()
    write.table(top500sig, file=paste0("top500.celltype.",DefaultAssay(obj),".sigs.txt"), sep="\t", col.names = T, row.names = F, quote=F)

    targets = names(table(Idents(obj)))

    dir.create(paste0(getwd(), '/enrich/'))
    kegg.list = list()
    gobp.list = list()

    for (i in 1:length(targets)){
      cell = targets[i]
      if (grepl("/", cell, fixed=TRUE)){
        cell.dir = gsub("/", ".", cell)
      } else {
        cell.dir = cell
      }
      cell.sigs <- subset(top500sig,cluster == cell)
      gostres <- gost(query = cell.sigs$gene,
                      organism = "hsapiens", ordered_query = FALSE,
                      multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                      measure_underrepresentation = FALSE, evcodes = TRUE,
                      user_threshold = 0.05, correction_method = "g_SCS",
                      domain_scope = "annotated", custom_bg = NULL,
                      numeric_ns = "", sources = NULL, as_short_link = FALSE)
      enrich <- gostres$result
      names(enrich)[1] <- "cell"
      enrich$cell <- cell
      names(enrich)[16] <- "genes"

      ##enrichment plot
      enrich.kegg <- subset(enrich,source == "KEGG", intersection_size > 3)
      if (nrow(enrich.kegg)>=10){
        enrich.kegg.top10 <- enrich.kegg[1:10,]
      } else {
        enrich.kegg.top10 <- enrich.kegg
      }

      enrich.kegg.top10$term_name <- factor(enrich.kegg.top10$term_name,levels = rev(enrich.kegg.top10$term_name))
      p1 <- ggplot(enrich.kegg.top10, aes(x = intersection_size, y = term_name)) + geom_bar(stat = "identity",aes(fill = -log10(p_value))) + scale_fill_gradientn("-log10(p-val)",colours = c("blue","red")) +theme_bw()+theme(axis.text.x = element_text(size = 12),axis.title=element_text(size=12),axis.text.y =element_text(size = 12),panel.border = element_rect(colour = "black", fill=NA, size=1)) +xlab("Count") + ylab("")
      p1 = p1 + ggtitle(paste0('KEGG Enrichment Analysis for Predicted <' , cell, '>'))
      kegg.list[[cell]] = p1

      enrich.bp <- subset(enrich,source == "GO:BP", intersection_size > 3)
      if (nrow(enrich.bp)>=10){
        enrich.bp.top10 <- enrich.bp[1:10,]
      } else {
        enrich.bp.top10 <- enrich.bp
      }

      enrich.bp.top10$term_name <- factor(enrich.bp.top10$term_name,levels = rev(enrich.bp.top10$term_name))
      p2 <- ggplot(enrich.bp.top10, aes(x = intersection_size, y = term_name)) + geom_bar(stat = "identity",aes(fill = -log10(p_value))) + scale_fill_gradientn("-log10(p-val)",colours = c("blue","red")) +theme_bw()+theme(axis.text.x = element_text(size = 12),axis.title=element_text(size=12),axis.text.y =element_text(size = 12),panel.border = element_rect(colour = "black", fill=NA, size=1)) +xlab("Count") + ylab("")
      p2 = p2 + ggtitle(paste0('GO:BP Enrichment Analysis for Predicted <' , cell, '>'))
      gobp.list[[cell]] = p2

      if (file.exists(paste0(getwd(),'/enrich/',cell.dir))){
        cell0 = gsub("/", ".", cell)
        write_xlsx(enrich,paste0(getwd(),'/enrich/',cell.dir,"/",cell0,".enrichment.result.xlsx"))
        ggsave(filename = paste0(getwd(),'/enrich/',cell.dir,"/",cell0,".enrichment.top10.kegg.tiff"),plot = p1,height = enrich.height, width = enrich.width, dpi = 300)
        ggsave(filename = paste0(getwd(),'/enrich/',cell.dir,"/",cell0,".enrichment.top10.GOBP.tiff"),plot = p2,height = enrich.height, width = enrich.width, dpi = 300)

      } else {
        dir.create(paste0(getwd(),'/enrich/',cell.dir))
        cell0 = gsub("/", ".", cell)
        write_xlsx(enrich,paste0(getwd(),'/enrich/',cell.dir,"/",cell0,".enrichment.result.xlsx"))
        ggsave(filename = paste0(getwd(),'/enrich/',cell.dir,"/",cell0,".enrichment.top10.kegg.tiff"),plot = p1,height = enrich.height, width = enrich.width, dpi = 300)
        ggsave(filename = paste0(getwd(),'/enrich/',cell.dir,"/",cell0,".enrichment.top10.GOBP.tiff"),plot = p2,height = enrich.height, width = enrich.width, dpi = 300)
      }
    }
  }


  #output PDF file
  dir.create('plots_pdf')
  pdf(file = "plots_pdf/1.umap.pdf", width = umap.width, height = umap.height)
  plot(umap_label)
  dev.off()
  pdf(file = 'plots_pdf/2.dotplot.pdf', width = dotplot.width, height = dotplot.height)
  plot(dotplot)
  dev.off()
  pdf(file = 'plots_pdf/3.violinplot.pdf', width = vlnplt.width, height = vlnplt.height)
  plot(vlnplt)
  dev.off()
  pdf(file = 'plots_pdf/4.heatmap.pdf', width = heatmap.width,height = heatmap.height)
  plot(hm)
  dev.off()
  pdf(file = 'plots_pdf/5.dotchart.pdf', width = dotchart.width, height = dotchart.height)
  plot(dotchar)
  dev.off()

  pdf(file = 'plots_pdf/6.enrich.pdf', width = enrich.width, height = enrich.height)
  for (i in 1:length(kegg.list)){
    grid.arrange(kegg.list[[i]], gobp.list[[i]] , ncol = 1)
  }
  dev.off()

  pdf_combine(input = paste0('plots_pdf/',list.files('plots_pdf')), output = 'AIO.pdf')
}
