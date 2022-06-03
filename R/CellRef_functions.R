#' Perform batch correction of single cell data
#'
#' This function loads a Seurat object that contains all data for batch correction and integration.
#'  It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param object A Seurat object
#' @param integration.batch a character string indicating the meta.data column that contains the batch for correction
#' @param ident.var a character string indicating the meta.data column that contains the identities of each cell
#' @param method the method used for batch correction
#' @param npcs the number of principal components used for batch correction and umap projection
#' @param vars.to.regression the covariates
#' @param umap.min_dist
#' @param umap.n_neighbors
#' @param do.clustering
#' @param clustering.method
#' @param clustering.resolution
#' @param rpca.normalization
#' @param rpca.k.anchor
#' @param rpca.k.weight
#' @param rpca.prune_anchors
#' @param verbose
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

  DefaultAssay(obj) = "RNA"
  object = NormalizeData(object)
  object = FindVariableFeatures(obj)
  object = ScaleData(object, vars.to.regress=vars.to.regress)

  object = RunPCA(obj, npcs=npcs)

  hmat = HarmonyMatrix(object@reductions$pca@cell.embeddings[, 1:npcs],
                       meta_data = object@meta.data, vars_use = integration.batch, do_pca = FALSE)
  object@reductions$harmony = CreateDimReducObject(embeddings = hmat, assay= DefaultAssay(object), key="harmony_")
  object = RunUMAP(object, dims=1:npcs, reduction = "umap",
                   reduction.name = "harmony", return.model = TRUE,
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
#' @param ctd a data.frame object that contains the cell type dictionary
#' @param cluster.var a character string indicating the meta.data column that contains cluster assignment of each cell
#' @param exp.thresh
#' @param use.scale.exp
#' @param score.thresh
#' @param seed.n.min
#' @param seed.n.max
#' @param group.size.min
#' @param precision.thresh
#' @param recall.thresh
#' @param verbose
#' @return a data.frame containing the mapping of candidate cell clusters and cell types
#' @export
findCandidateClusters <- function(object, ctd, cluster.var="integrated_clusters",
                                  exp.thresh = 0, use.scaled.exp=T,
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
        if (use.scaled.exp==T) {
          i.x = as.matrix(object@assays$RNA@scale.data[i.markers[i.markers.p, "Marker"], i.candidates])
        } else {
          i.x = as.matrix(object@assays$RNA@data[i.markers[i.markers.p, "Marker"], i.candidates])
        }

        i.x = i.x[, which(colSums(i.x>exp.thresh)>1)]

        i.cells = c()
        i.list <- list()
        for (j in 1:dim(i.x)[1]) {
          i.j.cells = NULL
          i.j.cells <- i.x[j,]
          i.j.cells <- i.j.cells[i.j.cells > exp.thresh]
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
#' @param ctd a data.frame object that contains the cell type dictionary for seed cell identification.
#' @param candidate.clusters a data.frame that contains the mapping between candidate cell clusters and cell types in the dictionary
#' @param cluster.var a character string indicating the meta.data column that contains cluster assignment of each cell
#' @param exp.thresh the number of principal components used for batch correction and umap projection
#' @param use.scaled.exp the covariates
#' @param score.thresh
#' @param seed.n.min
#' @param seed.n.max
#' @param verbose
#' @return A Seurat object of the identified seed cells
#' @export
# ctd: CellType, Marker, MarkerType
# candidate_clusters: Cluster, Type
findSeedCells <- function(object, ctd, candidate.clusters,
                          cluster.var="integrated_clusters",
                          exp.thresh = 0, use.scaled.exp=T,  score.thresh=Inf,
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
        if (use.scaled.exp==T) {
          i.x = as.matrix(object@assays$RNA@scale.data[i.markers[i.markers.p, "Marker"], i.candidates])
        } else {
          i.x = as.matrix(object@assays$RNA@data[i.markers[i.markers.p, "Marker"], i.candidates])
        }

        # cells with at least two markers
        i.x = i.x[, which(colSums(i.x>exp.thresh)>1)]

        i.cells = c()
        i.list <- list()
        for (j in 1:dim(i.x)[1]) {
          i.j.cells = NULL
          i.j.cells <- i.x[j,]
          i.j.cells <- i.j.cells[i.j.cells > exp.thresh]
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

#' Build a CellRef using the initial data integration and the seed cells
#'
#' @param object.all the Seurat object that contains all data for CellRef construction after batch correction
#' @param object.seed the Seurat object that contains the identified seed cells.
#' @param mapping.batch a character string indicating the meta.data column that contains data batch information to perform mapping to the seed cells
#' @param integration.npcs the number of principal components
#' @param rpca.prune_anchors
#' @param purity.pruning
#' @param purity.nn.k
#' @param purity.thresh
#' @param verbose
#' @return A Seurat object
#' @export
buildCellRef <- function(object.all, object.seed,
                         mapping.batch="Dataset",
                         integration.batch="Dataset", integration.npcs=200,
                         rpca.prune_anchors=T,
                         purity.pruning=T, purity.nn.k=20, purity.thresh=0.6,
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

  cellref_cells = obj_seed@meta.data
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
                             npcs=integration.npcs, rpca.prune_anchors = T,
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
                               verbose=F)

  }

  return(object)
}
