#devtools::install_github("immunogenomics/lisi")
library(lisi)
library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(rstatix)
library(ggpubr)

# load CellRef and HLCA objects
cellref = readRDS('.Data/LungMAP_HumanLung_CellRef.v1.rds')
hlca = readRDS('.Data/HLCA.rds')


# cLISI -----------------------------------------------------------------

  # calculate cLISI for CellRef cells
  coords = cellref@reductions$umap@cell.embeddings
  meta_data = cellref@meta.data[rownames(coords), c('DonorID','celltype_level3')]
  res <- compute_lisi(coords, meta_data, c("celltype_level3"))
  df = data.frame(row.names = rownames(res), cLISI = res$celltype_level3)
  colnames(df) = paste0('CellRef_', colnames(df))
  df$CellRef_Celltype = cellref@meta.data[rownames(df), 'celltype_level3']
  write.table(df,'clisi.cellref.txt', sep='\t')

  # calculate cLISI for HLCA cells
  coords = hlca@reductions$umap@cell.embeddings
  meta_data = hlca@meta.data[rownames(coords), c('subject_ID','ann_finest_level')]
  res <- compute_lisi(coords, meta_data, c("ann_finest_level"))
  df = data.frame(row.names = rownames(res), cLISI = res$ann_finest_level)
  colnames(df) = paste0('HLCA_', colnames(df))
  df$HLCA_Celltype = hlca@meta.data[rownames(df), 'ann_finest_level']
  write.table(df, 'clisi.hlca.txt',sep='\t')


# Celltype level iLISI ----------------------------------------------------------

  # Calculate iLISI for CellRef cells
  clusters = table(cellref$celltype_level3)
  clusters = names(clusters)
  myilisi = data.frame()
  mycount = data.frame()
  for (i in 1:length(clusters)) {
    tmp = subset(cellref, celltype_level3 == clusters[i])
    coords = tmp@reductions$umap@cell.embeddings
    meta_data = tmp@meta.data[rownames(coords), c("DonorID",'DataID','celltype_level3')]
    perplextiy.donor = length(table(tmp$DonorID))
    perplexity.data = length(table(tmp$DataID))

    df.count = data.frame(cluster= clusters[i], DonorID = length(table(tmp$DonorID)), DataID = length(table(tmp$DataID)),cellCount= dim(tmp@meta.data)[1])
    mycount = rbind(mycount, df.count)
    write.table(mycount, './cellref.count.txt', row.names = F, sep = '\t')

    # 3 * unique number of batches < cell number, perplexity remains same
    if (3* df.count$DonorID <= df.count$cellCount  ) {
      perplexity.donor = length(table(tmp$DonorID))
      df.count$donor.Normal = T
    }
    if (3* df.count$DataID <= df.count$cellCount) {
      perplexity.data = length(table(tmp$DataID))
      df.count$data.Normal = T
    }

    # 3 * unique number of batches > cell number, perplexity = cell number / 3
    if (3 * df.count$DonorID > df.count$cellCount) {
      perplexity.donor = floor(df.count$cellCount / 3)
      df.count$donor.Normal = F
    }

    if (3 * df.count$DataID > df.count$cellCount) {
      perplexity.data= floor(df.count$cellCount / 3)
      df.count$data.Normal = F
    }

    res1 <- compute_lisi(coords, meta_data, c("DonorID"),perplexity = perplexity.donor)
    res2 <- compute_lisi(coords, meta_data, c("DataID"),perplexity = perplexity.data)


    df = data.frame(row.names = rownames(res1), DonorID = res1, DataID = res2)
    colnames(df) = paste0('CellRef_', colnames(df))
    df$CellRef_Celltype = tmp@meta.data[rownames(df), 'celltype_level3']
    df$CellRef_DonorID.norm = df$CellRef_DonorID / length(table(tmp$DonorID))*100
    df$CellRef_DataID.norm = df$CellRef_DataID / length(table(tmp$DataID))*100
    myilisi = rbind(myilisi, df)


  }
  write.table(myilisi, './ilisi.cellref.txt', row.names = T, sep = '\t')


  # Calculate iLISI for HLCA cells
    clusters = table(hlca$ann_finest_level)
    clusters = names(clusters)

    myilisi = data.frame()
    mycount = data.frame()
    for (i in 1:length(clusters)) {
      tmp = subset(hlca, ann_finest_level == clusters[i])
      coords = tmp@reductions$umap@cell.embeddings
      meta_data = tmp@meta.data[rownames(coords), c("subject_ID",'sample','ann_finest_level')]
      perplexity.donor= length(table(tmp$subject_ID))
      perplexity.data = length(table(tmp$sample))

      df.count = data.frame(cluster= clusters[i], DonorID = length(table(tmp$subject_ID)), DataID = length(table(tmp$sample)), cellCount = dim(tmp@meta.data)[1])
      mycount = rbind(mycount, df.count)
      write.table(mycount, './hlca.count.txt', row.names = F, sep = '\t')


      if (3* df.count$DonorID <= df.count$cellCount  ) {
        perplexity.donor = length(table(tmp$subject_ID))
        df.count$donor.Normal = T
      }
      if (3* df.count$DataID <= df.count$cellCount) {
        perplexity.data = length(table(tmp$sample))
        df.count$data.Normal = T
      }

      if (3 * df.count$DonorID > df.count$cellCount) {
        perplexity.donor = floor(df.count$cellCount / 3)
        df.count$donor.Normal = F
      }
      if (3 * df.count$DataID > df.count$cellCount) {
        perplexity.data= floor(df.count$cellCount / 3)
        df.count$data.Normal = F
      }

      res1 <- compute_lisi(coords, meta_data, c("subject_ID"),perplexity = perplexity.donor)
      res2 <- compute_lisi(coords, meta_data, c("sample"), perplexity = perplexity.data)

      df = data.frame(row.names = rownames(res1), DonorID = res1, DataID = res2)
      colnames(df) = paste0('HLCA_', colnames(df))
      df$HLCA_Celltype = tmp@meta.data[rownames(df), 'ann_finest_level']
      df$HLCA_DonorID.norm = df$HLCA_subject_ID / length(table(tmp$subject_ID))*100
      df$HLCA_DataID.norm = df$HLCA_sample / length(table(tmp$sample))*100
      myilisi = rbind(myilisi, df)

    }
    write.table(myilisi, './ilisi.hlca.txt', row.names = T, sep = '\t')



#Session Info
sink('sessionInfo.txt')
sessionInfo()
sink()
