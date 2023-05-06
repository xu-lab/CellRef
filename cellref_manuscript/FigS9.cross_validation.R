library(rBayesianOptimization)
library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(mltools)

set.seed(1234)

cellref = readRDS(".Data/LungMAP_HumanLung_CellRef.v1.rds")
cellref = subset(cellref, celltype_level3 %in% filter(data.frame(table(cellref$celltype_level3)), Freq>500)$Var1)

transfer_labels = function(ref_obj, query_obj, md){
  transfer_anchors = FindTransferAnchors(
    reference = ref_obj,
    query = query_obj,
    reference.reduction = "pca",
    normalization.method = "SCT",
    dims = 1:50
  )

  predictions = TransferData(
    anchorset = transfer_anchors,
    reference = cellref_ref,
    refdata = ref_obj@meta.data[,md],
    dims = 1:50
  )
  query_obj$predicted.id = predictions$predicted.id
  query_obj
}

# Function to calculate metrics evaluating the performance of the classifier.
# This function was modified from evaluate.R in the https://github.com/tabdelaal/scRNAseq_Benchmark (PMID: 31500660)
# We added the calculation of Matthews Correlation Coefficient (MCC)
evaluate <- function(true_lab, pred_lab, Indices = NULL){

  if (! is.null(Indices)){
    true_lab <- true_lab[Indices]
    pred_lab <- pred_lab[Indices]
  }

  unique_true <- unlist(unique(true_lab))
  unique_pred <- unlist(unique(pred_lab))

  unique_all <- unique(c(unique_true,unique_pred))
  conf <- table(true_lab,pred_lab)
  pop_size <- rowSums(conf)

  pred_lab = gsub('Node..','Node',pred_lab)

  conf_F1 <- table(true_lab,pred_lab,exclude = c('unassigned','Unassigned','Unknown','rand','Node','ambiguous','unknown'))

  F1 <- vector()
  Matt <- vector()
  sum_acc <- 0

  for (i in c(1:length(unique_true))){
    findLabel = colnames(conf_F1) == row.names(conf_F1)[i]
    if(sum(findLabel)){
      prec <- conf_F1[i,findLabel] / colSums(conf_F1)[findLabel]
      rec <- conf_F1[i,findLabel] / rowSums(conf_F1)[i]
      Matt[i] = mcc(TP=conf_F1[i,findLabel],
                    FN=rowSums(conf_F1)[i]-conf_F1[i,findLabel],
                    FP=colSums(conf_F1)[findLabel]-conf_F1[i,findLabel],
                    TN=sum(conf_F1)-colSums(conf_F1)[findLabel]-rowSums(conf_F1)[i]+conf_F1[i,findLabel])
      if (prec == 0 || rec == 0){
        F1[i] = 0
      } else{
        F1[i] <- (2*prec*rec) / (prec + rec)
      }
      sum_acc <- sum_acc + conf_F1[i,findLabel]
    } else {
      Matt[i] = mcc(TP=0,
                    FN=rowSums(conf_F1)[i],
                    FP=0,
                    TN=sum(conf_F1)-rowSums(conf_F1)[i])
      F1[i] = 0
    }
  }

  pop_size <- pop_size[pop_size > 0]

  names(F1) <- names(pop_size)
  names(Matt) <- names(pop_size)

  med_F1 <- median(F1)

  total <- length(pred_lab)
  num_unlab <- sum(pred_lab == 'unassigned') + sum(pred_lab == 'Unassigned') + sum(pred_lab == 'rand') + sum(pred_lab == 'Unknown') + sum(pred_lab == 'unknown') + sum(pred_lab == 'Node') + sum(pred_lab == 'ambiguous')
  per_unlab <- num_unlab / total

  acc <- sum_acc/sum(conf_F1)

  result <- list(Conf = conf,  # confusion matrix
                 MedF1 = med_F1, # median F1 score
                 F1 = F1, # F1 score per class
                 Acc = acc, # accuracy
                 MCC=Matt, # Matthews correlation coefficient
                 PercUnl = per_unlab,  # percentage of unlabeled cells
                 PopSize = pop_size #  number of cells per cell type
                 )

  return(result)
}

results_df = data.frame(matrix(NA, nrow = length(unique(cellref$celltype_level3)), ncol = 0))

rownames(results_df) = unique(cellref$celltype_level3)

predictions_df = data.frame(matrix(NA, nrow = 0, ncol = 3))

colnames(predictions_df) = c("Run", "Predicted", "Original")

folds <- KFold(cellref$celltype_level3,nfolds = 10, stratified = TRUE)

for(i in 1:length(folds)){
  cellref_test = cellref[, c(folds[[i]])]
  cellref_test = transfer_labels(cellref[, -c(folds[[i]])], cellref_test, "celltype_level3")

  predictions_df = rbind(predictions_df, data.frame(Run=paste0("Run",i),Prediction=cellref_test$predicted.id,
                                                    Original=cellref_test$celltype_level3))
  evaluation = evaluate(cellref_test$celltype_level3, cellref_test$predicted.id)
  results_df[names(evaluation[["F1"]]),paste0("F1",i)]= evaluation[["F1"]]
  results_df[names(evaluation[["MCC"]]),paste0("MCC",i)]= evaluation[["MCC"]]
}

write.csv(results_df,"cross_validation.csv")
write.csv(predictions_df,"cross_validation_prediction.csv")
