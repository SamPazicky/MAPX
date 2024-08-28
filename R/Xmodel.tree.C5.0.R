#' Xmodel.tree.C5.0
#'
#' C5.0 MAP-X models.
#' 
#' @param data Data frame with pairwise differential values. Must contain columns "protein1" and "protein2", any number of columns with
#' predictors for modelling and a labels column with values 1 (for protein pairs that form a complex) and 0 (for protein pairs that
#' do not form a complex).
#' @param eval.metric Character string: How should the model be evaluated in cross-validation? Defalt is "prc" for area under the precision-recall curve.
#' Other options are "roc" for area under the receiver-operator curve and "kappa" for Cohen's kappa.
#' @param CF Numeric vector: Complexity parameters to be cross-validated. Default is 0.1.
#' @param costs Integer: A number specifying how much is the cost of falsely predicting non-interacting protein pairs as interacting higher than vice versa.
#' @param downsample Integer: How many times less of the non-interacting proteins should be used for the training? Default is 1. Applicable for C5.0 models.
#' 
#' @import tidyverse
#' @import C50
#' 
#' @return A list with two elements. $model contains the model and $predict.type contains a string that is used in
#' predict() to predict values using the model.
#' 
#' @export

Xmodel.tree.C5.0 <- function(
    data=NULL,
    costs=NA,
    winnowing=FALSE,
    noGlobalPruning=FALSE,
    CF=0.3,
    boost=1,
    downsample=NA
) {
  
  require(C50)
  
  cost_mat <- Xconvert_costs(costs)
  
  data <- Xdownsample(data,downsample)

  model <- C5.0(complex ~ ., data=data, na.action=na.omit,costs=cost_mat, trials=boost,
                control=C5.0Control(winnow=winnowing,CF=CF,noGlobalPruning=noGlobalPruning))
 
  return(list(model=model,predict.type="prob"))
  
}

