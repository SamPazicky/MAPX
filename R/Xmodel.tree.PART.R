#' Xmodel.tree.PART
#'
#' PART MAP-X models.
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
#' @import RWeka
#' 
#' @return A list with two elements. $model contains the model and $predict.type contains a string that is used in
#' predict() to predict values using the model.
#' 
#' @export

Xmodel.tree.PART <- function(
    data=NULL,
    CF=0.1
) {
  
  # options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx32g"))
  library(RWeka)
  
  model <- PART(complex ~ ., data=data, control = Weka_control(C = CF))
  return(list(model=model,predict.type="probability"))
  
}

