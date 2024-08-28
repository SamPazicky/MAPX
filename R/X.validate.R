#' X.validate
#'
#' Validate models.
#' 
#' @param data Data frame with measured parameters.
#' @param model The model.
#' @param predict.type The prediction type of the model (typically 'prob' or 'probability')
#' @param labels.col Character string: 'data' column name with labels (1 for complex-forming, 0 for others.)
#' @param eval.metric Character string: Method for predictor evaluation. 'roc' for area under the receiver-operator curve, 'prc' for area under the
#' precision-recall curve, 'kappa' for Cohen's kappa and 'F1' for F1 score.
#' @param plot Logical: Should the plots be plotted and outputted?
#' @param black Logical: Should the plotted curve be black or color gradient? Default is FALSE.
#' @param kappa.cutoff Numeric: What is the score cutoff for evaluating based on Cohen's kappa? Default is 0.5.
#' @param kappa.weight Character string: Same as argument weight in irr:kappa2
#' @param F1.cutoff Numeric: What is the score cutoff for evaluating based on F1? Default is 0.5.
#' @param evaluate Logical: Should the model be evaluated?
#' 
#' @import tidyverse
#' @import pROC
#' @import PRROC
#' @import irr
#' @import viridis
#' 
#' @return A list with three elements. $eval.metric is a number corresponding to the chosen evaluation method. $curvedata are data that allow plotting
#' a plot corresponding to the ealuation metric. $plot is the plot of the corresponding metric.
#' @examples 
#' evaluation <- X.validate(data,model,"score","complex","prc")
#' @export
#' 
X.validate <- function(
    data=NULL, # data frame with data
    model=NULL,
    predict.type="prob",
    labels.col=NA,
    eval.metric="prc", # roc, prc, kappa, F1 F1 is still not implemented
    plot=FALSE,
    black=FALSE,
    kappa.cutoff=0.5,
    kappa.weight="unweighted",
    F1.cutoff=0.5,
    evaluate=TRUE
) 
  
{
  predictions <- predict(model,data,type=predict.type) %>% as.data.frame() %>% pull("1")
  updata <- data %>% mutate(score=predictions)
  if(evaluate) {
    eval <- X.evaluate(updata, scores.col="score",labels.col="complex", eval.metric=eval.metric,
                       plot=plot,black=black,kappa.cutoff=kappa.cutoff,kappa.weight=kappa.weight,F1.cutoff=F1.cutoff)
  } else {
    eval <- list(data=updata)
  }
  return(eval)
  
}
