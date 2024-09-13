#' X.average.reps
#'
#' Average predictions for multiple models across replicates.
#' 
#' @param data List of data frame with columns for protein1, protein2 and a column with predictions.
#' @param scores.col Character string: Name of the column with prediction values.
#' @param weights Numeric vector: Typically values between 0-1, giving weights for each replicate for weighted averaging.
#' @param evaluate Logical: Should the averaged predictions be evaluated? Default is FALSE. If TRUE, standard.set must be included
#' @param standard.set Dataframe with columns protein1, protein2 and another column with labels.
#' @param labels.col Character string: Name of the standard.set column with labels.
#' @param eval.metric Character string: Evaluation metric to be used for evaluating the model performance. Options same as in X.evaluate.
#' 
#' @import tidyverse
#' 
#' @return A list with four elements: $data is the data.frame with averaged predictions, $eval.metric is the integer specifying the
#' resulting evaluation metric, $eval.data are underlying data for plotting the evaluation metric and $plot is the plot itself.
#' @examples 
#' final.model <- X.average.reps(data=comp.models.data, evaluate=TRUE,standard.set=GS)
#' @export
#' 
X.average.reps <- function(
    data=NULL,
    scores.col=NA,
    weights=NULL,
    evaluate=FALSE,
    standard.set=NULL,
    labels.col=NA,
    eval.metric="prc"
) {
  
  if(is.null(data)) {
    stop("Please include data")
  } else if(length(data)==1) {
    
    if(evaluate) {
      message("There is only one replicate in predictions, no need to average.")
      cat("Evaluating...\n")
      to.evaluate <- data[[1]] %>%
        cross_join(standard.set, vars=c("protein1","protein2"), mode="inner")
      
      eval <- X.evaluate(to.evaluate, scores.col=scores.col,labels.col=labels.col, eval.metric=eval.metric, plot=TRUE)
      
      output <- list(
        data=data,
        eval.metric=eval$eval.metric,
        eval.data=eval$curvedata,
        plot=eval$plot
      )
    } else {
      output <- list(data=av.data,
                     eval.metric=NA,
                     eval.data=NA,
                     plot=NA
      )
    }
    return(output)
    
  } else {
    predictions <- lapply(data, function(x) as.data.frame(x))
  }
  
  if(is.null(names(data))) {
    names(data) = paste0("rep",1:(length(data)))
  }
  
  if(evaluate) {
    if(is.null(standard.set)) {
      stop("If evaluate is TRUE, standard.set must be included.")
    } else {
      standard.set <- as.data.frame(standard.set)
    }
    if(is.na(labels.col)) {
      labels.col <- names(standard.set)[ncol(standard.set)]
      cat("No labels column selected. Selecting the last column:",labels.col,"\n")
    }
  }
  if(is.na(scores.col)) {
    scores.col <- names(data[[1]])[ncol(data[[1]])]
    cat("No scores column selected. Selecting the last column:",scores.col,"\n")
  }
  
  if(!is.null(weights)) {
    if(length(weights)!=length(data)) {
      stop("Weights must have the same length as data.")
    }
    w.frame <- data.frame(weight=weights) %>%
      mutate(rep=names(data))
  } else {
    w.frame <- data.frame(weight=rep(1,length(data))) %>%
      mutate(rep=names(data))
  }
  
  cat("Averaging...\n")
  av.data <- list_rbind(data, names_to="rep") %>%
    left_join(w.frame,by="rep") %>%
    rename(rep_score=!!sym(scores.col)) %>%
    group_by(protein1,protein2) %>%
    dplyr::summarise(score=weighted.mean(rep_score,weight), .groups="keep") %>%
    ungroup()
 
  if(evaluate) {
    cat("Evaluating...\n")
    to.evaluate <- av.data %>%
      cross_join(standard.set, vars=c("protein1","protein2"), mode="inner")
    
    eval <- X.evaluate(to.evaluate, scores.col="score",labels.col=labels.col, eval.metric=eval.metric, plot=TRUE)
    
    output <- list(
      data=av.data,
      eval.metric=eval$eval.metric,
      eval.data=eval$curvedata,
      plot=eval$plot
    )
  } else {
    output <- list(data=av.data,
                   eval.metric=NA,
                   eval.data=NA,
                   plot=NA
    )
  }
  return(output)
}
