#' X.average.crossmodels
#'
#' Average cross-trained models to get final prediction values for one dataset.
#' 
#' @param data Data frame with columns for protein1, protein2 and columns 'score', 'cycle' and 'split'.
#' @param scores.col Character string: Name of the columns with prediction values.
#' @param labels.col Character string: Name of the columns with labels.
#' @param evaluate Logical: Should the averaged predictions be evaluated? Default is FALSE. If TRUE, standard.set must be included
#' @param standard.set Dataframe with columns protein1, protein2 and another column with labels.
#' @param eval.metric Character string: Evaluation metric to be used for evaluating the model performance. Options same as in X.evaluate.
#' 
#' @import tidyverse
#'
#' @return A list with four elements: $data is the data.frame with averaged predictions, $eval.metric is the integer specifying the
#' resulting evaluation metric, $metric.plotdata are underlying data for plotting the evaluation metric and $plot is the plot itself.
#' @examples 
#' comp.models[[rep]] <- X.average.crossmodels(data=cross.model$data, eval.metric="prc", evaluate=TRUE, standard.set=GS_specific)
#' @export
#' 
X.average.crossmodels <- function(
  data=NULL,
  scores.col=NA,
  labels.col=NA,
  evaluate=FALSE,
  standard.set=NULL,
  eval.metric="prc"
) {
  

  if(is.null(data)) {
    stop("Please include data.")
  } else {
    data <- as.data.frame(data)
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
    scores.col <- names(data)[ncol(data)]
    cat("No scores column selected. Selecting the last column:",scores.col,"\n")
  }
  
  
  cat("Assessing data...\n")
  
  unique.data <- data %>% dplyr::select(starts_with("protein")) %>% distinct(.keep_all=TRUE)
  all.proteins <- unique.data %>% unlist() %>% unname() %>% unique()
  all.proteins.table <- expand.grid.unique(all.proteins,all.proteins) %>%
    as.data.frame() %>% setNames(paste0("protein",1:2)) %>% mutate(here=1)

  cat("Proteins in the dataset:",length(all.proteins),"\n",
      "Number of possible interactions:", nrow(all.proteins.table), "\n",
      "Number of predicted interactions:", nrow(unique.data),"\n")
  if(nrow(all.proteins.table)>nrow(unique.data)) {
    message("The missing predicted interactions are likely the interactions included in the golden standard that have by random never been
            selected during cross-training because of rounding down during splitting. Increase the number of cycles during cross-training.")
  }
  cat("Averaging scores...\n")
  av.data <- data %>%
    rename(subscore=score) %>%
    group_by(protein1,protein2) %>%
    summarise(score=mean(subscore,na.rm=TRUE)) %>%
    ungroup()
  
 if(evaluate) {
    cat("Evaluating...\n")
    to.evaluate <- av.data %>%
      cross_join(standard.set, vars=c("protein1","protein2"), mode="inner")
    
    eval <- X.evaluate(to.evaluate, scores.col="score",labels.col="complex", eval.metric=eval.metric, plot=TRUE)
    
    output <- list(
      data=av.data,
      eval.metric=eval$eval.metric,
      plot=eval$plot
    )
  } else {
    output <- list(data=av.data,
                   eval.metric=NA,
                   plot=NA
                   )
  }
  
  return(output)
}
