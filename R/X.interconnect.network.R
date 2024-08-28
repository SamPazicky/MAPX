#' X.interconnect.network
#'
#' With pairwise interaction table as an input, interconnect protein complexes such that each two subunits in the same cluster whose pairwise
#' interaction does not pass the cutoff is considered as an interaction.
#' 
#' @param data A data.frame with pair-wise interaction table.
#' @param scores.col Character string: Name of the columns with prediction values.
#' @return A data.frame with pair-wise interaction table.
#' @examples 
#' pairwise.table <- X.interconnect.network(averaged_predictions$predicted.dataset %>% filter(pred>=sug_cutoff))
#' @export
#' 

X.interconnect.network <- function(
    data=NULL,
    scores.col=NA 
)
{
  
  if(is.null(data)) {
    stop("Please specify data.")
  } else {
    if(!"protein1" %in% names(data) | !"protein2" %in% names(data)) {
      stop("The data must contain columns protein1, protein2 and a column with model scores.")
    }
    data <- as.data.frame(data)
    if(is.na(scores.col)) {
      scores.col <- names(data)[ncol(data)]
    }
    data <- data %>% rename(score=!!sym(scores.col)) %>% dplyr::select("protein1","protein2","score")
  }
  
  complexes <- X.pairwise.to.complexes(data)
  newdata <- X.complexes.to.pairwise(complexes)
  
  return(newdata)
  
}
