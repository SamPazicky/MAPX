#' X.assemble.traindata 
#'
#' Using available data and the gold standard, assemble training set for the model training.
#' 
#' @param data Data frame with columns 'protein1', 'protein2' and other predictor columns.
#' @param standard Data frame with columns 'protein1', 'protein2' and 'complex'. The column complex should have values of 1 or 0
#' for interacting and non-interacting proteins, respectively.
#' 
#' @import tidyverse
#' 
#' @return A data frame with complex labels from standard assigned to observed data.
#' @examples 
#' GS_specific <- X.assemble.traindata(all.predictors[[tp]][[rep]],GS)
#' @export

X.assemble.traindata <- function(
    data=NULL,
    standard=NULL
){

  if(is.null(data)|is.null(standard)) {
    stop("Input data and standard.")
  } else if(!"protein1" %in% names(data) | !"protein2" %in% names(data)) {
    stop("The data must contain columns protein1 and protein2")
  } else if(!"protein1" %in% names(standard) | !"protein2" %in% names(standard) | !"complex" %in% names(standard)) {
    stop("The standard must contain columns protein1, protein2 and complex.")
  } else if(paste(sort(unique(standard$complex)),collapse=";;")!="0;;1") {
    stop("The column complex in the standard can only contain values 0 and 1 for non-interacting and interacting proteins, respectively.")
  } else {
    data <- data %>% as.data.frame()
    standard <- standard %>% as.data.frame() %>% mutate(complex=as.numeric(complex))
  }

  output <- MAPX::cross_join(data,standard,vars=c("protein1","protein2"), mode="inner") %>%
    select(protein1,protein2,everything(),complex) %>%
    arrange(protein1,protein2)
  
  return(output)
  
}
