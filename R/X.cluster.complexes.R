#' X.cluster.complexes
#'
#' Internal function with simple clustering algorithm.
#' 
#' @param data Data frame with pair-wise interactions (columns 'protein1' and 'protein2') and a column with probabilities from ML model.
#' @param scores.col Character string: Name of the columns with prediction values.
#' 
#' @import tidyverse
#' 
#' @examples 
#' clustered <- X.cluster.complexes(data=cut_predictions, scores.col="score")

X.cluster.complexes <- function(
    data=NULL,
    scores.col=NA
) {
  
  if(is.null(data)) {
    stop("Please specify data.")
  } else {
    if(is.na(scores.col)) {
      scores.col <- names(data)[ncol(data)]
    }
    data <- as.data.frame(data) %>%
      rename(score=!!sym(scores.col))
  }
  
  # initiate the building
  grouped_data <- data %>% 
    select(starts_with("protein"),score) %>%
    arrange(desc(score))
  complexes=list()
  building="inprogress"
  
  while(building=="inprogress") {
    FMs <- grouped_data %>% slice_head(n=1) %>% select(-score) %>% unlist() %>% as.vector()
    GPs <- grouped_data %>% filter(protein1 %in% FMs | protein2 %in% FMs) %>% select(-score) %>%
      t() %>% c() %>% unique() %>% .[!. %in% FMs]
    if(length(GPs)==0) {
      complexes[[paste(FMs,collapse=";")]] <- FMs
      grouped_data <- grouped_data %>% slice(-1)
      if(nrow(grouped_data)==0) {
        building="Done"
      }
      next
    }

    complexes[[paste(FMs,collapse=";")]] <- c(FMs,GPs)
    grouped_data <- grouped_data %>%
      filter(!protein1 %in% c(FMs,GPs)) %>%
      filter(!protein2 %in% c(FMs,GPs))


    if(nrow(grouped_data)==0) {
      building="Done"
    }
  }

  # retrieve back the protein couples with prediction values
  grouped_data <- data %>% select(starts_with("protein"),score) %>% arrange(desc(score))


  complexes_table <- X.complexes.to.pairwise(complexes) %>%
    cross_join(grouped_data,vars=c("protein1","protein2"),mode="inner")

  return(complexes_table)
  
}
