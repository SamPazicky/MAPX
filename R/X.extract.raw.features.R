#' X.extract.raw.features
#'
#' Extract features for model training from raw data.
#' 
#' @param data Data frame with data containing a column 'protein' for protein IDs, 'length' for protein sequence length and
#' one column for abundances at each temperature labeled Ab1-AbN.
#' @param abundances Vector of integers: Which temperatures should be used to estimate protein abundances?.
#' @param remains Vector of integers: Which temperatures should be used to estimate protein remains?
#' 
#' @import tidyverse
#' 
#' @return A data frame with extracted features.
#' @examples 
#' features.raw <- X.extract.raw.features(MCdata.clean$data, abundances=c(1,2), remains=10)
#' @export

X.extract.raw.features = function(
  data=NULL,
  abundances=c(1,2),
  remains=10
) {
  
  if(is.null(data)) {
    stop("Load data")
  } else {
    data=as.data.frame(data)
    
    if(!"protein" %in% names(data) | !("length") %in% names(data) | !"Ab1" %in% names(data)) {
      stop("The loaded data must contain columns 'protein', 'length', and 'Ab1'...'AbN' where N is the number of temperatures used to generate the melting curves.")
    }
  }
  
  abcols <- paste0("Ab",abundances)
  recols <- paste0("Ab",remains)
  
  newdata <- data %>%
    mutate(Exp = rowMeans(select(., all_of(abcols)), na.rm = T)) %>%
    mutate(Remain = rowMeans(select(., all_of(recols)), na.rm = T)) %>%
    mutate(ABL=Exp/length) %>%
    mutate(Loss=Remain/Exp) %>%
    select(protein,Exp,Remain,ABL,Loss)
  
  return(newdata)
}
