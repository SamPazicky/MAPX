#' X.calculate.cors
#'
#' Using one of the network building algorithms, cluster initial network into individual protein complexes.
#' 
#' @param data Data frame with one column per temperature and protein IDs in row names.
#' @param cor.type 'type' argument in Hmisc::rcorr
#' 
#' @importFrom Hmisc rcorr
#' 
#' @return Data frame with pairwise protein interactions and R2 values
#' 
#' @examples 
#' all.predictors <- X.calculate.cors(MCdata.scaled$data)
#' @export
#' 

X.calculate.cors <- function(
    data=NULL,
    cor.type="pearson"
) {
  
  if(is.null(data)) {
    stop("Give data")
  } else {
    if(!"protein" %in% names(data)) {
      stop("Data must contain a column 'protein'.")
    }
    data <- as.data.frame(data) %>%
      column_to_rownames("protein") %>%
      mutate(across(everything(), as.numeric))
  }
  
  cormatrix <- (Hmisc::rcorr(t(as.matrix(data)),type=cor.type))
  
  cordata <- data.frame(row=rownames(cormatrix[["r"]])[row(cormatrix[["r"]])[upper.tri(cormatrix[["r"]])]],
                        col=colnames(cormatrix[["r"]])[col(cormatrix[["r"]])[upper.tri(cormatrix[["r"]])]],
                        corr=cormatrix[["r"]][upper.tri(cormatrix[["r"]])])
  
  names(cordata) <- c("protein1","protein2","R2")
  
  return(cordata)
  
}
