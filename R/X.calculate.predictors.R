#' X.calculate.predictors
#'
#' Using one of the network building algorithms, cluster initial network into individual protein complexes.
#' 
#' @param data Data frame with a column 'protein' and other columns with features.
#' @param features Vector of character strings: Names of columns with features. Default is NULL, in which case all but the protein column will be taken.
#' @param funs  Vector of character strings: Functions to be used for calculation of predictors, matching the order of predictors given in vars.
#' @param prefixes  Vector of character strings: Prefixes for naming of new columns in the new table. One for each var.
#' 
#' @import tidyverse
#' 
#' @return Data frame with all predictors for pairwise protein interactions.
#' @examples 
#' feature.predictors <- X.calculate.predictors(data=features, 
#'       features=c("Ti","logAUC","logABL","logLoss","Penalty_trans"), 
#'       funs=c(rep("absdif",4),"Xsum"),
#'       prefixes=c(rep("d",4),"sum")
#'     )
#' @export
#' 


X.calculate.predictors <- function(
    data=NULL,
    features=NULL,
    funs=c(),
    prefixes=NA
) {
  
  if(is.null(data)) {
    stop("Include data")
  } else {
    data <- as.data.frame(data)
    if(!"protein" %in% names(data)) {
      stop("The input data must contain a column 'protein'.")
    }
    for(f in features) {
      if(!f %in% names(data)) {
        stop(paste0("The input data must contain a column", f,"."))
      }
    }
  }
  
  # Custom function for expression differences
  expdif <- function(x, y) {            
    z <- log2(abs(x - y))
    return(z)
  }
  absdif <- function(x, y) {
    z <- abs(x-y)
    return(z)
  }
  Xproduct <- function(x, y) {
    z <- x*y
    return(z)
  }
  Xsum <- function(x, y) {
    z <- x+y
    return(z)
  }
  Xmean <- function(x,y) {
    z <- (x+y)/2
    return(z)
  }
  Xlower <- function(x,y) {
    z <- min(c(x,y))
    return(z)
  }
  
  add_data <- list()
  for(i in seq_along(features)) {
    var=features[i]
    fun=get(funs[i])
    if(length(prefixes)==1) {
      if(!is.na(prefixes)) {
        prefix=prefixes[i]
      } else {
        prefix=""
      }
    } else {
      prefix=prefixes[i]
    }
    
    add_matrix <- data.frame(outer(data%>%pull(!!sym(var)),data%>%pull(!!sym(var)),FUN = fun))
    colnames(add_matrix) <- data%>%pull(protein)
    rownames(add_matrix) <- data%>%pull(protein)
    add_data[[i]] <- data.frame(row=rownames(add_matrix)[row(add_matrix)[upper.tri(add_matrix)]],
                           col=colnames(add_matrix)[col(add_matrix)[upper.tri(add_matrix)]],
                           addvar=add_matrix[upper.tri(add_matrix)])
    
    names(add_data[[i]]) <- c("protein1","protein2",paste0(prefix,var))
  }
  
  output <- Reduce(function(x,y) cross_join(x,y,vars=paste0("protein",c(1,2)), mode="full"), 
                   add_data)
  
  return(output)
  
}
