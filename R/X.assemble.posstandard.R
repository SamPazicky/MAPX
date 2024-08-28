#' X.assemble.posstandard
#'
#' Using a database table with Entry names and GO terms, assemble a table of putatively non-interacting proteins.
#' 
#' @param standard Data frame with two columns: First for the complex name and second for comma-separated protein subunit entries.
#' @param sep Separator that separates protein IDs in the second column.
#' @import tidyverse
#' 
#' @return A data frame with pair-wise positive standard table for protein-protein interactions.
#' 
#' @examples 
#' GSpos <- X.assemble.posstandard(gold.standard, sep=";")
#' @export

X.assemble.posstandard <- function(
    standard=NULL,
    sep=";"
){
  
  if(is.null(standard)) {
    stop("Input gold standard table of interacting proteins")
  }
  # format into a normal pair-wise table
  gold_table <- standard %>% 
    select(1,2) %>% deframe() %>% as.list() %>% 
    sapply(function(x) strsplit(x, sep)) %>% 
    sapply(function(x) expand.grid.unique(x,x)) %>% 
    lapply(function(x) as.data.frame(x)) %>% bind_rows() %>%
    setNames(c("protein1","protein2")) %>%
    mutate(complex=1)
  print(paste0("Positive standard  constructed with ", nrow(gold_table), " protein pairs."))
  return(gold_table)
  
}
