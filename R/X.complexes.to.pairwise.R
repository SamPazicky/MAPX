#' X.complexes.to.pairwise
#'
#' Transform a list of protein complexes into a pairwise interaction table.
#' 
#' @param complexes A list of vectors where each vector represents subunits of a protein complex.
#' @import tidyverse
#' @return A data.frame with pair-wise interaction table.
#' @examples 
#' pairwise.table <- X.complexes.to.pairwise(list_of_complexes)
#' @export
#' 

X.complexes.to.pairwise <- function(
    complexes=list()
){
  
  if(length(complexes)==0) {
    stop("Please include a list of protein complexes")
  }
  
  if((sapply(names(complexes),function(x) (x)=="") %>% any())) {
    names(complexes)=as.character(1:length(complexes))
  }
  
  if(is.null(names(complexes))) {
    names(complexes)=as.character(1:length(complexes))
  }
  
  # rebuild table from complexes
  
  rebuilt_list <- list()
  for (complex in names(complexes)) {
    if(length(complexes[[complex]])<2) {
      next
    }
    all_combs <- combn(gtools::mixedsort(complexes[[complex]]), 2) %>% t() %>% as.data.frame() %>% setNames(c("protein1","protein2"))
      # cross_join(data,vars=c("protein1","protein2"),mode="left") %>%
      # select(starts_with("protein"),pred) %>%
      # mutate(pred=ifelse(is.na(pred),-1,pred))
    rebuilt_list[[complex]] <- all_combs %>% mutate(complex.names=complex)
  }
  rebuilt_data <- reduce(rebuilt_list, bind_rows) %>%
    group_by(protein1,protein2) %>%
    dplyr::summarise(complex.name=paste(complex.names,collapse=","), .groups="keep") %>%
    ungroup()
  
  
  data <- rebuilt_data
  
  return(data)

}
