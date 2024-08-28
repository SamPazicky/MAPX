#' X.assemble.negstandard
#'
#' Using a database table with Entry names and GO terms, assemble a table of putatively non-interacting proteins.
#' 
#' @param data Data frame with two columns: First for protein IDs and second for separated GO terms.
#' @param sep Separator for the GO terms in the second column
#' @param cutoff Numeric: How many GO terms can a pair of proteins share to still be included in the list of 
#' non-interacting proteins?
#' 
#' @import tidyverse
#' 
#' @return A data frame with pair-wise negative standard table for protein-protein interactions.
#' @examples 
#' GSneg <- X.assemble.negstandard(data=GOterms, sep=";")
#' @export

X.assemble.negstandard <- function(
  data=NULL,
  sep=";",
  cutoff=0
){
  

  if(is.null(data)) {
    stop("Input database.")
  } else {
    data <- data %>%
      dplyr::select(1,2) %>%
      setNames(c("protein","GOs")) %>%
      arrange(protein)
  }
  
  all_proteins <- data %>%
    pull(protein) %>%
    unique()
  
  all_GOs <- data %>%
    mutate(GO=str_remove_all(GOs,"\\s+")) %>%
    pull(GOs) %>%
    paste(collapse=sep) %>%
    str_split(sep) %>% 
    unlist() %>%
    unique() %>%
    .[(.)!=""]
  
  j=1
  cat("Assigning GO terms...\n")
  pb <- txtProgressBar(min=0, max=length(all_GOs), style=3, initial="") # progress bar
  for (GO in all_GOs) {
    data <- data %>% mutate(!!GO:=
                                    ifelse(grepl(GO, GOs),
                                           1,
                                           0))
    setTxtProgressBar(pb,j)
    j=j+1
  }
  close(pb)
  
  excluded <- data %>% filter(GOs=="") %>% pull(protein)
  included <- setdiff(all_proteins,excluded)
  
  mdata <- data %>% 
    column_to_rownames("protein") %>%
    filter(GOs!="") %>%
    dplyr::select(!GOs)
  
  intlist <- list()
  cat("Listing potential interactors...\n")
  pb <- txtProgressBar(min=0, max=length(colnames(mdata)), style=3, initial="") # progress bar
  
  for(go in colnames(mdata)) {
    intlist[[go]] <- mdata %>% 
      dplyr::filter((!!sym(go))==1) %>%
      rownames()
    setTxtProgressBar(pb,which(colnames(mdata)==go))
  }
  close(pb)
  
  cat("Assembling table of non-interacting proteins...\n")
  inttable <- intlist %>%
    MAPX::X.complexes.to.pairwise() %>%
    mutate(shared_GOs=str_count(complex.name,",")+1) %>%
    filter(shared_GOs>cutoff)
  
  alltable <- expand.grid.unique(included,included) %>%
    as.data.frame() %>%
    setNames(paste0("protein",1:2)) %>%
    mutate(complex=0)
  
  filtered_add_data <- inttable %>%
    dplyr::select(starts_with("protein")) %>%
    mutate(complex=1) %>%
    bind_rows(alltable) %>%
    group_by(protein1,protein2) %>%
    dplyr::summarise(scomplex=sum(complex,na.rm=TRUE)) %>%
    ungroup() %>%
    filter(scomplex==0) %>%
    dplyr::select(starts_with("protein"))
  
  cat("The resulting pairwise negative standard set is composed of", nrow(filtered_add_data), "pairs.\n")
  return(filtered_add_data)
  
}
