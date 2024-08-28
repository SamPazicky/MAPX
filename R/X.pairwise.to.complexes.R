#' X.pairwise.to.complexes
#'
#' Transform a pairwise interaction table into a list of protein complexes.
#' 
#' @param data A data.frame with pairwise interaction table with column names 'protein1' and 'protein2'.
#' @param names.col Character string: Name of the column with complex names.
#' @return A named list with complex names as the list names and subunits as vector elements.
#' @examples 
#' pairwise.table <- X.complexes.to.pairwise(list_of_complexes)
#' @export
#'
 
X.pairwise.to.complexes <- function(
    data=NULL,
    names.col=NA,
    scores.col=NULL
) {
  
  if(is.null(data)) {
    stop("Please specify data.")
  } else {
    if(!"protein1" %in% names(data) | !"protein2" %in% names(data)) {
      stop("The data must contain columns protein1 and protein2.")
    }
    if(is.na(names.col)) {
      data <- as.data.frame(data) %>%
        select(protein1,protein2)
    } else {
      data <- as.data.frame(data) %>%
        select(protein1,protein2,all_of(names.col))
    }

  }
  
  building=TRUE
  complexes=list()
  i=1
  
  # rebuild complexes from table
  while(building) {
    
    proteins <- data %>% slice_head(n=1) %>% select(c("protein1","protein2")) %>% unlist() %>% unname()
    if(!is.na(names.col)) {
      complexname = data %>% slice_head(n=1) %>% dplyr::select(all_of(c(names.col))) %>% unlist() %>% unname()
    } else {
      complexname = i
    }
    data <- data %>% slice(-1)
    
    adding=TRUE
    current_length=0
    while(adding) {
      
      group_proteins <- data %>% filter(protein1 %in% proteins | protein2 %in% proteins) %>%
        select(c("protein1","protein2")) %>% unlist() %>% unname() %>% union(proteins) %>% unique()
      if(length(group_proteins)!=0) {
        proteins <- group_proteins
        if(length(proteins)>current_length) {
          current_length = length(proteins)
        } else {
          adding=FALSE
        }
      } else {
        adding=FALSE
      }
    }
    
    data <- data %>% filter(!protein1 %in% proteins & !protein2 %in% proteins)
    
    
    complexes[[complexname]] <- proteins
    
    if(nrow(data)==0) {
      building=FALSE
    }
    i=i+1
  }
  
  return(complexes)
}
