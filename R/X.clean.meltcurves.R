#' X.clean.meltcurves
#'
#' Cleans the loaded mass spec data
#' @param data Data frame with the loaded mass spec data
#' @param keratin Logical: Should the keratins be removed?
#' @param trypsin Logical: Should trypsin be removed?
#' @param serum Logical: Should the serum proteins be removed?
#' @param average.duplicate.proteins Logical: Should the duplicate proteins be averaged?
#' @param remove.nas Logical: Should the rows with NA values be removed?
#' 
#' @import tidyverse
#' 
#' @return List of two elements. $data is the data frame with contaminants removed and $removed is the data frame listing the removed contaminants.
#' @examples 
#' MCdata.clean <- X.clean.meltcurves(MCdata.raw)
#' @export


X.clean.meltcurves <- function(
    data=NULL,
    keratins=TRUE,
    serum=TRUE,
    trypsin=TRUE,
    average.duplicate.proteins=TRUE,
    remove.nas=TRUE
) {
  
  if(is.null(data)) {
    stop("Please include data argument")
  } else if (!"description" %in% names(data) | !"protein" %in% names(data) | !"length" %in% names(data)) {
    stop("Data must contain columns 'protein', 'description' and 'length'.")
  }
  
  removed <- list()
  
  # remove keratins 
  if(keratins) {
    removed$keratins <- data %>% filter(grepl("Keratin",description)) %>% pull(protein)
    data <- data %>% filter(!grepl("Keratin", description))
    cat(paste0("Removed ",length(removed$keratins)," keratins from the data. \n"))
    
  }
  
  # remove serum albumin
  if(serum) {
    removed$albumin <-  data %>% filter(grepl("Serum albumin",description)) %>% pull(protein)
    data <- data %>% filter(!grepl("Serum albumin", description))
    cat(paste0("Removed ",length(removed$albumin)," albumin rows from the data. \n"))
    
  }
  
  # remove trypsin
  if(trypsin) {
    removed$trypsin <- data %>% filter(grepl("Trypsin",description)) %>% pull(protein)
    data <- data %>% filter(!grepl("Trypsin", description))
    cat(paste0("Removed ",length(removed$trypsin)," trypsin rows from the data. \n"))
    
  }
  
  # remove duplicate proteins
  if(average.duplicate.proteins) {
    removed$duplicates <- data$protein %>% table() %>% stack() %>% filter(values>1) %>% pull(ind) %>% as.character()
    oldata <- data
    data <- data %>%
      group_by(protein) %>%
      dplyr::summarise(across(where(is.numeric), ~ mean(.x))) %>%
      ungroup()
    rmd <- nrow(oldata)-nrow(data)
    cat(paste0("Averaged ",rmd," duplicate protein rows. \n"))
  }
  
  if(remove.nas) {
    removed$nas <- data %>% filter(if_any(everything(), is.na)) %>% pull(protein)
    oldata <- data
    data <- data %>%
      na.omit()
    rmd <- nrow(oldata)-nrow(data)
    cat(paste0("Removed ",rmd," rows with NAs.\n"))
  }
  output <- list("data"=data, "removed"=removed)
  return(output)
  
}
