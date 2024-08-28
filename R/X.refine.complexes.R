#' X.refine.complexes
#'
#' Using one of the network building algorithms, cluster initial network into individual protein complexes.
#' 
#' @param data Data frame with pair-wise interactions (columns 'protein1' and 'protein2') and a column with probabilities from ML model.
#' @param algo Character string: What algorithm should be used for refining? Options are 'SCBR' for shuffling-complexes-by-rows or
#' 'SCBP' for shuffling complexes by prediction value. Default is 'SCBP'.
#' @param scores.col Character string: Name of the columns with prediction values.
#' @param init.stats  Logical: Should the initial network stats be calculated? Default is FALSE.
#' @param final.stats  Logical: Should the final network stats be calculated? Default is FALSE.
#' @param standard.set Dataframe with columns protein1, protein2 and another column with labels.
#' @param labels.col Character string: Name of the columns with labels.
#' @param rep.steps Integer: How many times should the algorithm be repeated before averaging the results? Default is 100.
#' @param SP.shift Numeric: By how much can the predictions be randomly shifted in SCBP algorithm?
#' @param SR.rows Numeric: Maximum by how many rows can a protein pair shift in order predictions in SR algorithm? Default is 20.
#' 
#' @import tidyverse
#' 
#' @return A list with three elements: $data with columns 'protein1', 'protein2' and score contains all data for plotting or further processing
#' of the network, 'stats_initial' contains statistics of the network before clustering. and 'stats_final' contains statistics of the network after
#' clustering.
#' @examples 
#' X.refine.complexes(network.built)
#' @export
#' 
X.refine.complexes <- function(
    data=NULL,
    algo="SCBP",
    scores.col=NA,
    init.stats=FALSE,
    final.stats=FALSE,
    standard.set=NULL,
    labels.col=NA,
    rep.steps=100,
    SP.shift=0.2,
    SR.rows=20
) {

  if(is.null(data)) {
    stop("Please specify data.")
  } else {
    if(is.na(scores.col)) {
      scores.col <- names(data)[ncol(data)]
      cat("No scores column selected. Selecting the last column:",scores.col,"\n")
      
    }
    data <- as.data.frame(data) %>%
      rename(score=!!sym(scores.col))
  }
  
  
  if(init.stats|final.stats) {
    if(is.null(standard.set)) {
      stop("If evaluate is TRUE, include standard.set.")
    } else {
      standard.set <- as.data.frame(standard.set)
      if(is.na(labels.col)) {
        labels.col <- names(standard.set)[ncol(standard.set)]
        cat("No labels column selected. Selecting the last column:",labels.col,"\n")
      }
    }
  }
  
  # prior stats
  
  if(init.stats) {
    cat("Calculating initial network statistics... ")
    stats_prior <- X.network.stats(data, standard.set=standard.set,labels.col=labels.col, complex.stats=FALSE)
    cat("\rCalculating initial network statistics... done. \n")
  } else {
    stats_prior <- NULL
  }
  
  cat("Refining protein complex network...\n")
  
  # initiate the building
  grouped_data <- data %>% 
    select(starts_with("protein"),score) %>%
    arrange(desc(score))
  complexes=list()
  building="inprogress"
  
  if(algo=="SCBR") {
    
    complexes <- X.pairwise.to.complexes(data)
    
    all.complexes.table.rearranged <- list()
    for(cmplx in seq_along(complexes)) {
      
      cat("Iterating over complex",cmplx,"/",length(complexes),"...\n")
      
      grouped_data_cmplx <- data %>% 
        filter(if_all(starts_with("protein"), ~ . %in% complexes[[cmplx]])) %>%
        dplyr::select(starts_with("protein"),score)
      
      complexes.table.rearranged <- list()
      pb <- txtProgressBar(min=0, max=rep.steps, style=3, initial="") # progress bar
      for (itr in 1:rep.steps) {
        building="inprogress"
        
        grouped_data_mod <- grouped_data_cmplx %>%
          arrange(desc(score)) %>% 
          mutate(row_id=row_number()) %>% 
          mutate(random_no=round(runif(nrow(.),min=(-1)*SR.rows,max=SR.rows))) %>% 
          mutate(newrow_id=row_id+random_no) %>% mutate(pseudoscore=nrow(.)-newrow_id) %>% 
          dplyr::select(starts_with("protein"),pseudoscore)
        
        complexes.table.rearranged[[itr]] <- X.cluster.complexes(data=grouped_data_mod, scores.col="pseudoscore")
      }
      #average
      all.complexes.table.rearranged[[cmplx]] <- purrr::reduce(complexes.table.rearranged,bind_rows,.id="itr") %>%
        group_by(protein1, protein2) %>%
        rename(itr_score=score) %>%
        dplyr::summarise(count=n(), .groups="keep") %>%
        ungroup() %>%
        cross_join(data%>%dplyr::select(starts_with("protein"),score), vars=paste0("protein",1:2),mode="left") %>%
        mutate(prediction=count/rep.steps*score) %>%
        dplyr::select(starts_with("protein"),count,score,prediction)
      
      setTxtProgressBar(pb, itr)
      cat("\n")
    }
    close(pb)
    
    complexes_table <- purrr::reduce(all.complexes.table.rearranged, bind_rows)
    grouped_data <- complexes_table %>% 
      dplyr::select(starts_with("protein"),prediction) %>% rename(score=prediction) %>% arrange(desc(score))

    #
    # SHUFFLE COMPLEXES BY scores. Parameters: rep.steps, SP.shift
    # 
    
  } else if(algo=="SCBP") {
    
    complexes=X.pairwise.to.complexes(data,scores.col="score")

    all.complexes.table.rearranged <- list()
    for(cmplx in seq_along(complexes)) {
      cat("Iterating over complex",cmplx,"/",length(complexes),"...\n")
      grouped_data_cmplx <- data %>% 
        filter(if_all(starts_with("protein"), ~ . %in% complexes[[cmplx]])) %>%
        dplyr::select(starts_with("protein"),score)
      
      complexes.table.rearranged <- list()
      pb <- txtProgressBar(min=0, max=rep.steps, style=3, initial="") # progress bar
      for (itr in 1:rep.steps) {
        building="inprogress"
        
        grouped_data_mod <- grouped_data_cmplx %>%
          mutate(score=score+runif(nrow(.),(-1)*SP.shift,SP.shift)) %>%
          arrange(desc(score))
        
        #cluster
        complexes.table.rearranged[[itr]] <- X.cluster.complexes(grouped_data_mod,scores.col="score")
        
        setTxtProgressBar(pb, itr)
      }
      close(pb)
      
      #average
      all.complexes.table.rearranged[[cmplx]] <- purrr::reduce(complexes.table.rearranged,bind_rows,.id="itr") %>%
        group_by(protein1, protein2) %>%
        rename(itr_score=score) %>%
        dplyr::summarise(count=n(), score=mean(itr_score,na.rm=TRUE), .groups="keep") %>%
        ungroup() %>%
        mutate(prediction=count/rep.steps*score) %>%
        dplyr::select(starts_with("protein"),count,score,prediction)
    }
    
    complexes_table <- purrr::reduce(all.complexes.table.rearranged, bind_rows)
    grouped_data <- complexes_table %>% 
      dplyr::select(starts_with("protein"),prediction) %>% rename(score=prediction) %>% arrange(desc(score))
    
  }
  
  #final building step
  
  complexes_table <- X.cluster.complexes(data=grouped_data,scores.col="score")
  
  cat("Refining protein complex network...done.\n")
  if(final.stats) {
    cat("Calculating final network statistics...")
    stats_final <- X.network.stats(data=complexes_table,standard.set=standard.set,labels.col=labels.col,complex.stats=FALSE)
    cat("\rCalculating final network statistics... done. \n")
    
  } else {
    stats_final <- NULL
  }

  output=list(data=complexes_table, stats_initial=stats_prior,stats_final=stats_final)
  return(output)
}



