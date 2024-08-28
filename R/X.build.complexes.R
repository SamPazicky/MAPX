#' X.build.complexes
#'
#' Using one of the network building algorithms, cluster initial network into individual protein complexes.
#' 
#' @param data Data frame with pair-wise interactions (columns 'protein1' and 'protein2') and a column with scores from ML model.
#' @param algo Character string: What algorithm should be used for clustering? Options are 'SP' for shuffle-predictions, 'SR' for shuffle-rows,
#' 'TR' for takeout-rows. Any other option (e.g. 'S') means that simple algorithm will be used. Default is 'S'.
#' @param scores.col Character string: Name of the columns with model scores.
#' @param init.stats  Logical: Should the initial network stats be calculated? Default is FALSE.
#' @param final.stats  Logical: Should the final network stats be calculated? Default is FALSE.
#' @param standard.set Data frame with columns protein1, protein2 and another column with labels.
#' @param labels.col Character string: Name of the columns with labels.
#' @param rep.steps Integer: How many times should the algorithm be repeated before averaging the results? Default is 100.
#' @param SP.shift Numeric: By how much can the predictions be randomly shifted in SP algorithm?
#' @param SP.finalpreds character: Which predictions should be use in final prediction calculation in SP algorithm?
#' 'randomized' for means of randomly changed predictions upon iteration, 'original' for original predictions.
#' @param SR.rows Numeric: Maximum by how many rows can a protein pair shift in order predictions in SR algorithm? Default is 20.
#' 
#' @import tidyverse
#' 
#' @return A list with three elements: $data with columns 'protein1', 'protein2' and score contains all data for plotting or further processing
#' of the network, 'stats_initial' contains statistics of the network before clustering. and 'stats_final' contains statistics of the network after
#' clustering.
#' @examples 
#' network.built <- X.build.complexes(data=cut_data, algo="SP",scores.col="score",init.stats=TRUE,final.stats=TRUE,
#'                   standard.set=GS,labels.col="complex")
#' @export
#' 
X.build.complexes <- function(
    data=NULL,
    algo="S",
    scores.col=NA,
    init.stats=FALSE,
    final.stats=FALSE,
    standard.set=NULL,
    labels.col=NA,
    rep.steps=100,
    SP.shift=0.2,
    SP.finalpreds="randomized",
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
  
  
  cat("Building the network of complexes...\n")
  
  # initiate the building
  grouped_data <- data %>% 
    dplyr::select(starts_with("protein"),score) %>%
    arrange(desc(score))
  complexes=list()
  building="inprogress"
  
  
  #
  # SHUFFLE COMPLEXES BY PREDS. Parameters: rep.steps, SP.shift
  # 
    
  if(algo=="SP") {
    
    pb <- txtProgressBar(min=0, max=rep.steps, style=3, initial="") # progress bar
    complexes_table <- list()
    
    for (itr in 1:rep.steps) {
      complexes=list()
      building="inprogress"
      
      
      # shuffle
      if(itr==1) {
        grouped_data_mod <- grouped_data
      } else {
        grouped_data_mod <- grouped_data %>%
          mutate(score=score+runif(nrow(.),(-1)*SP.shift,SP.shift)) %>%
          arrange(desc(score))
      }
      
      # cluster
      complexes_table[[itr]] <- X.cluster.complexes(data=grouped_data_mod,scores.col="score") 
      
      setTxtProgressBar(pb, itr)
    }
    close(pb)
    
    #average
    complexes_table_ref <- purrr::reduce(complexes_table,bind_rows,.id="itr") %>%
      group_by(protein1, protein2) %>%
      rename(itr_score=score)
    if(SP.finalpreds=="randomized") {
      complexes_table_ref <- complexes_table_ref %>%
        dplyr::summarise(count=n(), score=mean(itr_score,na.rm=TRUE), .groups="keep") %>%
        ungroup()
    } else if (SP.finalpreds=="original") {
      complexes_table_ref <- complexes_table_ref %>%
        dplyr::summarise(count=n(), .groups="keep") %>%
        ungroup() %>%
        cross_join(data%>%dplyr::select(starts_with("protein"),score), vars=paste0("protein",1:2),mode="left")
    } else {
      stop(paste0(SP.finalpreds," is not an option for SP.finalpreds."))
    }
    complexes_table_ref <- complexes_table_ref %>%
      mutate(prediction=count/rep.steps*score) %>%
      dplyr::select(starts_with("protein"),count,score,prediction)
    
    # final building step with simple algorithm
    grouped_data <- complexes_table_ref %>% dplyr::select(starts_with("protein"),prediction) %>% rename(score=prediction) %>% arrange(desc(score))
    complexes=list()
    building="inprogress"
    # code follows after the if-elses.
  
    #
    # SHUFFLE ROWS. Parameters: rep.steps, SR.rows
    #
    
    
  } else if (algo=="SR") {
    
    pb <- txtProgressBar(min=0, max=rep.steps, style=3, initial="") # progress bar
    complexes_table <- list()
    
    for (itr in 1:rep.steps) {
      complexes=list()
      building="inprogress"
      
      
      # shuffle
      if(itr==1) {
        grouped_data_mod <- grouped_data %>% rename(pseudoscore=score)
      } else {
        grouped_data_mod <- grouped_data %>% arrange(desc(score)) %>% 
          mutate(row_id=row_number()) %>% 
          mutate(random_no=round(runif(nrow(.),min=(-1)*SR.rows,max=SR.rows))) %>% 
          mutate(newrow_id=row_id+random_no) %>% mutate(pseudoscore=nrow(.)-newrow_id) %>% 
          dplyr::select(starts_with("protein"),pseudoscore)
      }
      
      #cluster
      complexes_table[[itr]] <- X.cluster.complexes(data=grouped_data_mod, scores.col="pseudoscore")
      
      setTxtProgressBar(pb, itr)
    }
    close(pb)
    
    #average
    complexes_table_ref <- purrr::reduce(complexes_table,bind_rows,.id="itr") %>%
      group_by(protein1, protein2) %>%
      rename(itr_score=score) %>%
      dplyr::summarise(count=n(), .groups="keep") %>%
      ungroup() %>%
      cross_join(data%>%dplyr::select(starts_with("protein"),score), vars=paste0("protein",1:2),mode="left") %>%
      mutate(prediction=count/rep.steps*score) %>%
      dplyr::select(starts_with("protein"),count,score,prediction)
    
    # final building step with simple algorithm
    grouped_data <- complexes_table_ref %>% dplyr::select(starts_with("protein"),prediction) %>% rename(score=prediction) %>% arrange(desc(score))
    complexes=list()
    building="inprogress"
    # code follows after the if-elses.
    
    
    #
    # TAKE OUT ROWS. Parameters: rep.steps
    #

  }  else if (algo=="TR") {
    
    if (rep.steps>nrow(grouped_data)) {
      stop(paste0("Number of iteration steps (rep.steps) is too high. The number can be max ", nrow(grouped_data), "for this dataset."))
    }
    
    pb <- txtProgressBar(min=0, max=rep.steps, style=3, initial="") # progress bar
    complexes_table <- list()
    folds=round(runif(nrow(grouped_data),1,rep.steps))
    cat("Iterating over left-out rows...\n")
    for (itr in 1:rep.steps) {
      
      #shuffle
      grouped_data_mod <- grouped_data %>% mutate(Fold=folds) %>% filter(Fold==itr) %>% dplyr::select(-Fold)
      
      #cluster
      complexes_table[[itr]] <- (X.cluster.complexes(data=grouped_data_mod,scores.col="score"))
      
      setTxtProgressBar(pb, itr)
    }
    close(pb)
    
    #average
    complexes_table_ref <- purrr::reduce(complexes_table,bind_rows,.id="itr") %>%
      group_by(protein1, protein2) %>%
      rename(itr_score=score) %>%
      dplyr::summarise(count=n(), score=mean(itr_score,na.rm=TRUE), .groups="keep") %>%
      ungroup() %>%
      mutate(prediction=count/rep.steps*score) %>%
      dplyr::select(starts_with("protein"),count,score,prediction)
    
    # final building step with simple algorithm
    grouped_data <- complexes_table_ref %>% dplyr::select(starts_with("protein"),prediction) %>% rename(score=prediction) %>% arrange(desc(score))
    complexes=list()
    building="inprogress"
    # code follows after the if-elses.
  }
  
  # simple algo 
  complexes_table <- X.cluster.complexes(data=grouped_data, scores.col="score")

  cat("Building the network of complexes... done. \n")
  
  #final statistics
  if(final.stats) {
    cat("Calculating final network statistics...")
    stats_final <- X.network.stats(data=complexes_table,standard.set=standard.set,labels.col=labels.col,complex.stats=TRUE)
    cat("\rCalculating final network statistics... done. \n")
    
  } else {
    stats_final <- NULL
  }
  
  output <- list("data"=complexes_table,"stats_initial"=stats_prior, "stats_final"=stats_final)
  return(output)
}
