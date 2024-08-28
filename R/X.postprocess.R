#' X.postprocess
#'
#' Using one of the network building algorithms, cluster initial network into individual protein complexes.
#' 
#' @param data Data frame with pair-wise interactions (columns 'protein1' and 'protein2') and a column with probabilities from ML model.
#' @param mode Character string:  'trim' for trimming of not-well-connected proteins and 'split' for separating of two loosely connected subunits.
#' Default is 'trim'.
#' @param scores.col Character string: Name of the columns with prediction values.
#' @param weight Logical: Should the prediction values be taken into account for the postprocessing? Default is FALSE
#' @param zcut Numeric: Sub-units with z-score above this value will be removed. Default is 2.
#' @param init.stats  Logical: Should the initial network stats be calculated? Default is FALSE.
#' @param final.stats  Logical: Should the final network stats be calculated? Default is FALSE.
#' @param standard.set Data frame with columns protein1, protein2 and another column with labels.
#' @param labels.col Character string: Name of the columns with labels.
#' 
#' @import tidyverse
#' 
#' @return A list with three elements: $data with columns 'protein1', 'protein2' and score contains all data for plotting or further processing
#' of the network, 'stats_initial' contains statistics of the network before clustering. and 'stats_final' contains statistics of the network after
#' clustering.
#' @examples 
#' network.post.split <- X.postprocess(data=network.refined$data, mode="split", scores.col="score", final.stats=TRUE, 
#'         standard.set=GS,labels="complex", weighted=FALSE)
#'
#' network.post.trim <- X.postprocess(data=network.post.split$data, mode="trim", scores.col="score", final.stats=TRUE, 
#'         standard.set=GS,labels="complex", weighted=TRUE)
#' @export
#
X.postprocess = function(
    data=NULL, 
    mode="trim", # trim or split
    scores.col=NA,
    weighted=FALSE,
    zcut=2,
    init.stats=FALSE,
    final.stats=FALSE,
    standard.set=NULL,
    labels.col=NA
)
{
  

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
  
  #initiate variables
  
  complexes <- list()
  building=TRUE
  i=1
  
  
  if(init.stats) {
    cat("Calculating initial network statistics...")
    stats_prior <- X.network.stats(data, standard.set=standard.set,labels.col=labels.col, complex.stats=FALSE)
    cat("\rCalculating initial network statistics... done. \n")
  } else {
    stats_prior <- NULL
  }
  
  
  cat("Post-processing complexes...\n")
  
  if(mode=="trim") {
    
    complexes <- X.pairwise.to.complexes(data,names="complex.name")
    
    connections_table <- data.frame()
    
    for (comp in seq_along(complexes)) {
      complex_data <- data %>% filter(protein1 %in% complexes[[comp]] & protein2 %in% complexes[[comp]])
      proteins <- complex_data %>% select(c("protein1","protein2")) %>% unlist() %>% unname() %>% unique()
      
      for (p in proteins) {
        protein_no <- complex_data %>% filter(protein1==p | protein2==p) %>% nrow()
        protein_prob <- complex_data %>% filter(protein1==p | protein2==p) %>% summarise(prob=mean(score)) %>% unlist() %>% unname()
        protein_no_subtable <- data.frame("Complex"=comp,"Protein"=p,"Connections"=protein_no,"Prob"=protein_prob)
        connections_table <- bind_rows(connections_table,protein_no_subtable)
      }
      
    }
    
    looped_table <- connections_table %>%
      mutate(ConProb=Connections*Prob)
    looping=TRUE
    if(weighted) {
      evalvar="ConProb"
    } else {
      evalvar="Connections"
    }
    rown=nrow(looped_table)
    
    while(looping) {
      looped_table <- looped_table %>%
        group_by(Complex) %>%
        mutate(z_score=as.numeric(scale(!!sym(evalvar)))) %>%
        ungroup() %>%
        mutate(z_score=ifelse(is.nan(z_score),0,z_score)) %>%
        filter(z_score>=(-1)*zcut)
      if(nrow(looped_table)==rown) {
        looping=FALSE
      }
      rown = nrow(looped_table)
    }
    
    filtered_proteins <- anti_join(connections_table,looped_table,by="Protein") %>% pull(Protein)
    
    filtered_complexes <- list()
    splitd_complexes <- list()
    j=1
    for (i in seq_along(complexes)) {
      filtered_complexes[[i]] <- setdiff(complexes[[i]], filtered_proteins)
      if( length(complexes[[i]]) - length(filtered_complexes[[i]]) >= 2 ) {
        splitd_complexes[[j]] <- intersect(filtered_proteins,complexes[[i]])
        j=j+1
      }
    }
    
    filtered_complexes <- c(filtered_complexes,splitd_complexes)
    
    complexes_table <- X.complexes.to.pairwise(filtered_complexes) %>%
      cross_join((data%>%dplyr::select(!complex.name)),vars=c("protein1","protein2"),mode="left") %>%
      filter(!is.na(score)) %>%
      dplyr::select(starts_with("protein"),complex.name,all_of(scores.col)) %>%
      as.data.frame()
      
  
    
  } else if(mode=="split") {
    
    # rebuild complexes from table
    complexes <- X.pairwise.to.complexes(data,scores.col="score",names="complex.name")

    orig_complexes <- complexes
    trimming=TRUE
    do_not_filter <- c()
    while (trimming) {
      
      # constitute connections_table: for each complex, each protein connects to some number of proteins (column Connections)
      # and has some average probability ("Prob") and when I multiply those, I get column "ConProb".
      # If weighted=TRUE than ConProb is taken into account, otherwise just Connections.
      
      connections_table <- data.frame()
      
      for (comp in seq_along(complexes)) {
        complex_data <- data %>% filter(protein1 %in% complexes[[comp]] & protein2 %in% complexes[[comp]])
        proteins <- complex_data %>% select(c("protein1","protein2")) %>% unlist() %>% unname() %>% unique()
        
        for (p in proteins) {
          protein_no <- complex_data %>% filter(protein1==p | protein2==p) %>% nrow()
          protein_prob <- complex_data %>% filter(protein1==p | protein2==p) %>% summarise(prob=mean(score)) %>% unlist() %>% unname()
          protein_no_subtable <- data.frame("Complex"=comp,"Protein"=p,"Connections"=protein_no,"Prob"=protein_prob)
          connections_table <- bind_rows(connections_table,protein_no_subtable)
        }
      }
      
      looped_table <- connections_table %>%
        mutate(ConProb=Connections*Prob)
      looping=TRUE
      if(weighted) {
        evalvar="ConProb"
      } else {
        evalvar="Connections"
      }
      
      
      # now for each complex calculate a z score for each protein - either from number of Connections (non-weighted)
      # or from ConProb (weighted). Proteins with z-score above zcut will be further processed.
      
      rown=nrow(looped_table)
      leftover_table <- data.frame()
      while(looping) {
        looped_table <- looped_table %>%
          group_by(Complex) %>%
          mutate(z_score=as.numeric(scale(!!sym(evalvar)))) %>%
          ungroup() %>%
          mutate(z_score=ifelse(is.nan(z_score),0,z_score))
        

        leftover_table <- bind_rows(leftover_table, looped_table %>% filter(z_score>zcut) )
        looped_table <- looped_table %>% filter(z_score<=zcut)
        
        if(nrow(looped_table)==rown) {
          looping=FALSE
        }
        rown = nrow(looped_table)
      }
      

      filtered_proteins <- leftover_table %>% arrange(desc(z_score)) %>% pull(Protein)
      filtered_proteins <- setdiff(filtered_proteins, do_not_filter)
      
      cat("\rPost-processing complexes...",length(filtered_proteins),"proteins remaining...\n")
      if(length(filtered_proteins)==0) {
        trimming=FALSE
        break
      }
      
      # taking out highest-rated outlier protein, rebuild that particular complex
      
      i=filtered_proteins[1]
      
      which.complex <- leftover_table %>% filter(Protein==i) %>% pull(Complex)
      complex <- complexes[[which.complex]]
      filtered_complex <- setdiff(complex,i)
      grouped_data <- data %>% 
        filter(protein1 %in% filtered_complex & protein2 %in% filtered_complex) %>%
        select(starts_with("protein"),all_of(scores.col))
      split_complexes<-list()
      building=TRUE
      while(building) {
        FMs <- grouped_data %>% slice_head(n=1) %>% select(-score) %>% unlist() %>% as.vector()
        GPs <- grouped_data %>% filter(protein1 %in% FMs | protein2 %in% FMs) %>% select(-score) %>%
          t() %>% c() %>% unique() %>% .[!. %in% FMs]
        if(length(GPs)==0) {
          split_complexes[[paste(FMs,collapse=";")]] <- FMs
          grouped_data <- grouped_data %>% slice(-1)
          if(nrow(grouped_data)==0) {
            building=FALSE
          }
          next
        }
        
        split_complexes[[paste(FMs,collapse=";")]] <- c(FMs,GPs)
        grouped_data <- grouped_data %>%
          filter(!protein1 %in% c(FMs,GPs)) %>%
          filter(!protein2 %in% c(FMs,GPs))
        
        
        if(nrow(grouped_data)==0) {
          building=FALSE
        }
      }
      
      # if the complex was not split, it is unsplittable by this protein and this protein should not be considered anymore.
      if(length(split_complexes)==1) {
        do_not_filter <- c(do_not_filter, i)
        next
      }
      
      # now add the taken-out protein to the complex which it has higher average probability to
      meanvalues<-c()
      for (sc in seq_along(split_complexes)) {
        sclist <- unlist(split_complexes[sc]) %>% unname()
        meanvalues[sc] <- data %>%
          filter( (protein1 %in% sclist & protein2==i) | (protein2 %in% sclist & protein1==i) )%>%
          summarise(mean=mean(score)) %>% unlist() %>% unname()
        
      }
      split_complexes[[which.max(meanvalues)]] <- c(split_complexes[[which.max(meanvalues)]],i)
      complexes <- complexes[-which.complex]
      complexes <- c(complexes,split_complexes)
      
    }
    
    complexes_table <- X.complexes.to.pairwise(complexes) %>%
      cross_join((data%>%dplyr::select(!complex.name)),vars=c("protein1","protein2"),mode="left") %>%
      filter(!is.na(score)) %>%
      dplyr::select(starts_with("protein"),complex.name,all_of(scores.col)) %>%
      as.data.frame()
    
  }
  
  cat("\rPost-processing complexes... done.                 \n")
  
  
  if(final.stats) {
    cat("Calculating final network statistics...")
    stats_final <- X.network.stats(data=complexes_table,standard.set=standard.set,labels.col=labels.col,complex.stats=TRUE)
    cat("\rCalculating final network statistics... done. \n")
    
  } else {
    stats_final <- NULL
  }
  
  output <- list("data"=complexes_table, "stats_initial"=stats_prior, "stats_final"=stats_final)
  return(output)
  
}
