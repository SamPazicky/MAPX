#' X.network.stats
#'
#' Evaluate a network of protein complexes
#' 
#' @param data Data frame with columns 'protein1' and 'protein2' and pairwise prediction values.
#' @param scores.col Character string: Name of the data column with model scores.
#' @param standard.set Data frame with columns protein1, protein2 and another column with labels.
#' @param labels.col Character string: Name of the columns with labels.
#' @param complex.stats Logical: Should the statistics specific for protein complexes be calculated? Do not use with unclustered networks. Default is FALSE.
#' @param interconnected Logical: If TRUE, the statistics will be calculated on interconnected network, meaning that even the subunits
#' that are present in the same cluster, but did not pass the prediction cut-off, will be considered interconnected. Default is FALSE.
#' @import tidyverse
#' 
#' @return A list with two elements: $data with underlying data and $plot with the plotted statistics.
#' @examples 
#' network.stats <- X.network.stats(data)
#' @export
#' 

X.network.stats = function(
    data=NULL,
    scores.col=NA,
    standard.set=NULL,
    labels.col=NA,
    complex.stats=FALSE,
    interconnected=FALSE)
{

  customPlot <- list(
    theme_bw(base_size = 12),
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.margin=ggplot2::margin(5,5,5,5, "pt")
    )
  )
  
  if(is.null(data)) {
    stop("Please specify data.")
  } else {
    if(!"protein1" %in% names(data) | !"protein2" %in% names(data)) {
      stop("The data must contain columns protein1, protein2 and a column with model scores.")
    }
    data <- as.data.frame(data)
    if(is.na(scores.col)) {
      scores.col <- names(data)[ncol(data)]
    }
    data <- data %>% rename(score=!!sym(scores.col))
  }
  
  
  if(is.null(standard.set)) {
    stop("Include standard.set.")
  } else {
    standard.set <- as.data.frame(standard.set)
    if(is.na(labels.col)) {
      labels.col <- names(standard.set)[ncol(standard.set)]
    }
    
  }

  
  if(interconnected) {
    intdata <- X.interconnect.network(data, scores.col="score")
    data <- cross_join(data,intdata, vars=c("protein1","protein2"), mode="left")
  }
  
  stats_table <- data %>% 
    cross_join(standard.set,vars=c("protein1","protein2"),mode="left")
  
  
  TPs <- stats_table %>% filter(!!sym(labels.col)==1) %>% nrow()
  Ps <- stats_table %>% filter(!is.na(!!sym(labels.col))) %>% nrow()
  TNs <- Ps-TPs
  New <- stats_table %>% filter(is.na(complex)) %>% nrow()
  
  stats <- data.frame(
    Interaction=c("TPs","TNs","new"),
    Number=c(TPs,TNs,New)
  )
 
  stats1data <- stats %>%
    column_to_rownames("Interaction")%>% t() %>% as.data.frame() %>%
    setNames(stats$Interaction) %>%
    mutate(recall=TPs/Ps) %>%
    mutate(precision=TPs/(TPs+TNs)) %>%
    mutate(stringency=new/(TPs+TNs+new))
  
  
  statplot1 <- stats %>% 
    mutate(Interaction=factor(Interaction,levels=c("TNs","TPs","new"))) %>%
    ggplot(aes(x=Interaction,y=Number)) + geom_bar(mapping=aes(fill=Interaction), stat="identity") +
    geom_text(aes(label=Number), position=position_dodge(width=0.9), vjust=-0.25) + 
    scale_y_continuous(limits=c(0,max(stats$Number)*1.1), expand=c(0,0)) +
    customPlot +
    theme(legend.position="none")

  # to fill in complex statistics
  if(complex.stats) {

    complexes <- X.pairwise.to.complexes(data)
    stand.complexes <- standard.set %>% filter(!!sym(labels.col)==1) %>% X.pairwise.to.complexes()

    #fragmentation
    # how it works: for each standard complex, the algo detects to how many cluster it is divided (e.g. 20 subunits to c(10,5,5)).
    # each of these numbers is divided by number of subunits detected c(0.5,0.25,0.25), of that, power of two is calculated so that
    # the larger subunits have relatively more weight (c(0.25,0.0625,0.0625))
    
    fragment.data <- data.frame()
    
    # for each standard complex
    for(stand.cpx in 1:length(stand.complexes)) {
      
      # check how many subunits of that standard complex are in each detected complex
      cur.vec <- c()
      for(cpx in 1:length(complexes)) {
        cur.vec <- c(cur.vec,length(intersect(stand.complexes[[stand.cpx]],complexes[[cpx]])))
      }
      
      # put that vector in data frame (fragment.data) - one standard complex per column
      if(nrow(fragment.data)==0) {
        fragment.data <- as.data.frame(cur.vec)%>%rename(!!sym(as.character(stand.cpx)):=cur.vec)
      } else {
        fragment.data <- bind_cols(fragment.data, 
                                   as.data.frame(cur.vec)%>%rename(!!sym(as.character(stand.cpx)):=cur.vec)
        )
      }
      
    }

    fragment.scores <- c()
    complex.sizes <- stand.complexes %>% sapply(length)
    real.sizes <- c()
    
    #go through each column of fragment.data (each column is fragmentation of a complex)
    for(frag.col in 1:ncol(fragment.data)) {
      # isolate non-zero values of the column
      vector <- fragment.data[[frag.col]] %>% .[(.)!=0]
      if(length(vector)==0) {
        next
      }
      # calculate fragment.score for the column: each fragment divided by total subunits, and then from those calculate weighted mean
      fragment.scores[as.character(frag.col)] <- vector %>% sapply(function(x) x/sum(vector)) %>% weighted.mean(w=(.))
      real.sizes[as.character(frag.col)] <- sum(fragment.data[[frag.col]])
    }
    complex.sizes <- complex.sizes[as.numeric(names(real.sizes))]
    results.table <- data.frame(
      fragment.score=fragment.scores,
      size=complex.sizes,
      real.size=real.sizes
    ) %>% 
      na.omit() %>%
      mutate(retrieval=real.size/size)
    
    
    
    
    #purity
    
    purity.scores <- c()
    #loop through standard complexes
    for(stand.cpx in 1:length(stand.complexes)) {
      
      # calculate raw purity for standard complex: for each calculated complex that contains subunits of standard complex, what percentage of nodes are standard complex?
      raw.purity <- lapply(complexes, function(x) intersect(x,stand.complexes[[stand.cpx]])%>%length()/length(x))%>%unlist() %>% .[(.)!=0]
      # what portion of the standard complex is in each calculated complex
      fragment.ratio <- lapply(complexes, function(x) intersect(x,stand.complexes[[stand.cpx]])%>%length()/length(stand.complexes[[stand.cpx]]))%>%unlist() %>% .[(.)!=0]
      
      # weight the raw purity scores by fragment.ratio: The smaller the number of standard subunits in a complex, the less important it is for purity
      # score, so powering by the fragment ratio will result in higher scores for lower-fragment-ratio calculated complexes.
      purity.scores[as.character(stand.cpx)] <- raw.purity^fragment.ratio %>% mean()
    }
 
    results.table$purity.score <- purity.scores[!is.na(purity.scores)]
    results.table <- results.table %>%
      filter(real.size>1) %>%
      filter(!(size==2&real.size==2&fragment.score!=1))
    
    fragment.score <- mean(results.table$fragment.score)
    purity.score <- mean(results.table$purity.score)
    retrieval.score <- mean(results.table$retrieval)
    
    network.data=data.frame(fragmentation=fragment.score,purity=purity.score,retrieval=retrieval.score)
    
    output=list(data=stats1data,plot=statplot1,network.stats=network.data, complex.stats=results.table)
  } else {
    output=list("data"=stats1data, "plot"=statplot1)
    
  }
  
  return(output)
  
}
