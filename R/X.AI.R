#' X.AI
#'
#' Calculate assembly indeces from complex features reduce into a PCA plot.
#' 
#' @param data Data frame with features reduce by PCA. Columns are 'protein', and PC columns starting with 'PC'. Additional columns
#' 'condition' and 'replicate' can be used.
#' @param complexes Named list: A list where each element is vector with proteins and the name of that element is the name of the protein complex.
#' @param quality.data Data.frame with a column for protein and 
#' @param noPCs Integer: How many principal components should be used for the assembly index calculation?
#' @param trials Integer: How many random Euclidean distances should be calculated for the statistics?
#' @param cluster.cut Numeric: At what height should the clustered data from each complex be cut? The largest cluster will be used to calculate
#' the complex centroid, distance from which is further calculated and gives
#' @param weights Numeric vector: Weights for averaging principal components for Euclidean distance calculations.
#' @param y.trans Logical: Should the assembly index be log2 transformed? Default is TRUE.
#' 
#' @import tidyverse
#' @import caret
#' 
#' @return PCA data and different supporting information.
#' @examples 
#' all.fitdata.frame <- data.frame()
#' for(tp in timepoints) {
#'   for(rep in reps) {
#'     all.fitdata.frame <- bind_rows(all.fitdata.frame,
#'                                    all.fitdata[[tp]][[rep]] %>% mutate(replicate=as.character(rep)) %>% mutate(condition=tp)
#'     )
#'   }
#' }
#' quality.data <- all.fitdata.frame %>%
#'   dplyr::select(protein,condition,replicate,R2) %>%
#'   rename(quality=R2)
#'
#' data.AIs <- X.AI(data.PCAs$data,complexes, quality.data)
#' 
#' @export
#' 

X.AI <- function(
    data=NULL,
    complexes=list(),
    quality.data=NULL,
    noPCs=0,
    trials=500,
    cluster.cut=0.65,
    weights=NULL,
    y.trans=TRUE
    
) {
  
  customPlot <- list(
    theme_bw(base_size = 12),
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.margin=ggplot2::margin(5,5,5,5, "pt")
    )
  )
  
  if(is.null(data)) {
    stop("Please include data")
  } else {
    data <- data %>% as.data.frame()
    if(!any(grep("PC",names(data)))) {
      stop("The data must contain at least one PC column. Name the PC columns PC1, PC2....PCn")
    }
    if(!"protein" %in% names(data)) {
      stop("The data must contain a 'protein' column.")
    }
    if(!"replicate" %in% names(data)) {
      data <- data %>% mutate(replicate="1")
    }
    replicates <- unique(data$replicate)
    if(!"condition" %in% names(data)) {
      data <- data %>% mutate(condition="1")
    }
    conditions <- unique(data$condition)
  }
  if(length(complexes)<1) {
    stop("Please include a named list with complex subunits")
  }
  if(is.null(quality.data)) {
    quality.data <- data %>% dplyr::select(protein,condition,replicate) %>% mutate(quality=1)
  } else {
    if(!"quality" %in% names(quality.data)) {
      stop("Quality data must contain the column 'quality.")
    } else {
      if(any(is.na(quality.data$quality))) {
        cat("Some quality values are NA. Rewriting them to 0.001.\n")
        quality.data <- quality.data %>% mutate(quality=ifelse(is.na(quality),0.001,quality))
      }
      if(any(quality.data$quality<=0)) {
        cat("Some quality values are lower than 0. Rewriting them to 0.001.\n")
        quality.data <- quality.data %>% mutate(quality=ifelse(quality<=0,0.001,quality))
      }
      if(!"replicate" %in% names(quality.data)) {
        quality.data <- quality.data %>% mutate(replicate="1")
      }
      if(!"condition" %in% names(quality.data)) {
        quality.data <- quality.data %>% mutate(condition="1")
      }
    }
  }
  
  if(length(setdiff(replicates, unique(quality.data$replicate)))>0 | 
     length(setdiff(conditions, unique(quality.data$condition)))>0 |
     length(setdiff(unique(quality.data$condition), conditions))>0 |
     length(setdiff(unique(quality.data$replicate), replicates))>0
  ) {
    stop("The condition and replicate columns in the data and the quality.data must contain the same values.")
  }
  
  zstats <- list()
  for(cond in conditions) {
    zstats[[as.character(cond)]] <- list()
    
    for(rep in replicates) {
      cat(paste0("WORKING ON CONDITION ",cond,", REPLICATE ",rep,".\n"))
      
      cur.data <- data %>%
        filter(replicate==rep & condition==cond) %>%
        dplyr::select(protein,starts_with("PC"))
      
      cur.quality.data <- quality.data %>%
        filter(replicate==rep & condition==cond) %>%
        dplyr::select(protein,quality)
      
      # Calculate complex Euclidean distances
      cat("Calculating Euclidean distances for protein complexes...")
      complexEDs <- list()
      for(cpx in names(complexes)) {
        cpx.data <- cur.data %>%
          filter(protein %in% complexes[[cpx]])
        if(nrow(cpx.data)<2) {
          complexEDs[[cpx]] <- data.frame(log2mean.complexED=NA,n=NA,complex=cpx)
          next
        }
        cpx.quality.data <- cur.quality.data %>%
          filter(protein %in% complexes[[cpx]])
        
        complexEDs[[cpx]] <- data.frame(
          log2mean.complexED=log2(X.calculate.ED(cpx.data,cpx.quality.data,weights=weights,PCs=noPCs,cluster.cut=cluster.cut)$average.edist),
          n=nrow(cpx.data),
          complex=cpx
        )
      }
      complexEDs <- purrr::reduce(complexEDs,bind_rows) %>%
        na.omit()
      if(nrow(complexEDs)<1) {
        cat("\rCalculating Euclidean distances for protein complexes... undetected, skipped\n")
        next
      }
      cat("\rCalculating Euclidean distances for protein complexes...done\n")
      
      
      # Calculate random Euclidean distances
      cat("Calculating random Euclidean distances...\n")
      randomEDs <- X.random.ED(cur.data,cur.quality.data,members=sort(unique(complexEDs$n)),trials=trials,cluster.cut=cluster.cut,weights=weights,PCs=noPCs)
      
      zstats[[as.character(cond)]][[as.character(rep)]] <- complexEDs %>%
        left_join(randomEDs$summary,by="n") %>%
        mutate(nosd=abs(log2mean.complexED-log2mean.randomED)/log2sd.randomED)
    }
  }
  cat("Calculating assembly indeces...")
  zstats.table <- lapply(zstats, function(x) data.table::rbindlist(x,idcol ="replicate")) %>%
    data.table::rbindlist(idcol="condition") %>%
    rename(AI=nosd) %>%
    group_by(complex,condition) %>%
    mutate(mean.AI=mean(AI)) %>%
    ungroup() %>%
    group_by(complex) %>%
    mutate(rel.n=n/max(n)) %>%
    ungroup()
  zstats.plots <- list()
  plot.complexes <- sort(unique(zstats.table$complex))
  for(pc in plot.complexes) {
    if(y.trans) {
      zstats.plots[[pc]] <- zstats.table %>%
        mutate(condition=str_remove(condition,"tp")) %>%
        mutate(condition=factor(condition, levels=gtools::mixedsort(unique(condition)))) %>%
        filter(complex==pc) %>%
        mutate(across(ends_with("AF"), log2)) %>%
        ggplot(aes(x=condition,y=AF)) +
        geom_line(inherit.aes=FALSE, aes(x=condition,y=mean.AF,group=complex)) +
        geom_point(aes(color=replicate, size=rel.n*6)) +
        scale_y_continuous(name="Assembly index", limits=c(min(c(-5,min(zstats.table$AF))),max(5,max(zstats.table$AF)))) +
        scale_x_discrete(name="Timepoint (hpi)",drop=FALSE) +
        scale_size_continuous(range=c(0.5,4), limits=c(0,6), guide="none") +
        scale_color_manual(values=c("red3","blue2","green3")) +
        customPlot +
        guides(color = guide_legend(nrow = 1)) +
        theme(legend.position="bottom",
              axis.text=element_text(size=7))
    } else {
      zstats.plots[[pc]] <- zstats.table %>%
        mutate(condition=str_remove(condition,"tp")) %>%
        mutate(condition=factor(condition, levels=gtools::mixedsort(unique(condition)))) %>%
        filter(complex==pc) %>%
        ggplot(aes(x=condition,y=AF)) +
        geom_line(inherit.aes=FALSE, aes(x=condition,y=mean.AF,group=complex)) +
        geom_point(aes(color=replicate, size=rel.n*6)) +
        scale_y_continuous(name="Assembly index", limits=c(0,16)) +
        scale_x_discrete(name="Timepoint (hpi)",drop=FALSE) +
        scale_size_continuous(range=c(0.5,4), limits=c(0,6), guide="none") +
        scale_color_manual(values=c("red3","blue2","green3")) +
        customPlot +
        guides(color = guide_legend(nrow = 1)) +
        theme(legend.position="bottom",
              axis.text=element_text(size=7))
    }
  }
  cat("\rCalculating assembly indeces...\n")
  output <- list(
    data=zstats.table,
    plots=zstats.plots
  )
  return(output)
  
  
}
