#' X.predict.moonlighters
#'
#' Calculate assembly factors from complex features reduce into a PCA plot.
#' 
#' @param data Data frame with features reduce by PCA. Columns are 'protein', and PC columns starting with 'PC'. Additional columns
#' 'condition' and 'replicate' can be used.
#' @param complexes Named list: A list where each element is vector with proteins and the name of that element is the name of the protein complex.
#' @param plot Logical: Should the plots be plotted?
#' @param sig.cutoff Numeric: Significance level cutoff, 0.05 by default.
#' @param moon.proportion.cutoff Upper proportion of moonlighting subunits cutoff. Default is 1/3.
#' @param min.reps Integer: At least in how many replicates should a subunit be significantly moonlighting to appear in hits?
#' @param min.conditions Integer: At least in how many conditions should a subunit be significantly moonlighting to appear in hits?
#' @param noPCs Integer: How many principal components should be used for the clustering?
#' @param weights Numeric vector: weights for each principal component for k.means clustering.
#' the complex centroid, distance from which is further calculated and gives
#' @param kmeans.maxiter Integer: Maximum interation for kmeans clustering.
#' @param kmeans.nstart Integer: How many starting sets should be tried for kmeans clustering?
#' @param p.adj.method Character string: One of the p.adj methods, see ?p.adjust. Default is "none".
#' 
#' @import tidyverse
#' @importFrom gtools mixedsort
#' @import ggforce
#' 
#' @return PCA data and different supporting information.
#' @examples 
#' data.moonlighters <- X.predict.moonlighters(data=data.PCAs$data, min.reps=2)
#' @export
#'
#'

X.predict.moonlighters <- function(
    data=NULL,
    complexes=list(),
    plot=TRUE,
    sig.cutoff=0.05,
    moon.proportion.cutoff=1/3,
    min.reps=1,
    min.conditions=1,
    noPCs=0,
    weights=NULL,
    kmeans.maxiter=10,
    kmeans.nstart=5,
    p.adj.method="none"
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
    noPCs <- ifelse(noPCs==0,data%>%dplyr::select(starts_with("PC"))%>%names()%>%length(),noPCs)
    data <- data %>% dplyr::select(protein,condition,replicate,all_of(paste0("PC",1:noPCs)))
  }
  
  if(length(complexes)<1) {
    stop("Please include a named list with complex subunits")
  } else {
    if(length(names(complexes))<1) {
      stop("The list of complexes must be named.")
    }
    complexes <- lapply(complexes,unique)
  }
  
  if(is.null(weights)) {
    weights=rep(1,noPCs)
  } else if(length(weights)<lengths(noPCs)) {
    stop("weights must be a numeric vector of the same length as the number of principal components used.")
  }
  
  cat("Predicting moonlighters...")
  zdata <- data %>%
    mutate(across(starts_with("PC"), scale)) %>%
    dplyr::select(starts_with("PC")) %>%
    purrr::map2_dfc(weights, `*`)
  zdata <- bind_cols(data%>%dplyr::select(!starts_with("PC")), zdata)
  
  looptable <- expand.grid(names(complexes),conditions,replicates) %>%
    as.data.frame() %>% setNames(c("complex","condition","replicate"))
  complextable <- stack(complexes) %>% setNames(c("protein","complex"))
  result.data <- list()
  for(i in 1:nrow(looptable)) {
    all.clustdata <- looptable[i,] %>%
      left_join(complextable,by="complex") %>%
      left_join(zdata,by=c("protein","condition","replicate")) %>%
      na.omit() %>%
      remove_rownames()
    clustdata <- all.clustdata %>%
      column_to_rownames("protein") %>%
      dplyr::select(starts_with("PC"))
    if(nrow(clustdata)<3) {
      next
    }
    km <- kmeans(clustdata,centers=2,iter.max=kmeans.maxiter,nstart=kmeans.nstart)
    result.data[[i]] <- stack(km$cluster) %>% 
      setNames(c("cluster","protein")) %>%
      left_join(all.clustdata,by="protein") %>%
      mutate(significance=km$tot.withinss/km$totss)
  }
  
  result.table <- purrr::reduce(result.data,bind_rows) %>%
    group_by(complex,condition,replicate) %>%
    mutate(no.subunits.total=n()) %>%
    add_count(cluster,name="no.subunits.in.cluster") %>%
    mutate(cluster=ifelse(no.subunits.in.cluster==max(no.subunits.in.cluster),1,2)) %>%
    mutate(weighted.significance=significance/no.subunits.total) %>%
    ungroup() %>%
    mutate(proportion.of.moonlighters=no.subunits.in.cluster/no.subunits.total) %>%
    mutate(proportion.of.moonlighters=ifelse(cluster==1,NA,proportion.of.moonlighters))
    
  moonlighter.data <- result.table %>%
    na.omit() %>%
    filter(proportion.of.moonlighters<=moon.proportion.cutoff) %>%
    group_by(complex,protein,condition) %>%
    dplyr::summarise(sig=prod(significance),w.sig=prod(weighted.significance),replicates=paste(replicate,collapse=","), no.replicates=n(), .groups="keep") %>%
    ungroup() %>%
    mutate(w.sig.adj=p.adjust(w.sig,method=p.adj.method)) %>%
    filter(no.replicates>=min.reps) %>%
    group_by(complex,protein) %>%
    mutate(no.conditions.sig=n()) %>%
    ungroup() %>%
    filter(no.conditions.sig>=min.conditions) %>%
    filter(w.sig.adj<=sig.cutoff) %>%
    arrange(complex,protein,condition)
  
  cat("\rPredicting moonlighters... done.\n")
  
  if(plot) {
    cat("Plotting...")
    
    moonlighter.plot <- moonlighter.data %>%
      mutate(condition=factor(condition,levels=gtools::mixedsort(unique(moonlighter.data$condition)))) %>%
      ggplot(aes(y=protein,x=condition)) +
      ggforce::facet_col(facets = vars(complex), 
                         scales = "free_y", 
                         space = "free",
                         strip.position="top") +
      geom_line(aes(group=protein),size=0.1,color="gray50") +
      geom_point(aes(color=w.sig.adj)) +
      scale_color_gradientn(colours=c("red4","pink"),limits=c(0,0.05)) +
      customPlot +
      theme_bw() +
      theme(axis.text.y=element_text(size=7,color="black"),
            panel.spacing = unit(0.1, "lines"),
            strip.text=element_text(size=8,hjust=0,margin=ggplot2::margin(0,0,0,0,"pt")),
            axis.title=element_blank(),
            strip.background=element_rect(color="transparent",fill="transparent"),
            legend.title=element_blank(),
            legend.position="bottom")
    
    PCA.plots <- vector("list",length=length(unique(looptable$complex)))%>%setNames(names(complexes)) %>%
      lapply(function(x) x <- vector("list",length(conditions))%>%setNames(conditions) %>%
               lapply(function(y) y <- vector("list",length(replicates))%>%setNames(paste0("rep",replicates))
               )
      )
    for(i in 1:nrow(looptable)) {
      PCA.plots[[looptable[i,"complex"]]][[looptable[i,"condition"]]][[paste0("rep",looptable[i,"replicate"])]] <- looptable[i,] %>%
        left_join(complextable, by="complex") %>%
        left_join(data,by=c("protein","condition","replicate")) %>%
        left_join(moonlighter.data,by=c("complex","protein","condition")) %>%
        mutate(moonlighting=ifelse(is.na(no.conditions.sig),"No","Yes")) %>%
        ggplot(aes(x=PC1,y=PC2,color=moonlighting)) +
        geom_point() +
        scale_color_manual(values=c("black","red")) +
        theme_bw() +
        scale_x_continuous(limits=c(-2.5,2.5)) + scale_y_continuous(limits=c(-2.5,2.5))
    }
    cat("\rPlotting... done.\n")
  } else {
    moonlighter.plot=ggplot()
    PCA.plots=NA
  }
  
  
  output <- list(
    data=result.table,
    moonlighters=moonlighter.data,
    plot=moonlighter.plot,
    PCA.plots=PCA.plots
  )
  return(output)
}