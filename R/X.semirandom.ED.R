X.semirandom.ED <- function(
    data=NULL,
    quality.data=NULL,
    complex.data=NULL,
    trials=1000,
    weights=NULL,
    PCs=0,
    stats=TRUE
) {
  
  
  complex.data <- complex.data %>%
    group_by(complex) %>%
    mutate(n=n_distinct(protein)) %>%
    ungroup() %>%
    left_join(data, by="protein")
  
  dists <- summaries <- vector("list",length=nrow(complex.data))
  
  for(cpd in 1:nrow(complex.data)) {
    mno <- complex.data[cpd,] %>% pull(n)
    p <- complex.data[cpd,] %>% pull(protein)
    cpx <- complex.data[cpd,] %>% pull(complex)
    cat(paste0("Semirandom Euclidean distances for complexes with protein ",p," in complex ", cpx,". (", cpd,"/",nrow(complex.data),")\n"))
    
    dists[[cpd]] <- vector("list",length=trials)
    pb <- txtProgressBar(min=0, max=trials, style=3, initial="") # progress bar
    for(trial in 1:trials) {
      pdata <- data %>% filter(protein==p)
      tdata <- data %>% filter(protein!=p) %>% slice_sample(n=(mno-1)) %>% arrange(protein)
      tquality.data <- quality.data %>%
        filter(protein %in% tdata$protein) %>%
        arrange(protein)
      dists[[cpd]][[trial]] <- vector("integer",length=nrow(tdata))
      for(trow in 1:nrow(tdata)) {
        dists[[cpd]][[trial]][trow] <- e.dist(pdata%>%dplyr::select(!protein) %>% unlist %>% unname,
                                              tdata[trow,]%>%dplyr::select(!protein) %>% unlist %>% unname,
                                              weights)
      }
      dists[[cpd]][[trial]] <- weighted.mean(dists[[cpd]][[trial]], tquality.data$quality)
      setTxtProgressBar(pb,trial)
    }
    mean.dists <- mean(log2(unlist(dists[[cpd]])))
    sd.dists <- sd(log2(unlist(dists[[cpd]])))
    summaries[[cpd]] <- data.frame(complex=cpx,protein=p,n=mno,log2mean.semirandomED=mean.dists,log2sd.semirandomED=sd.dists)
    
    close(pb)
  }
  summary <- purrr::reduce(summaries,bind_rows)
  
  if(stats) {
    nstats <- list()
    for(no_subunits in names(dists)) {
      aemin <- dists[[no_subunits]] %>% min()
      ae5 <- dists[[no_subunits]] %>% quantile(0.05)
      ae10 <- dists[[no_subunits]] %>% quantile(0.10)
      aemax <- dists[[no_subunits]] %>% max()
      nstats[[as.character(no_subunits)]] <- data.frame(
        aemin=aemin, ae5=ae5, ae10=ae10,aemax=aemax
      ) %>% mutate(no_subunits=no_subunits)
    }
    stats_all <- purrr::reduce(nstats,bind_rows) %>% remove_rownames()
  } else {
    stats_all=NA
  }
  
  output <- list(data=dists,summary=summary,stats=stats_all)
  return(output)
}
