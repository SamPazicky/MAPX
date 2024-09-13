X.random.ED <- function(
  data=NULL,
  quality.data=NULL,
  trials=1000,
  members=2:100,
  cluster.cut=0.65,
  weights=NULL,
  PCs=0,
  stats=TRUE
) {
  
  dists <- summaries <- vector("list",length=length(members)) %>% setNames(as.character(members))
 
  for(mno in members) {
    cat(paste0("Random Euclidean distances for complexes with ",mno," subunits:\n"))
    dists[[as.character(mno)]] <- vector("integer",length=trials)
    
    pb <- txtProgressBar(min=0, max=trials, style=3, initial="") # progress bar
    for(trial in 1:trials) {
      tdata <- data %>%
        slice_sample(n=mno) %>%
        arrange(protein)
      tquality.data <- quality.data %>%
        filter(protein %in% tdata$protein) %>%
        arrange(protein)
      
      dists[[as.character(mno)]][trial] <- (X.calculate.ED(tdata,tquality.data,weights=weights,PCs=PCs,cluster.cut=cluster.cut))$average.edist

      setTxtProgressBar(pb,trial)
    }
    mean.dists <- mean(log2(dists[[as.character(mno)]]))
    sd.dists <- sd(log2(dists[[as.character(mno)]]))
    summaries[[as.character(mno)]] <- data.frame(n=mno,log2mean.randomED=mean.dists,log2sd.randomED=sd.dists)
    
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
