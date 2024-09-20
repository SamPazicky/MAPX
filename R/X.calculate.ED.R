X.calculate.ED <- function(
  data=NULL,
  quality.data=NULL,
  centroid=TRUE,
  weights=NULL,
  PCs=0,
  cluster.cut=0.65
) {
  
  if(PCs==0) {
    PCs=data%>%dplyr::select(starts_with("PC")) %>% names() %>% length()
  }
  if(is.null(weights)) {
    weights=rep(1,PCs)
  }
  
  edata <- data %>%  
    column_to_rownames("protein") %>%
    dplyr::select(starts_with("PC")) %>%
    dplyr::select(order(names(.)))
  equality.data <- quality.data %>%
    column_to_rownames("protein") %>%
    .[rownames(edata),]
  
  if(centroid) {
    hc <- hclust(dist(edata))
    hc.height <- max(hc$height, na.rm=TRUE)
    hc.cut <- cutree(hc, h=hc.height*cluster.cut)
    biggest_cluster <- table(hc.cut) %>% sort(decreasing=TRUE) %>% .[1] %>% names()
    
    center <- data %>%
      # add quality data
      left_join(quality.data,by="protein") %>%
      # add the clustering results to each protein (hc = which cluster?)
      mutate(hc=cutree(hc,h=hc.height*cluster.cut)) %>%
      # calculate mean PC coordinates - this is the centroid of the complex
      # filter for the biggest cluster available for the particular condition
      filter(hc==biggest_cluster) %>%
      summarise(across(starts_with("PC"), function(x) weighted.mean(x, quality, na.rm=TRUE)))
   
      
    edists <- vector()
    for(m in 1:nrow(edata)) {
      edists[m] <- e.dist(edata[m,][1:PCs],center[1:PCs],weights[1:PCs])
    }
    av.edist <- weighted.mean(edists,equality.data)
    edists=data.frame(
      protein=rownames(edata),
      edist=edists
    )
    pp.table <- NULL
    
  } else {
    pptable <- expand.grid.unique(data$protein, data$protein) %>% as.data.frame() %>% setNames(paste0("protein",1:2))
    pp.edists <- vector()
    for(i in 1:nrow(pptable)) {
      pp.edists[i] <- e.dist(edata[pptable$protein1[i],], edata[pptable$protein2[i],],weights[1:PCs])
    }
    pptable <- pptable %>% mutate(pp.edist=pp.edists)
    edists <- bind_rows(
      pptable%>%dplyr::select(protein1,pp.edist) %>% setNames(c("protein","pp.edist")),
      pptable%>%dplyr::select(protein2,pp.edist) %>% setNames(c("protein","pp.edist")),
    ) %>%
      group_by(protein) %>%
      dplyr::summarise(edist=mean(pp.edist,na.rm=TRUE)) %>%
      ungroup()
    av.edist <- edists %>% left_join(quality.data, by="protein") %>% dplyr::summarise(weighted.mean(edist,quality)) %>% unlist() %>% unname()
  }
  
  output <- list(
    average.edist=av.edist,
    edists=edists,
    pp.table=pptable
  )
  
  return(output)
}
