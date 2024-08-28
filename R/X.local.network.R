#' X.local.network
#'
#' Calculate and plot a local subnetwork of a specific protein or proteins of interest.
#' 
#' @param data List of data frame with columns for protein1, protein2 and a column with scores
#' @param proteins Vector: proteins whose local network should be plotted.
#' @param cutoff Numeric: Probabilities cutoff applied to the data.
#' @param scores.col Character string: Name of the columns with prediction values.
#' @param condition.col Character string: If the data has more conditions, what is the name of the column containing conditions?
#' @param min.conditions Numeric: At least how many conditions must contain the predicted interaction to not be filtered out?
#' @param plot Character string: "cytoscape" for plotting using cytoscape and "visnetwork" for plotting using visNetwork package. NULL for no plotting.
#' @param plot.annotation Data frame with column 'id' with protein IDs and any other annotation columns.
#' @param plot.design.params List: as design.params in X.plot.network.
#' @param plot.cytoscape.path String: As in X.plot.network.
#' @param plot.cytoscape.waittime Numeric: How many second should R wait for the cytoscape to start.
#' @param plot.name String: Plot name for saved cytoscape network.
#' @param plot.path String: Output destination for the plot.
#' 
#' @import tidyverse
#' @import visNetwork
#' 
#' @return A data frame with network edges.
#' @examples 
#' local.subnetwork <- X.local.network(calibrated.model$data%>%dplyr::select(!score),"PF3D7_0412200",
#'             plot="cytoscape",
#'             plot.annotation=plasmoDB_data%>%rename(protein=Accession),
#'             plot.cytoscape.path="C://Program Files/Cytoscape_v3.9.1/Cytoscape.exe")
#' @export
#' 
#' 
X.local.network <- function(
    data=c(),
    proteins=c(),
    cutoff=0.5,
    scores.col=NA,
    condition.col=NA,
    min.conditions=1,
    plot=NULL,
    plot.annotation=NULL,
    plot.design.params=list(),
    plot.cytoscape.path=NA,
    plot.cytoscape.waittime=40,
    plot.name=paste(proteins,collapse="_"),
    plot.path="./"
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
    if(is.na(condition.col)) {
      data <- data %>%
        mutate(condition=1)
      min.conditions <- 1
    } else {
      data <- data %>%
        rename(condition=!!sym(condition.col))
    }
    wdata <- data %>%
      dplyr::select(protein1,protein2,score,condition)
  }
  
  pdata <- wdata %>%
    filter(protein1 %in% proteins | protein2 %in% proteins) %>%
    filter(score>=cutoff)
  
  int.proteins <- unique(c(pdata$protein1,pdata$protein2)) %>% gtools::mixedsort()
  if(length(int.proteins)<=length(proteins)) {
    message("No interacting proteins with proteins ", paste(proteins, collapse=", "))
    return(data.frame())
  }
  combs <- combn(int.proteins, 2, simplify = FALSE) %>% 
    lapply(as.data.frame) %>% lapply(t) %>% lapply(as.data.frame) %>% rbindlist() %>% 
    remove_rownames() %>% setNames(c("protein1","protein2"))
  
  ndata <- wdata %>%
    filter(protein1 %in% int.proteins & protein2 %in% int.proteins) %>%
    MAPX::cross_join(combs, vars=c("protein1","protein2"), mode="right") %>%
    filter(score>=cutoff) %>%
    group_by(protein1,protein2) %>%
    filter(n_distinct(condition) >= min.conditions) %>%
    ungroup() %>%
    arrange(protein1,protein2,condition) 
  
  int.proteins <- ndata %>%
    filter(protein1 %in% proteins | protein2 %in% proteins) %>%
    dplyr::select(protein1,protein2) %>% unlist() %>% unname()
  
  ndata <- ndata %>%
    filter(protein1 %in% int.proteins & protein2 %in% int.proteins)
  
  ndata_complexes <- MAPX::X.pairwise.to.complexes(ndata)
  which_complex <- list()
  for(protein in proteins) {
    which_complex[[protein]] <- ndata_complexes %>% sapply(function(x) protein %in% x) %>% which()
    if(length(which_complex)==0) {
      message("No interactions retrieved for protein ",protein)
    }
  }
  if(sum(sapply(which_complex,length))==0) {
    message("No complexes retrieved for all inputed proteins.")
    return(data.frame())
  }
  which_complex <- unique(which_complex) %>% unlist()
  
  ndata <- ndata %>%
    filter(protein1 %in% unlist(ndata_complexes[[which_complex]]) & protein2 %in% unlist(ndata_complexes[[which_complex]]))
  
  
  nproteins <- unique(c(ndata$protein1,ndata$protein2))
  
  if(is.null(plot)) {
    
  } else if (plot=="visnetwork") {
    
    if(!"protein" %in% names(plot.annotation)) {
      stop("annotation must contain a column 'protein'")
    }
    
    nodes <- data.frame(protein=nproteins) %>%
      left_join(plot.annotation) %>%
      rename(id=protein) %>%
      # left_join(plasmoDB%>%setNames(c("id","description"))) %>%
      mutate(shape="circle") %>%
      mutate(color=ifelse(id==protein,"red","gray")) %>%
      mutate(title = paste0("<b><u>",id,"</b></u><br>"))
    
    edges <- ndata %>%
      group_by(protein1,protein2) %>%
      summarise(scores=sum(score),conditions=paste0(condition,collapse=", "),n=n()) %>%
      ungroup() %>%
      mutate(width=(scores)) %>%
      # left_join(data.frame(condition=timepoints,color=my.rainbow)) %>%
      dplyr::select(starts_with("protein"),width,conditions) %>%
      setNames(c("from","to","width","conditions")) %>%
      mutate(length=100) %>%
      mutate(title=paste0("Found in timepoints: ", conditions))
    
    vis <- visNetwork(nodes,edges) %>%
      visIgraphLayout(layout = "layout_with_fr") %>%
      visNodes(scaling = list(label = list(enabled = T)))
    
    plot.path <- ifelse(endsWith(plot.path,"/"),plot.path,paste0(plot.path,"/"))

    visSave(vis,file=paste0(plot.path,"/",protein,".html"), selfcontained=TRUE)
    
  } else if (plot=="cytoscape") {
  
    MAPX::X.plot.network(data=ndata,scores.col="score",annotation=plot.annotation,design.params=plot.design.params,
                         export.destination=plot.path,cytoscape.path=plot.cytoscape.path,network.collection=protein, network.name=plot.name,
                         cytoscape.waittime=plot.cytoscape.waittime)
  }
  if(length(proteins)==1) {
    message("Complex retrieved for protein ", proteins,".")
  } else {
    message("Complex retrieved for proteins ", paste(proteins, collapse=", "),".")
  }
  return(ndata)
}
