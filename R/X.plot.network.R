#' X.plot.network
#'
#' Using one of the network building algorithms, cluster initial network into individual protein complexes.
#' 
#' @param data Data frame with pair-wise interactions (columns 'protein1' and 'protein2') and a column with probabilities from the ML model.
#' @param scores.col Character string: Name of the columns with prediction values.
#' @param standard.set Data frame with columns protein1, protein2 and another column with labels.
#' @param labels.col Character string: Name of the columns with labels.
#' @param annotation Data frame with column 'id' with protein IDs and any other annotation columns.
#' @param design.params List of lines of pieces of code to set up the design of the network plot.
#' @param export.destination Folder where the files will be saved. Default is the current directory.
#' @param network.name Title of the network in Cytoscape.
#' @param network.collection Collection of the network in Cytoscape.
#' @param cytoscape.path Path to Cytoscape executable file. If not given, the program will try to find it, which may take a long time.
#' @param cytoscape.waittime Integer: How long should the program wait for the cytoscape to initiate.
#' @import tidyverse
#' @import RCy3
#' 
#' @return Nothing to return, the network is saved in the export.destination in various formats.
#' @examples 
#' design=list('setNodeShapeDefault("ELLIPSE")','setNodeColorDefault("#5A5A5A")','lockNodeDimensions(TRUE)','setNodeLabelOpacityDefault(1)','setEdgeColorMapping(table.column="complex",list("1","0"), list("#00FF00","#FF0000"),mapping.type="d")','setEdgeLineWidthMapping(mapping.type="c",table.column="pred", widths=c(1,20))')
#' X.plot.network(final_data,scores.col="score",standard.set=GS,labels.col="complex",annotation=plasmoDB_data%>%rename(protein=Accession), 
#'               design.params=design, network.name="MAPX_network", cytoscape.path="C://Program Files/Cytoscape_v3.9.1/Cytoscape.exe")
#' @export
#' 
X.plot.network=function (
  data=NULL, # data
  scores.col=NA,
  standard.set=NULL,
  labels.col=NA,
  annotation=NULL,
  design.params=list(),
  export.destination=".",
  network.name="PPIN",
  network.collection="PPINs",
  cytoscape.path=NA,
  cytoscape.waittime=40
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
    if(!"protein1" %in% names(data) | !"protein2" %in% names(data)) {
      stop("Data must contain columns protein1 and protein2.")
    }
    data <- data %>%
      rename(source=protein1) %>% rename(target=protein2)
  }
  
  if(!is.null(standard.set)) {
    standard.set <- as.data.frame(standard.set)
    if(is.na(labels.col)) {
      labels.col <- names(standard.set)[ncol(standard.set)]
      cat("No labels column selected. Selecting the last column:",labels.col,"\n")
    }
    if(!"protein1" %in% names(standard.set) | !"protein2" %in% names(standard.set)) {
      stop("Standard set must contain columns protein1 and protein2.")
    }
    standard.set <- standard.set %>%
      rename(source=protein1) %>% rename(target=protein2)
  }

  if(!"protein" %in% names(annotation)) {
    stop("annotation must contain a column 'protein'")
  }
  
  #check whether cytoscape is on 
  
  is.on <- system2( 'tasklist' , stdout = TRUE) %>%
    sapply(function(x) str_extract(x,"^.*exe")) %>%
    unname() %>% na.omit()  %>% unique() %>% sort() %>%
    str_detect("Cytoscape") %>% any()
  
  if(is.na(cytoscape.path)) {
    cat("Cytoscape path not given. Searching for cytoscape. This might take some time.\n")
    cytoscape.path <- list.files(path="C:/",pattern="Cytoscape.exe",recursive=TRUE,full.names=TRUE) %>% str_subset(pattern="Program Files") %>% .[1]
  }
  
  # start Cytoscape
  if(!is.on) {
    cat("Starting Cytoscape...\n")
    system2(cytoscape.path,stdout=FALSE, wait=FALSE, minimized=TRUE)
    Sys.sleep(cytoscape.waittime) # wait 40 second to launch Cytoscape
    
    #check whether R communicates with Cytoscape
    cytoscapePing()
  }
  
  # define data - this is a table of pairwise interactions, where the two protein columns are rename to "source" and "target"
  # an additional column called "interaction" needs to be added, in this case the value should be "interacts" in all cases
  data <- data %>%
    mutate(interaction="interacts")
  
  if(!is.null(standard.set)) {
    data <- data %>% cross_join(standard.set,vars=c("source","target"),mode="left")
  }
  
  # define nodes - this is just a table of proteins that are used in data as well. The Accession codes need to be called "id".
  # in the first line, I use the list of data to define which proteins should be listed in nodes
  nodes <- data.frame(id=unique(c(data%>%pull(source),data%>%pull(target)))) %>%
    left_join(annotation%>%rename(id=protein))
  
  # print(paste0("Building protein network for Plasmodium falciparum: timepoint ", tp, ", replicate ", rp))
  # the following code creates the network
  createNetworkFromDataFrames(nodes,data, title=network.name, collection=network.collection)
  
  # NETWORK DESIGN
  
  if(length(design.params)!=0) {
    design.code <- design.params %>% paste(collapse="\n")
    eval(parse(text=design.code))
  }
  
  #The following lines of code save the network into png file:
  initial_network_png_file_name <- file.path(getwd(),"230_Isserlin_RCy3_intro", "images","initial_example_network.png")
  if(file.exists(initial_network_png_file_name)){
    #cytoscape hangs waiting for user response if file already exists.  Remove it first
    file.remove(initial_network_png_file_name)
  }
  
  #export the network as image
  exportImage(paste0(export.destination,"/",network.name,".png"), type = "png")
  exportImage(paste0(export.destination,"/",network.name,".svg"), type = "svg")
  
  #export the network as sif file
  exportNetwork(filename=paste0(export.destination,"/",network.name,".sif"))
  
  #save session
  saveSession(filename = paste0(export.destination,"/",network.name,".cys"))
  closeSession(FALSE)
  system2("taskkill", args = "/im CYTOSC~1.EXE")
  Sys.sleep(15)
}