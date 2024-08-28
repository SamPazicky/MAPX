#' X.plot.meltcurves
#'
#' Plots fitted melting curves.
#' @param data Data frame or list of data frames: Scaled data with removed outliers and fitting statistics, ideally $data output from X.scale.meltcurves.
#' @param fits List containing all fit objects, ideally $fits output from X.meltcurve.fit. The list must be named and the names must correspond to id column
#' values of the data. Alternatively, a list of such lists of same length as the lenght of the list of data.
#' @param samples Character vector: names of samples in the same order as the list of data. NA if only one data frame.
#' @param replicates Character vector: names of replicates in the same order as the list of data. NA if only one data frame.
#' @param color.col Character string: Which column of the data frame(s) should be used for color scaling?
#' @param pointshape.col Character string: Which column of the data frame(s) should be used for point shape definition?
#' @param colors Character vector: "rainbow" will plot curves from red to violet across sample. "distinct" will use preset distinct color scheme.
#' User can also define a vector of colors.
#' @param points Integer vector: numbers that define point types used for plotting.
#' @param fit.length Integer: How many points should be used for fitting curves. Default is 100 which is sufficient for plotting.
#' @param pdf.export Logical: If TRUE (default), a pdf with all plots will be exported.
#' @param pdf.folder Character string: Name of the directory for pdf export. Default is the working directory.
#' @param pdf.name Character string: Name of the exported pdf file. Default is 'MCcurves'
#'
#' @import tidyverse
#' @import patchwork
#' 
#' @return A list of melting curve plots.
#' @examples 
#' MCdata.plotted <- X.plot.meltcurves(all.sdata %>% unlist(recursive=FALSE),all.fits %>% unlist(recursive=FALSE), samples=rep(c("tp28","tp40"),each=3),replicates=c(1,2,3,1,2,3))
#' @export


X.plot.meltcurves <- function(
    data=NULL,
    fits=NULL,
    samples=NA, # vector of samples in the same order as list of data. NA if only one data frame.
    replicates=NA, # vector of replicates in the same order of list of data.NA if only one data frame.
    color.col="Sample", # or a column in the data frames
    pointshape.col="Replicate", # or a column in the data frames
    colors="rainbow", # rainbow, discrete or vector of custom colors
    points=c(16,17,18,8,4,3,9,10,11,12,13,14),
    fit.length=100,
    pdf.export=TRUE,
    pdf.folder=".",
    pdf.name="MCcurves"
) 
{
  
  
  get_legend <- function(myggplot) { # to extract legend from a single ggplot
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  
  if(is.null(data)|is.null(fits)) {
    stop("Please include the data (data frame or list of data frames) and fits (list of curve fit objects")
  }
  if(!is.null(nrow(data))) {
    data <- list(data)
    fits <- list(fits)
  }
  if(is.na(replicates[1])) {
    if(!is.na(samples[1])) {
      replicates=rep(NA,length(samples))
    } 
  } else {
    replicates <- as.character(replicates)
  }
  if(is.na(samples[1])) {
    if(!is.na(replicates[1])) {
      samples=seq(NA,length(replicates))
    } 
  } else {
    samples <- as.character(samples)
  }
  
  if(length(samples)!=length(data)|length(replicates)!=length(data)|length(fits)!=length(data)) {
    stop("There must be the same number of samples, replicates, data frames and lists of fits.")
  }
  
  # prepare data for plotting
  ldata <- list()
  for(i in seq_along(data)) {
    ldata[[i]] <- data[[i]] %>%
      pivot_longer(cols=starts_with("T"), names_to="Temperature",values_to="Fraction_soluble") %>%
      mutate(Temperature=as.numeric(str_remove(Temperature,"T"))) %>%
      mutate(Sample=as.character(samples[i])) %>%
      mutate(Replicate=as.character(replicates[i]))
  }
  ldata <- Reduce(bind_rows,ldata)
 
  # read temperatures
  temperatures <- ldata$Temperature %>% unique() %>% sort()
  
  # prepare linedata
  fakedata=data.frame(x=seq(min(temperatures),max(temperatures), length.out=fit.length))
  ilinedata <- list()
  for(i in seq_along(fits)) {
    curfits <- fits[[i]]
    plinedata <- list()
    cat(paste0("Calculating line curves, round ",i,"/",length(fits),"...\n"))
    pb <- txtProgressBar(min=0, max=length(curfits), style=3, initial="") # progress bar
    for(p in names(curfits)) {
      pfits <- curfits[[p]]
      plinedata[[p]] <- data.frame(y=predict(curfits[[p]],newdata=fakedata)) %>%
        bind_cols(fakedata) %>%
        setNames(c("Fraction_soluble","Temperature")) %>%
        mutate(protein=p) %>%
        mutate(Sample=as.character(samples[i])) %>%
        mutate(Replicate=as.character(replicates[i]))
      setTxtProgressBar(pb, which(names(curfits)==p))
    }
    close(pb)
    ilinedata[[i]] <- Reduce(bind_rows,plinedata)
  }
  linedata <- Reduce(bind_rows,ilinedata)
  cat("Lines calculated.\n")
  
  # set up plot designs
  if(length(colors)==1) {
    if(colors=="discrete") {
      plot.colors <- c("black","#E31A1C","green4","#6A3D9A","#FF7F00","dodgerblue2","gold1","skyblue2","#FB9A99","palegreen2","#CAB2D6","#FDBF6F","gray70",
                       "khaki2","maroon","orchid1","deeppink1","blue1","steelblue4","darkturquoise","green1","yellow4","yellow3","darkorange4","brown"
      )[1:length(unique(samples))] %>% as.character()
    }
    if(colors=="rainbow") {
      plot.colors <- rainbow(length(unique(samples)))
      if(length(plot.colors)==1) {
        plot.colors="black"
      }
    }
  } else {
    plot.colors <- colors
  }
  names(plot.colors) <- unique(samples)
  
  plot.shapes <- points[1:length(unique(replicates))]
  names(plot.shapes) <- unique(replicates)
  
  # plotting
  proteins <- sort(unique(linedata$protein))
  MCplots <- list()
  
  cat("Drawing plots...\n")
  pb <- txtProgressBar(min=0, max=length(proteins), style=3, initial="") # progress bar
  for(p in proteins) {
    
    ldata_p <- ldata %>% filter(protein==p)
    linedata_p <- linedata %>% filter(protein==p)
    
    MCplots[[p]] <- 
      ggplot(mapping=aes(Temperature,Fraction_soluble, color=!!sym(color.col), shape=!!sym(pointshape.col))) +
      ggtitle(p) +
      geom_point(data=ldata_p) +
      geom_line(data=linedata_p) +
      scale_shape_manual(values=plot.shapes, na.value=16, drop=FALSE) +
      scale_color_manual(values=plot.colors, na.value="black", drop=FALSE) +
      scale_y_continuous(limits=c(0,1), name="Fraction soluble") +
      scale_x_continuous(limits=c(min(temperatures),max(temperatures))) +
      theme_bw(base_size = 12) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.margin=ggplot2::margin(5,5,5,5, "pt"),
            plot.title = element_text(hjust=0.5, size=10), 
            legend.position="bottom",
            legend.title=element_text(size=6),
            legend.text=element_text(size=8),
            axis.text=element_text(size=8))
    
    setTxtProgressBar(pb, which(proteins==p))
  }
  close(pb)
  
  if(pdf.export) {
    cat("Preparing for pdf export...\n")
    # get a legend
    legenddata <- expand.grid(samples,replicates) %>% as.data.frame() %>% 
      setNames(c("Sample","Replicate")) %>%
      mutate(Temperature=1) %>% mutate(Fraction_soluble=1)
    legendplot <- ggplot(legenddata,aes(x=Temperature,y=Fraction_soluble,color=Sample,shape=Replicate)) +
      geom_point() +
      scale_shape_manual(values=plot.shapes, na.value=16, drop=FALSE) +
      scale_color_manual(values=plot.colors, na.value="black", drop=FALSE) +
      theme(legend.position="bottom")
    legend <- get_legend(legendplot)
    
    # list all plots
    pdfplots <- lapply(MCplots, function(x) x + theme(axis.text=element_text(size=6),
                                                      plot.title=element_text(hjust=0.5, size=10),
                                                      legend.position="none",
                                                      axis.title=element_blank(),
                                                      legend.spacing.y=unit(-2,"cm"),
                                                      legend.margin = ggplot2::margin(0,0,0,0, unit="cm"),
                                                      legend.background=element_rect(fill="transparent")) +
                         guides(color = guide_legend(nrow = 1),fill=guide_legend(nrow=1))
    )
    
    cat("Saving pdf file. This may take several minutes, depending on the number of detected proteins...\n")
    
    # patchwork assembly
    
    layout <- "
    ABBBBBBBBB
    #CCCCCCCCC
    #DDDDDDDDD
    #EEEEEEEEE
    #FFFFFFFFF"
    glob_lab <- "Fraction soluble"
    
    emptyplot <- ggplot() + theme_void() + theme(panel.border=element_blank())
    eplist=list()
    for (ep in 1:20) {
      eplist[[ep]]=emptyplot
    }
    y_lab <- 
      ggplot() + 
      annotate(geom = "text", x = 1, y = 1, label = glob_lab, angle = 90) +
      coord_cartesian(clip = "off")+
      theme_void()
    
    glob_lab <- "Temperature (°C)"
    x_lab <- 
      ggplot() + 
      annotate(geom = "text", x = 1, y = 1, label = glob_lab, angle = 0) +
      coord_cartesian(clip = "off")+
      theme_void()
    
    pdf.folder <- str_remove(pdf.folder,"/$")
    pdf(paste0(pdf.folder,"/",pdf.name,".pdf"), width = 10 , height = 14.5)
    no.pages=ceiling(length(pdfplots)/20)
    cat("Printing the plots into a pdf document...\n")
    pb <- txtProgressBar(min=0, max=no.pages, style=3, initial="") # progress bar
    for (page in 1:no.pages) {
      if(page==no.pages) {
        print(
          patchwork::wrap_elements(y_lab) +
            patchwork::wrap_plots(c(pdfplots[((20*page)-19):(length(pdfplots))],eplist[0:(20-length(((20*page)-19):(length(pdfplots))))]), ncol = 4, nrow=5 ) + 
            patchwork::plot_spacer() +
            patchwork::wrap_elements(x_lab) +
            patchwork::plot_spacer() +
            patchwork::wrap_elements(legend) +
            patchwork::plot_layout(guides = 'collect', design=layout, heights=c(1,-0.05,0.1,-0.08,0.1))
          
        )
      } else {
        print(
          patchwork::wrap_elements(y_lab) + 
            patchwork::wrap_plots(pdfplots[((20*page)-19):(20*page)], ncol = 4, nrow=5) + 
            patchwork::plot_spacer() +
            patchwork::wrap_elements(x_lab) +
            # patchwork::guide_area() + 
            patchwork::plot_spacer() +
            patchwork::wrap_elements(legend) +
            patchwork::plot_layout(guides = 'collect', design=layout, heights=c(1,-0.05,0.1,-0.08,0.1))
          
        )
      }
      setTxtProgressBar(pb, page)
    }
    close(pb)
    dev.off()
  }
  
  return(MCplots)
}
