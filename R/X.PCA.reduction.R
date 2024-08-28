#' X.PCA.reduction
#'
#' Examine state of an assumed complex on a PCA plot using extracted features
#' 
#' @param data Data frame with features. Columns are 'protein', 'condition', 'replicate' and feature columns.
#' @param average Logical: Should the replicates be averaged before the PCA? Default is FALSE.
#' @param plot Logical: Should the results be plotted?
#' @param separate.plots Logical: Should the plots be plotted separately for different conditions? Default is FALSE.
#' @param plot.complexes Named list of vectors of complexes.
#' @param plot.colors Vector: Colors to be used for plotting
#' @param plot.colors.by Character string: By which column should the plots be colored? Default is 'condition'.
#' @param plot.path Logical: Should the same proteins by connected between the conditions?
#' @param plot.animate Character string: 'Time' for time-based animation and 'State' for state-based animation. Default is NA for no animation.
#' 
#' @import tidyverse
#' @import caret
#' 
#' @return List of five elements: $data contains PCA coordinates for each protein, $PCA contains information about the PCA,
#' var.explained shows what percentage of variation is explained by which PC, $plots are PCA plots and $separate.plots is
#' a list of separate PCA plots. 
#' @examples 
#' data.PCAs <- X.PCA.reduction(all.features.data,plot=TRUE,plot.complexes=complexes)
#' @export
#' 

X.PCA.reduction <- function(
    data=NULL,
    average=FALSE,
    plot=FALSE,
    separate.plots=TRUE,
    plot.complexes=NA,
    plot.colors=NA,
    plot.color.by="condition",
    plot.path=FALSE,
    plot.animate=NA #time or state
    
) {
  
  if(is.null(data)) {
    stop("Input data")
  } else {
    data <- as.data.frame(data)
    if(!"protein" %in% names(data)) {
      stop("The data must contain a column named 'protein'.")
    }
    if(!"condition" %in% names(data)) {
      data <- data %>%
        mutate(condition="1")
    }
    if(!"replicate" %in% names(data)) {
      data <- data %>%
        mutate(replicate="1")
    }
  } 
  
  if(any(is.na(plot.colors))) {
    now.rainbow <- rainbow(length(unique(data$condition)))
  } else {
    if(length(plot.colors)<length(unique(data$condition))) {
      stop("Colors must be as long as the number of conditions in the data.")
    }
    now.rainbow <- plot.colors
  }
  
  customPlot <- list(
    theme_bw(base_size = 12),
    # scale_fill_brewer(palette = "Set1"),
    # scale_colour_brewer(palette = "Set1"),
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())
  )
  
  PCA_alldata <- data.frame()
  
  if(average) {
    feed_data <- data %>%
      group_by(protein,condition) %>%
      dplyr::select(!replicate) %>%
      dplyr::summarise(across(where(is.numeric), function(x) mean(x, na.rm=TRUE)), .groups="keep") %>%
      ungroup() %>%
      mutate(replicate="averaged")
  } else {
    feed_data <- data
  }
  
  feed_data <- feed_data %>% mutate(sample=paste0(condition,";_;",replicate))
  cat("Calculating PCA...")
  # for each condition, the centering and scaling is done separately
  for(cond in unique(feed_data$sample)) {
    
    # modify data for PCA analysis
    PCA_partdata <- feed_data %>% 
      filter(sample==cond) %>%
      dplyr::select(!sample) %>%
      group_by(protein) %>%
      summarise(across(!c(replicate,condition), ~ mean(.x, na.rm=TRUE))) %>%
      ungroup() %>%
      mutate(across(where(is.numeric), ~ifelse(is.na(.x), min(.x,na.rm=TRUE), .x)))
    
    # center and scale
    PCA_preprocess <- preProcess(PCA_partdata,method=c("center","scale"), n.comp=2)
    
    PCA_partdata <- predict(PCA_preprocess, PCA_partdata) %>%
      mutate(protein=paste0(protein,";;",cond)) %>%
      column_to_rownames("protein")
    
    
    PCA_alldata <- bind_rows(PCA_alldata, PCA_partdata)
    
  }
  
  # Calculate PCA
  all_PCAs_together <- PCA_alldata %>% as.matrix() %>% prcomp
  
  # What percentage of variability is explained by which PCA?
  perc_var <-
    round(100 * all_PCAs_together$sdev ^ 2 /
            sum(all_PCAs_together$sdev ^ 2), 1)
  
  # PCA table
  all_PCAs_together_data <- all_PCAs_together %>%
    .$x %>%
    as.data.frame() %>% 
    rownames_to_column("protein;sample") %>%
    separate(col="protein;sample", into=c("protein","sample"), sep=";;") %>%
    separate(col="sample",remove=FALSE,into=c("condition","replicate"),sep=";_;")
  
  cat("\rCalculating PCA... done.\n")
  
  # plotting
  if(plot) {
    cat("Plotting...")
    
    PCA.plots <- list()
    separate_plots <- list()
    
    if(is.na(plot.complexes[1])) {
      stop("For plotting, you need to include a named list of vectors that represents plotted protein complexes.")
    }
    complexes <- plot.complexes
    
    if(length(complexes)==1) {
      if(is.na(complexes)) {
        stop("For plotting,include list of protein complexes")
      }
    }
    
    # loop through each complex
    for(cpx in seq_along(complexes)) {
      
      # isolate data for the complex in loop
      nowdata <- all_PCAs_together_data %>%
        mutate(sample=factor(sample, levels=gtools::mixedsort(unique(data$condition)))) %>%
        mutate(condition=factor(condition, levels=gtools::mixedsort(unique(data$condition)))) %>%
        filter(protein %in% complexes[[cpx]])
      
      if(nrow(nowdata)<2) {
        PCA.plots[[names(complexes)[cpx]]] <- ggplot() + ggtitle("Less than 2 subunit detected.") +
          theme_void()
        next
      }
      
      # choose colors
      if(plot.color.by=='condition') {
        now.colors <- now.rainbow[which(gtools::mixedsort(unique(data$condition)) %in% gtools::mixedsort(unique(nowdata$condition)) )]
        now.shapes=c(16,15,17,18,8,4,3,9,10,11,12,13,14)[1:length(nowdata$replicate)]
      } else {
        now.colors=rainbow(length(nowdata$protein))
        now.shapes= (c(16,15,17,18,8,4,3,9,10,11,12,13,14))[1:length(nowdata$protein)]
      }
      
      # just for animation - this needs to be finished XXX
      if(!is.na(plot.animate)) {
        if(plot.animate=="time") {
          nowdata <- nowdata %>%
            mutate(condition=as.integer(as.character(condition)))
          PCA.plots[[names(complexes)[cpx]]] <- nowdata %>%
            ggplot(aes(x=PC1,y=PC2)) +
            geom_point(data=all_PCAs_together_data%>%dplyr::select(!condition), size=0.5, alpha=0.2, color="gray50") +
            geom_point(color="black", size=3.5) +
            ggtitle(names(complexes[cpx])) +
            scale_x_continuous(name=paste0("PC1 (explains ",perc_var[1],"% variation)")) +
            scale_y_continuous(name=paste0("PC2 (explains ",perc_var[2],"% variation)")) +
            customPlot +
            theme(legend.position="none")
        } else if(plot.animate=="states") {
          
          PCA.plots[[names(complexes)[cpx]]] <- nowdata %>%
            ggplot(aes(x=PC1,y=PC2,color=condition,group=protein)) +
            geom_point(data=all_PCAs_together_data%>%dplyr::select(!condition), size=0.5, alpha=0.2, color="gray50") +
            geom_point(size=3.5) +
            ggtitle(names(complexes[cpx])) +
            scale_x_continuous(name=paste0("PC1"),limits=c(-4.5,4.5)) +
            scale_y_continuous(name=paste0("PC2"),limits=c(-4.5,4.5)) +
            scale_color_manual(values=now.colors) +
            customPlot +
            theme(legend.position="none")
          
        }
        
        
      } else {
        
        PCA.plots[[names(complexes)[cpx]]] <- nowdata %>%
          ggplot(aes(x=PC1,y=PC2)) +
          geom_point(data=all_PCAs_together_data, size=0.5, alpha=0.2, color="gray50")
        
        if(plot.path) {
          PCA.plots[[names(complexes)[cpx]]] <- PCA.plots[[names(complexes)[cpx]]] +
            geom_path(aes(group=protein), color="black",alpha=0.7)
        }
        
        PCA.plots[[names(complexes)[cpx]]] <- PCA.plots[[names(complexes)[cpx]]] +
          geom_point(aes(color=as.character(!!sym(plot.color.by)), shape=replicate), size=3.5) + ## this shape depends on color.by or I should add shape.by XXX
          ggtitle(names(complexes[cpx])) +
          scale_color_manual(values=now.colors) +
          scale_shape_manual(values=now.shapes) +
          scale_x_continuous(name=paste0("PC1 (explains ",perc_var[1],"% variation)")) +
          scale_y_continuous(name=paste0("PC2 (explains ",perc_var[2],"% variation)")) +
          customPlot +
          theme(legend.position="bottom")
        
        if(separate.plots) {
          cur_conds <- unique(nowdata$condition)
          cur_reps <- unique(nowdata$replicate)
          for(cc in cur_conds) {
            separate_plots[[paste0(names(complexes)[cpx],"_cond",cc)]] <- list()
            for(cr in cur_reps) {
              separate_data <- nowdata %>%
                filter(condition==cc) %>% filter(replicate==cr)
              if(nrow(separate_data)<2) {
                separate_plots[[paste0(names(complexes)[cpx],"_cond",cc)]][[paste0("rep_",cr)]] <- ggplot() + ggtitle("Less than 2 subunit detected.") +
                  theme_void()
                next
              } else {
                separate_plots[[paste0(names(complexes)[cpx],"_cond",cc)]][[paste0("rep_",cr)]] <-separate_data %>%
                  ggplot(aes(x=PC1,y=PC2,color=as.character(!!sym(plot.color.by)),group=protein)) +
                  geom_point(data=all_PCAs_together_data%>%dplyr::select(!condition), size=0.5, alpha=0.2, color="gray50") +
                  geom_point(size=3.5) +
                  ggtitle(names(complexes[cpx])) +
                  scale_x_continuous(name=paste0("PC1 (explains ",perc_var[1],"% variation)")) +
                  scale_y_continuous(name=paste0("PC2 (explains ",perc_var[2],"% variation)")) +
                  scale_color_manual(values=now.colors,drop=FALSE) +
                  customPlot +
                  theme(legend.position="none")
              }
            }
          }
        }
      }
    }
    cat("\rPlotting... done.\n")
  } else {
    PCA.plots=NULL
    separate_plots<-"Not plotted."
  }
  
  if(!separate.plots) {
    separate_plots<-"Not plotted."
  }
  
  output <- list(
    data=all_PCAs_together_data,
    PCA=all_PCAs_together,
    var.explained=perc_var,
    plots=PCA.plots,
    separate.plots=separate_plots
  )
  return(output)
}
