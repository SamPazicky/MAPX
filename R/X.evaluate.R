#' X.evaluate
#'
#' Evaluate classification prediction results.
#' 
#' @param data Data frame with columns for model scores and labels.
#' @param scores.col Character string: 'data' column name with probs that is to be evaluated.
#' @param labels.col Character string: 'data' column name with labels (1 for complex-forming, 0 for others.)
#' @param eval.metric Character string: Method for predictor evaluation. 'roc' for area under the receiver-operator curve, 'prc' for area under the
#' precision-recall curve, 'kappa' for Cohen's kappa and 'F1' for F1 score.
#' @param plot Logical: Should the plots be plotted and outputted?
#' @param black Logical: Should the plotted curve be black or color gradient? Default is FALSE.
#' @param kappa.cutoff Numeric: What is the score cutoff for evaluating based on Cohen's kappa? Default is 0.5.
#' @param kappa.weight Character string: Same as argument weight in irr:kappa2
#' @param F1.cutoff Numeric: What is the score cutoff for evaluating based on F1? Default is 0.5.
#' 
#' @import tidyverse
#' @import pROC
#' @import PRROC
#' @import irr
#' @import viridis
#' 
#' @return A list with three elements. $eval.metric is a number corresponding to the chosen evaluation method. $curvedata are data that allow plotting
#' a plot corresponding to the evaluation metric. $plot is the plot of the corresponding metric.
#' @examples 
#' evaluation <- X.evaluate(data)
#' @export
#' 
X.evaluate <- function(
    data=NULL, # data frame with data
    scores.col=NA, # column name with probs to be validated
    labels.col=NA,
    eval.metric="prc",
    plot=FALSE,
    black=FALSE,
    kappa.cutoff=0.5,
    kappa.weight="unweighted",
    F1.cutoff=0.5
) 

{
  
  # plot settings
  
  customPlot <- list(
    theme_bw(base_size = 12),
    # scale_fill_brewer(palette = "Set1"),
    # scale_colour_brewer(palette = "Set1"),
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())
  )
  
  if(is.null(data)) {
    stop("Include data")
  } else {
    data <- as.data.frame(data)
    if(is.na(labels.col)) {
      cat("Column labels not specified. Taking the last column.\n")
      data <- data %>%
        rename(labels=names(.)[ncol(.)])
    } else {
      data <- data %>%
        rename(labels=!!sym(labels.col))
    }
    if(is.na(scores.col)) {
      cat("Column with scores not specified. Taking the first column.\n")
      data <- data %>%
        rename(score=names(.)[1])
    } else {
      data <- data %>%
        rename(score=!!sym(scores.col))
    }
    data <- data %>% dplyr::select(score,labels)
    if(length(setdiff(unique(data$labels),c(0,1)))>0) {
      stop("The label must contain values 0 and 1 for non-interacting and interacting protein pairs, respectively.")
    }
  }
  
  if(length(unique(data$score))==1) {
    output <- list(eval.metric=NA,curvedata=NA,plot=NA)
    return(output)
  }

  plotting=NULL
  
  if(eval.metric=="roc") {
    eval <- pROC::roc(labels ~ as.numeric(score), data, levels=c(0,1),direction="<")
    metric <- c("Area under the receiver-operator curve"=as.numeric(eval[["auc"]]))
    curvedata <- data.frame(Sensitivity=eval$sensitivities,Specificity=eval$specificities,Probability=eval$thresholds)
    
    if(plot) {
      if(black) {
        plotting = ggplot(data=curvedata, aes(y=Sensitivity,x=Specificity)) +
          geom_line(linewidth=1.5, color="black") +
          scale_x_reverse(expand=c(0.001,0.001)) +
          scale_y_continuous(expand=c(0.001,0.001)) +
          geom_abline(intercept=1, slope=1) +
          customPlot +
          ggtitle(paste0("AUROC=",round(metric,2)))
      } else {
        plotting = ggplot(data=curvedata, aes(y=Sensitivity,x=Specificity, color=Probability)) +
          geom_line(linewidth=1.5) +
          scale_color_viridis(option="turbo") +
          scale_x_reverse(expand=c(0.001,0.001)) +
          scale_y_continuous(expand=c(0.001,0.001)) +
          geom_abline(intercept=1, slope=1) +
          customPlot +
          ggtitle(paste0("AUROC=",round(metric,2)))
      }
    }
    
  } else if (eval.metric=="prc") {
    
    eval <- pr.curve((data %>% filter(labels==1)%>%pull(score)%>%na.omit()),
                     (data %>% filter(labels==0)%>%pull(score)%>%na.omit()),
                     curve=TRUE)
   
    metric=c("Area under the precision-recall curve"=eval[[2]])
    curvedata <- eval$curve %>% as.data.frame() %>% setNames(c("Recall","Precision","Probability"))
    
    if(plot) {
      if(black) {
        plotting <- ggplot(data=curvedata, aes(x=Recall,y=Precision)) +
          geom_line(linewidth=1.5, color="black") +
          scale_x_continuous(expand=c(0.001,0.001), limits=c(0,1)) +
          scale_y_continuous(expand=c(0.001,0.001), limits=c(0,1)) +
          customPlot +
          ggtitle(paste0("AUPRC=",round(metric,2)))
      } else {
        plotting <- ggplot(data=curvedata, aes(x=Recall,y=Precision,color=Probability)) +
          geom_line(linewidth=1.5) +
          scale_color_viridis(option="turbo") +
          scale_x_continuous(expand=c(0.001,0.001),limits=c(0,1)) +
          scale_y_continuous(expand=c(0.001,0.001),limits=c(0,1)) +
          customPlot +
          ggtitle(paste0("AUPRC=",round(metric,2)))
      }
    }
  } else if(eval.metric=="kappa") {
    
    eval <- data %>% select(score,labels) %>% mutate(score=ifelse(score>=kappa.cutoff,1,0)) %>% 
      irr::kappa2(weight=kappa.weight)
    metric <- c("Cohen's kappa"=eval$value,"p-value"=eval$p.value)
    curvedata=NULL
    if(plot) {
      plotting=NULL
    }
  } else if(eval.metric=="F1") {
    TP <- data %>% filter(score>=F1.cutoff) %>% filter(labels==1) %>% nrow()
    FP <- data %>% filter(score>=F1.cutoff) %>% filter(labels==0) %>% nrow()
    FN <- data %>% filter(score<F1.cutoff) %>% filter(labels==1) %>% nrow()
    prec <- TP/(TP+FP)
    rec <- TP/(TP+FN)
    metric <- c("F1 score"=2/(1/prec+1/rec),"threshold"=F1.cutoff)
    curvedata=NULL
    if(plot) {
      plotting=NULL
    }
  }
  
  output <- list(eval.metric=metric,curvedata=curvedata,plot=plotting)
  return(output)
}

