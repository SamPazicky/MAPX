#' X.calibrate.scores
#'
#' Calibrate model scores with sigmoidal model to approximate interaction probabilities
#' 
#' @param data Data frame with columns for protein1, protein2 and a column with scores.
#' @param standard.set Data frame with column for protein1, protein2 and a column with labels.
#' @param scores.col Name of the data column with model scores.
#' @param labels.col Name of the standard.set column with labels with values 0 for non-interacting and 1 for interacting proteins.
#' 
#' @import tidyverse
#' 
#' @return A list of three: $data with the ata frame with the original data and an additional column 'probabilities' with the calibrated probabilities,
#' $plot with the calibration plot and $model with the model fitted into the original data.
#' 
#' @examples
#' calibrated.model <- X.calibrate.scores(final.model$data, GS)
#' 
#' @export
#' 

X.calibrate.scores <- function(
    data=NULL,
    standard.set=NULL,
    scores.col=NA,
    labels.col=NA)
{
  
  customPlot <- list(
    theme_bw(base_size = 12),
    # scale_fill_brewer(palette = "Set1"),
    # scale_colour_brewer(palette = "Set1"),
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())
  )
  
  if(is.null(data)) {
    stop("Please specify data.")
  } else {
    if(is.na(scores.col)) {
      scores.col <- names(data)[ncol(data)]
      cat("No scores column selected. Selecting the last column:",scores.col,"\n")
    }
    data <- as.data.frame(data) %>%
      rename(score=!!sym(scores.col))
  }
  if(!is.null(standard.set)) {
    standard.set <- as.data.frame(standard.set)
    if(is.na(labels.col)) {
      labels.col <- names(standard.set)[ncol(standard.set)]
      cat("No labels column selected. Selecting the last column:",labels.col,"\n")
    }
  }
  calibration_data <- standard.set %>% 
    cross_join(data, vars=c("protein1","protein2"), mode="left") %>%
    na.omit() %>%
    mutate(bin=cut(score, breaks=0:49/50, right=TRUE)) %>%
    arrange(bin)  %>%
    group_by(bin) %>%
    summarise(sum=sum(complex), count=n(), prob=sum(complex)/n()) %>%
    mutate(midpoint=(0:49/50+0.01)[1:nrow(.)])
  tofit.data <- calibration_data%>%dplyr::select(midpoint,prob)%>%setNames(c("x","y"))
  prob_model <- fit_sigmoid01(tofit.data)
  # prob_model <- nls(formula=y~(0+(1/(1+exp(b*(log(x)-log(e)))))), data=tofit.data,
  #                   start=list(b=-4,e=0.1))
  predict_line <- data.frame(x=1:10000/10000)
  predict_frame <- data.frame(x=predict_line, y=predict(prob_model,predict_line)) %>%
    setNames(c("midpoint","prob"))
  
  calibration_plot <- calibration_data %>%
    ggplot(aes(x=midpoint,y=prob)) +
    geom_point() + geom_line() +
    geom_line(inherit.aes=F, data=predict_frame, aes(x=midpoint,y=prob), color="red") +
    customPlot +
    geom_abline(intercept=0)
  
  data$probability <- predict(prob_model, data.frame(x=data$score))
  
  calibration_data2 <- standard.set %>% 
    cross_join(data, vars=c("protein1","protein2"), mode="left") %>%
    na.omit() %>%
    mutate(bin=cut(probability, breaks=0:49/50, right=TRUE)) %>%
    arrange(bin)  %>%
    group_by(bin) %>%
    summarise(sum=sum(complex), count=n(), prob=sum(complex)/n()) %>%
    mutate(midpoint=(0:49/50+0.01)[1:nrow(.)])
  
  calibrated_plot <-
    ggplot() +
    geom_abline(intercept=0) +
    geom_point(data=calibration_data, aes(x=midpoint,y=prob),color="red") + 
    geom_line(data=calibration_data, aes(x=midpoint,y=prob),color="red") + 
    geom_point(data=calibration_data2, aes(x=midpoint,y=prob),color="blue2") + 
    geom_line(data=calibration_data2, aes(x=midpoint,y=prob),color="blue2") +
    customPlot +
    xlab("Bin midpoint") + ylab("Observed event percentage")

  output <- list(data=data,plot=calibrated_plot,model=prob_model)
  return(output)
}
