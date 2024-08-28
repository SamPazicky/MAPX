#' X.scale.meltcurves
#'
#' Scale measured protein abundances across temperature gradient.
#' 
#' @param data Data frame with data containing a column 'id' and one 'Ab' column for each temperature.
#' @param temperatures Vector of integers: Temperatures used for the generation of the meltcurves.
#' 
#' @import tidyverse
#' @import minpack.lm
#' 
#' @return A list with two elements: $data with scaled data and $plot showing the data distribution before and after the scaling.
#' @examples 
#' MCdata.scaled <- X.scale.meltcurves(MCdata.clean$data)
#' @export


X.scale.meltcurves <- function(
    data=NULL,
    temperatures=c(37,41,45,49,53,57,61,65,69,73)
) {
  
  customPlot <- list(
    theme_bw(base_size = 12),
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.margin=ggplot2::margin(5,5,5,5, "pt")
    )
  )
  
  if(is.null(data)) {
    stop("Please load the data.")
  } else {
    if(!"protein" %in% names(data)) {
      stop("The data has to contain the column 'protein' with protein IDs.")
    }
    for(temp in seq_along(temperatures)) {
      if(!paste0("Ab",temp) %in% names(data)) {
        stop("The data must contain as many 'Ab' columns as number of temperatures")
      }
    }
    data <- data %>% dplyr::select(protein,starts_with("Ab"))
  }
  
  
  data <- as.data.frame(data) %>%
    dplyr::select(protein,everything()) %>%
    setNames(c("protein",paste0("T",temperatures)))
  
  scaled.data <- data %>%
    column_to_rownames("protein") %>%
    (t) %>% as.data.frame() %>%
    mutate(across(everything(), range.scale)) %>%
    t() %>% as.data.frame()
  
  medians <- data %>%
    select(!protein) %>%
    summarise(across(everything(), ~ median(.x,na.rm=TRUE)))
  
  fitdata <- medians %>% t() %>% as.data.frame() %>% 
    rownames_to_column("T") %>% mutate(T=str_remove(T,"T")) %>% mutate(across(everything(), as.numeric)) %>% 
    setNames(c("x","y"))
  
  fit <- fit_median_sigmoid(fitdata)
  if(class(fit)=="try-error") {
    stop("Median of raw data cannot be fitted with sigmoidal curve. Check the quality of your data!")
  }
  
  fitted.y <- predict(fit,data$x)
  fitting.factor <- fitted.y/medians
  scaling.factor <- 1/fitted.y[1]
  

  fitscaled.data <- data.frame()
  rownames(data) = data$protein
  for(i in 1:nrow(data)) {
    fitscaled.data <- bind_rows(
      fitscaled.data,
      (data%>%dplyr::select(!protein)%>%slice(i))*fitting.factor*scaling.factor
    )
  }
  
  fitscaled.data <- fitscaled.data %>%
    (t) %>% as.data.frame() %>%
    mutate(across(everything(), range.scale)) %>%
    t() %>% as.data.frame()
  
  boxplot <- bind_rows(
    scaled.data %>% rownames_to_column("protein") %>% mutate(Scaling="Before"),
    fitscaled.data %>% rownames_to_column("protein") %>% mutate(Scaling="After")
  ) %>%
    na.omit() %>%
    pivot_longer(cols=!c(protein,Scaling), names_to="Temperature", values_to="Fraction_soluble") %>%
    mutate(Temperature=str_remove(Temperature,"T")) %>%
    mutate(Scaling=factor(Scaling,levels=c("Before","After"))) %>%
    ggplot(aes(x=Temperature,y=Fraction_soluble)) +
    geom_boxplot(outlier.shape=NA) +
    facet_wrap(~Scaling) +
    customPlot +
    scale_y_continuous(name="Fraction soluble")
  boxplot
  
  fitscaled.data <- fitscaled.data %>% rownames_to_column("protein") %>% dplyr::select(protein,everything())
  
  output<-list(data=fitscaled.data,
               plot=boxplot
  )
  return(output)
  
}
