#' X.fit.meltcurves
#'
#' Fit scaled protein melting curves to output fitting parameters.
#' 
#' @param data Data frame with one column for each temperature, starting with 'T' and followed by number and a column 'protein' with protein IDs.
#' 
#' @import tidyverse
#' @import e1071
#' @import minpack.lm
#' 
#' @return A list with four elements: $data with fitted data parameters, $all.data with connected table of input and output data,
#' $fits with nls fit objects for each fitted protein and $plot showing the distribution of these parameters across dataset.
#' @examples 
#' MC.data.fitted <- X.fit.meltcurves(MC.data.scaled$data)
#' @export

X.fit.meltcurves=function(
    data=NULL
) {
  
  if(is.null(data)) {
    stop("Please include data")
  } else {
    data <- as.data.frame(data)
    if(!"protein" %in% names(data)) {
      stop("Data must contain a column 'protein'")
    }
    data <- data %>% column_to_rownames("protein") %>%
      dplyr::select(starts_with("T"))
    if(nrow(data)<6) {
      stop("Less than five temperature columns found. The temperature column names must start with 'T' followed by a number.")
    }
    temperatures <- names(data) %>% str_remove("T") %>% as.numeric()
  }
      
  #initiate variables:
 
  outlier<-list()
  all_fits <- list()
  
  result.table <- data.frame()
  pb <- txtProgressBar(min=0, max=nrow(data), style=3, initial="") # progress bar
  cat(paste0("Fitting the melting curves...\n"))
  
  for(i in 1:nrow(data)) {
    
    # prepare data for the row i
    y <- data[i,] %>% as.list() %>% unname() %>% unlist()
    valueindex = which(!is.na(y)) # identify which temperatures contain non-NA values (some proteins are missing information for some temperatures)
    tempvector=temperatures[valueindex] # tempvector are values of temperatures that are not NA
    yvector <- y[valueindex]
    y_dif_saved=0
    
    # fit
    fitdata=data.frame(x=tempvector,y=yvector)
    sigmoid <- fit_sigmoid(fitdata)
    initial_sigmoid <- sigmoid
    
    # outlier removal
    R2Tm_table <- data.frame("point"=0, "R2"=R2nls(sigmoid,fitdata)$R2,"Tm"=calcTm(sigmoid))
    sigmoids<-list()
    
    if(R2Tm_table$R2[1]>=0.95) {
      outlier[[i]] <-c(NA,NA)
      penalty_score <- NA
      R2<-R2Tm_table$R2[1]
      R2_orig <- 1 - sum((residuals(sigmoid)^2))/sum((yvector-mean(yvector))^2)
    } else {
      for(outpar in seq_along(tempvector)) {
        yOut <- range.scale(yvector[-outpar])
        tOut <- tempvector[-outpar]
        outdata <- data.frame(x=tOut,y=yOut)
        sigmoid <- fit_sigmoid(outdata)
        sigmoids[[outpar]] <- sigmoid
        R2Tm_table <- rbind(R2Tm_table,c(outpar,R2nls(sigmoid,outdata)$R2,calcTm(sigmoid)))
      }
      R2Tm_table <- R2Tm_table %>% na.omit() 
      if(nrow(R2Tm_table)>0) {
        excluded_line <- grep(max(R2Tm_table$R2),R2Tm_table$R2)
        excluded_point <- R2Tm_table$point[excluded_line]
        R2<-R2Tm_table$R2[excluded_line]
        if(excluded_point==0) {
          sigmoid <- initial_sigmoid
        } else {
          sigmoid <- sigmoids[[excluded_point]]
        }
        maxR2_Tm <- R2Tm_table$Tm[excluded_line]
        R2_orig <- 1 - sum((residuals(sigmoid)^2))/sum((range.scale(y[-excluded_point])-mean(range.scale(y[-excluded_point])))^2)
        
        R2Tm_stattable <- R2Tm_table %>% 
          filter(point!=excluded_point) %>%
          mutate(R2=ifelse(R2<0,0,R2))
        
        penalty_score <- abs(weighted.mean(R2Tm_stattable$Tm,R2Tm_stattable$R2)/maxR2_Tm-1)
        if(is.nan(penalty_score) | is.na(penalty_score)) {
          penalty_score=1
        }
        if(penalty_score>1) {
          
          penalty_score <- 1-(penalty_score-1)
        }
        
        outlier[[i]] <- c(tempvector[excluded_point],yvector[excluded_point])
        
        # if the removed point is lowest or highest y value, the meltcurve needs to be re-scaled.
        if(excluded_point!=0) {
          if( yvector[excluded_point]==max(seq_along(tempvector)) | yvector[excluded_point]==min(seq_along(tempvector)) ) {
            yvector[which(tempvector==outlier[[i]][1])]<-NA
            yvector <- range.scale(yvector)
            for(val in seq_along(yvector)) {
              data[i,val] <- yvector[val]
            }
            outlier[[i]] <- c(NA,NA)
          }
        } else {
          outlier[[i]] <- c(NA,NA)
        }
        
      } else {
        outlier[[i]] <- c(NA,NA)
      }
      
    }
    
    # save the fitting coefficients and fits of each protein into result lists:
    if (class(sigmoid)=="nls") {
      
      all_fits[[rownames(data)[i]]] <- sigmoid
      
      coeffs <- data.frame(coefficients(sigmoid))
      fit.coefs <- c(
        "slope"= coeffs["b",1],
        "plateau"=coeffs["c",1],
        "top"=coeffs["d",1],
        "tilt"=coeffs["a",1],
        "Ti"=coeffs["e",1],
        "Tm"=calcTm(sigmoid),
        "origR2"=R2_orig,
        "R2"=R2,
        "Penalty"=penalty_score
      )
      
      
    } else {
      fit.coefs=rep(NA,9) %>%setNames(c("slope","plateau","top","tilt","Ti","Tm","origR2","R2","Penalty"))
    }
    
    result.table <- bind_rows(
      result.table,
      t(fit.coefs) %>% as.data.frame() %>% add_column("colname"=rownames(data)[i]) %>% column_to_rownames("colname")
      )
    
    setTxtProgressBar(pb, i)
    
  }
  close(pb) # to close the progress bar
  
  #amend outliers
  result.table$outlier_x <- unlist(lapply(outlier, function(x) x[1]))
  result.table$outlier_y <- unlist(lapply(outlier, function(x) x[2]))
  
  #set unrealistic Tm and Ti out
  maxT <- max(temperatures)
  minT <- min(temperatures)
  result.table$Tm[result.table$Tm<minT] <- NA
  result.table$Ti[result.table$Ti<minT] <- NA
  result.table$Tm[result.table$Tm>maxT] <- NA
  result.table$Ti[result.table$Ti>maxT] <- NA
  
  #calculate AUC
  result.table$AUC <- apply(data, 1, function(x) log10(MESS::auc(temperatures, x, type="linear")))
  data <- data %>% rownames_to_column("protein") %>% dplyr::select(protein,everything())
  
  #prepare plots
  predictors <- c("slope","plateau","top","tilt","Ti","Tm","origR2","R2","Penalty","AUC")
  skew_table <- result.table%>% select(where( ~!all(is.na(.x)))) # to remove if some predictor only has NA values
  predictors <- intersect(predictors,names(skew_table))
  distribution_plots <- list()
  for (i in predictors) {
    predvector <- skew_table%>%pull(!!sym(i))%>%na.omit()%>%as.vector()
    if(length(predvector)<3) {
      next
    }
    skewness <- e1071::skewness(predvector,type=2)
    if(is.na(skewness)) {
      next
    }
    if(skewness<(-1)) {
      trafo="probit"
    } else if(skewness>1) {
      trafo="sqrt"
    } else {
      trafo="identity"
    }
    
    distribution_plots[[i]] <-  result.table %>% 
      rownames_to_column("id") %>%
      ggplot(aes(x=!!sym(i))) + geom_density() + 
      theme_bw(base_size = 12) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank()
      ) +
      scale_x_continuous(trans=trafo)
    
  }
  
  # construct a table with scaled abundance values and fit parameters:
  all_data <- bind_cols(data,result.table)
  result.table <- result.table %>%
    rownames_to_column("protein") %>%
    dplyr::select(protein,everything())
  
  output <- list(data=result.table,all.data=all_data,fits=all_fits,plots=distribution_plots)
  return(output)
}

