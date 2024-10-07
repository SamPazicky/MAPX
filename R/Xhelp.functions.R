# simple functions for easier life with R

expdif <- function(x, y) {            # Custom function for expression differences
  z <- log2(abs(x - y))
  return(z)
}

absdif <- function(x, y) {
  z <- abs(x-y)
  return(z)
}

myproduct <- function(x, y) {
  z <- x*y
  return(z)
}

mysum <- function(x, y) {
  z <- x+y
  return(z)
}

mymean <- function(x,y) {
  z <- (x+y)/2
  return(z)
}

#' mylower
mylower <- function(x,y) {
  z <- min(c(x,y))
  return(z)
}
#' interleave
interleave <- function(x, y) { # to re-arrange two vectors such that their values alter proportionally.
  m <- length(x)
  n <- length(y)
  xi <- yi <- 1
  len <- m + n
  err <- len %/% 2
  res <- vector()
  for (i in 1:len)
  {
    err <- err - m
    if (err < 0)
    {
      
      res[i] <- x[xi]
      xi <- xi + 1
      err <- err + len
    } else
    {
      res[i] <- y[yi]
      yi <- yi + 1
    }
  }
  res
}

#' range.scale
range.scale <- function(x) {
  z <- (x-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))
  return(z)
}


#' my_get_legend
my_get_legend<-function(myggplot) { # to extract legend from a single ggplot
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


#' robust_dist
# robust distance metric for heatmap that ignores outliers:
robust_dist = function(x, y) {
  qx = quantile(x, c(0.1, 0.9))
  qy = quantile(y, c(0.1, 0.9))
  l = x > qx[1] & x < qx[2] & y > qy[1] & y < qy[2]
  x = x[l]
  y = y[l]
  sqrt(sum((x - y)^2))
}

#' localMaxima
localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

#' addspace

# add space
# function for case-insentive search
`%qsin%` <- function(str, vec) {
  tolower(str) %chin% na.omit(tolower(vec))
}

#' ffit
# function for calculating root of sigmolinear function
ffit <- function (tempvectorOut,a,b,d,e) {
  a*tempvectorOut+((d)/(1+exp(b*(log(tempvectorOut)-log(e)))))-d/2
}

#' list2df_tibble
list2df_tibble <- function(x) {
  tmp <- purrr::map(x, tibble::as_tibble)
  dplyr::bind_rows(tmp, .id = "name")
}

#' cross_join
#' @export
#' 

#I think this is the best I can have

cross_join <- function(
    data1,data2,vars,mode="full"
) {
  
  names1 <- data1 %>% dplyr::select(-all_of(vars)) %>% names()
  names2 <- data2 %>% dplyr::select(-all_of(vars)) %>% names()
  if(length(intersect(names1,names2))>0) {
    stop("Two cross-joined datasets cannot share column names.")
  }
  
  if(mode=="inner") {
    
    output1 <- inner_join(data1,data2,by = structure(names = vars, .Data =vars))
    output2 <- inner_join(data1,data2,by = structure(names = rev(vars), .Data =vars))
    output <- rbind(output1,output2)
  }
  
  if(mode=="left") {
    
    output1 <- inner_join(data1,data2,by = structure(names = vars, .Data =vars))
    output2 <- inner_join(data1,data2,by = structure(names = rev(vars), .Data =vars))
    output_inner <- rbind(output1,output2)
    
    output3 <- anti_join(data1,output_inner,by = structure(names = vars, .Data = vars))
    output <- bind_rows(output_inner,output3)
    
  }
  
  if(mode=="right") {
    
    output1 <- inner_join(data2,data1,by = structure(names = vars, .Data =vars))
    output2 <- inner_join(data2,data1,by = structure(names = rev(vars), .Data =vars))
    output_inner <- rbind(output1,output2)
    
    output3 <- anti_join(data2,output_inner,by = structure(names = vars, .Data = vars))
    output <- bind_rows(output_inner,output3)
    
  }
  
  if(mode=="full") {
    
    output1 <- inner_join(data1,data2,by = structure(names = vars, .Data =vars))
    output2 <- inner_join(data1,data2,by = structure(names = rev(vars), .Data =vars))
    output_inner_left <- rbind(output1,output2)
    
    output_anti_left <- anti_join(data1,output_inner_left,by = structure(names = vars, .Data = vars))
    
    output1 <- inner_join(data2,data1,by = structure(names = vars, .Data =vars))
    output2 <- inner_join(data2,data1,by = structure(names = rev(vars), .Data =vars))
    output_inner_right <- rbind(output1,output2)
    
    output_anti_right <- anti_join(data2,output_inner_right,by = structure(names = vars, .Data = vars))
    
    output_anti <- bind_rows(output_anti_left,output_anti_right)
    output <- bind_rows(output_inner_left,output_anti)
    
  }
  return(output)
}

#' expand.grid.unique
expand.grid.unique <- function(x, y, include.equals=FALSE) {
  x <- unique(x)
  y <- unique(y)
  g <- function(i) {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}

#' moderate.inf
moderate.inf <- function(data) {
  for (i in names(data)) {
    data <- data %>% mutate(!!i:=ifelse(is.infinite(get(i)), 
                                      ifelse(get(i)<0,
                                             min(get(i)[which(is.finite(get(i)))]),
                                             max(get(i)[which(is.finite(get(i)))])),
                                      get(i)))
  }
  
} 


#' fit.sigmoid
fit.sigmoid <- function(data,lower=c(-Inf,-Inf,-Inf,-Inf,-Inf),upper=c(Inf,Inf,Inf,Inf,Inf)) {
    if(length(data$x)!=length(data$y)) {
      stop("Vectors x and y do not have the same length.")
    }
    # guess the initial value of inflection point for fitting. 
    # The guess is made when difference between subsequent scaled abundances is >0.25
    for (point in 1:(length(data$y)-1)) {
      y_dif <- abs(data$y[point+1] - data$y[point])
      if(y_dif>0.25) {
        e_guess <- data$x[point+1]
        e_vec <- point+1
        break
      } else {
        y_dif_saved <- y_dif
        e_guess <- median(data$x)
        e_vec <- 5
      }
    }
    #guess slope
    b_guess <- ((data$y[e_vec-1]-data$y[e_vec+1])*10)^3/7
    d_guess <- max(data$y)
    c_guess <- min(data$y)
    # form=as.formula(paste0(y,"~(c+a*",x,"+((d-c)/(1+exp(b*(log(",x,")-log(e))))))"))
    assign("my.fit.dat",
           try(minpack.lm::nlsLM(formula=y~(c+a*x+((d-c)/(1+exp(b*(log(x)-log(e)))))),
                                 data=data,
                                 start=list(a=0,b=b_guess,c=c_guess,d=d_guess, e=e_guess),
                                 lower=c(-0.004,0,-0.35,0.85,37),
                                 upper=c(0.004,150,0.26,1.1,73),
                                 control=list(maxiter=100)),
               silent=TRUE)
    )
    
  return(my.fit.dat)
}


#' fit_model
fit_model <- function(data) {
  if(length(data$x)!=length(data$y)) {
    stop("Vectors x and y do not have the same length.")
  }
  # guess the initial value of inflection point for fitting. 
  # The guess is made when difference between subsequent scaled abundances is >0.25
  data_piv <- data %>% group_by(x) %>% summarise_at(vars(y), mean) %>% mutate(y=range.scale(y))
  e_table <- data.frame(x=data_piv$x[-length(data_piv$x)]+diff(data_piv$x)/2,y=diff(data_piv$y))
  e_fit <-try(minpack.lm::nlsLM(formula = y ~ b*x+a*x^2+c, data = e_table, start=list(a=1,b=1,c=0)),silent=TRUE)
  if(class(e_fit)!="trial-error") {
    coefs=coefficients(e_fit)
    e_guess=(-1)*coefs["b"]/(2*coefs["a"])
    e_vec=which.min((abs(data_piv$x-e_guess)))
  } else {
    for (point in 1:(length(data_piv$y)-1)) {
      y_dif <- abs(data_piv$y[point+1] - data_piv$y[point])
      if(y_dif>0.25) {
        e_guess <- data_piv$x[point+1]
        e_vec <- point+1
        break
      } else {
        y_dif_saved <- y_dif
        e_guess <- median(data_piv$x)
        e_vec <- 5
      }
    }
  }
 
  
  #guess slope
  b_guess <- ((data_piv$y[e_vec-1]-data_piv$y[e_vec+1])*10)^3/7
  #guess min and max
  data_piv <- data %>% group_by(x) %>% summarise_at(vars(y), mean)
  d_guess <- max(data_piv$y)
  c_guess <- min(data_piv$y)
  c_min <- min(data$y)
  d_max <- max(data$y)
  # guess tilt by checking first points
  # first_inc <- min(which(diff(data %>% arrange(x) %>% pull(y)) > 0))
  # a_guess <- data %>% arrange(x) %>% slice_head(n=2) %>% lm(formula='y~x') %>% .[["coefficients"]] %>% .[["x"]]
  
  # form=as.formula(paste0(y,"~(c+a*",x,"+((d-c)/(1+exp(b*(log(",x,")-log(e))))))"))
  assign("my.fit.dat",
         try(minpack.lm::nlsLM(formula=y~(c+a*x+((d-c)/(1+exp(b*(log(x)-log(e)))))),
                               data=data,
                               start=list(a=0,b=b_guess,c=c_guess,d=d_guess, e=e_guess),
                               lower=c(-0.004,-150,c_min,NA,4),
                               upper=c(0.004,150,NA,d_max,40),
                               control=list(maxiter=100)),
             silent=TRUE)
  )
  
  if(class(my.fit.dat)=="try-error") {
    assign("my.fit.dat",
           try(minpack.lm::nlsLM(formula = y ~ b*x+a*x^2+c, data = data, start=list(a=1,b=1,c=0)),
               silent=TRUE
           ),
           
    )
  }
  
  
  return(my.fit.dat)
}



#' list.dirs.wp
list.dirs.wp <- function(path=".", 
                         pattern=NULL,
                         recursive=FALSE) {
  
  # Load necessary package

  # Get a list of all directories
  all_folders <- list.dirs(path, recursive = recursive)
  
  # Filter the list of folders using the specified pattern
  if(is.null(pattern)) {
    matching_folders <- all_folders
  } else {
    matching_folders <- all_folders[grep(pattern, all_folders)]
  }
  return(matching_folders)
}


#' fit sigmoid function into meltcurve
fit_sigmoid <- function(data) {
  
  require(minpack.lm)
  
  if(length(data$x)!=length(data$y)) {
    stop("Vectors x and y do not have the same length.")
  }
  # guess the initial value of inflection point for fitting. 
  # The guess is made when difference between subsequent scaled abundances is >0.25
  for (point in 1:(length(data$y)-1)) {
    y_dif <- abs(data$y[point+1] - data$y[point])
    if(y_dif>0.25) {
      e_guess <- data$x[point+1]
      e_vec <- point+1
      break
    } else {
      y_dif_saved <- y_dif
      e_guess <- median(data$x)
      e_vec <- 5
    }
  }
  #guess slope
  b_guess <- ((data$y[e_vec-1]-data$y[e_vec+1])*10)^3/7
  d_guess <- max(data$y)
  c_guess <- min(data$y)
  # form=as.formula(paste0(y,"~(c+a*",x,"+((d-c)/(1+exp(b*(log(",x,")-log(e))))))"))
  assign("my.fit.dat",
         try(minpack.lm::nlsLM(formula=y~(c+a*x+((d-c)/(1+exp(b*(log(x)-log(e)))))),
                               data=data,
                               start=list(a=0,b=b_guess,c=c_guess,d=d_guess, e=e_guess),
                               lower=c(-0.004,0,-0.35,0.85,min(data$x)),
                               upper=c(0.004,150,0.26,1.1,max(data$x)),
                               control=list(maxiter=100)),
             silent=TRUE)
  )
  
  return(my.fit.dat)
}

#' fit median sigmoid function into meltcurve
fit_median_sigmoid <- function(data) {
  
  require(minpack.lm)
  
  if(length(data$x)!=length(data$y)) {
    stop("Vectors x and y do not have the same length.")
  }
  # guess the initial value of inflection point for fitting. 
  # The guess is made when difference between subsequent scaled abundances is 1/4 of expected range
  
  dif_gauge <- (max(data$y)-min(data$y))/4
  
  for (point in 1:(length(data$y)-1)) {
    y_dif <- abs(data$y[point+1] - data$y[point])
    if(y_dif>dif_gauge) {
      e_guess <- data$x[point+1]
      e_vec <- point+1
      break
    } else {
      y_dif_saved <- y_dif
      e_guess <- median(data$x)
      e_vec <- 5
    }
  }
  #guess slope
  d_guess <- max(data$y)
  c_guess <- min(data$y)
  # form=as.formula(paste0(y,"~(c+a*",x,"+((d-c)/(1+exp(b*(log(",x,")-log(e))))))"))
  assign("my.fit.dat",
         try(minpack.lm::nlsLM(formula=y~(c+((d-c)/(1+exp(b*(log(x)-log(e)))))),
                               data=data,
                               start=list(b=10,c=c_guess,d=d_guess, e=e_guess),
                               lower=c(0,NA,NA,min(data$x)),
                               upper=c(150,NA,NA,max(data$x)),
                               control=list(maxiter=100)),
             silent=TRUE)
  )
  
  return(my.fit.dat)
}


#' calcTm
calcTm <- function(fit) {
  if(class(fit)=="try-error") {
    return(NA)
  } else {
    coeffs <- data.frame(coefficients(fit))
    slope <- coeffs["b",1]#b
    plateau <- coeffs["c",1]#c
    maxprot <- coeffs["d",1] #d
    tilt <- coeffs["a",1] #a
    Ti <- coeffs["e",1] #e
    yvar <- formula(fit) %>% as.character() %>% .[2]
    
    eval(parse(text=paste0('
      fffit <- function(',yvar,',a,b,d,e) {
        a*',yvar,'+((d)/(1+exp(b*(log(',yvar,')-log(e)))))-d/2
      }')
      )
    )
    # ffit <- function (y,a,b,d,e) {
    #   a*y+((d)/(1+exp(b*(log(y)-log(e)))))-d/2
    # }
    
    Tm <- try((uniroot(fffit,interval=c(37,73),a=tilt,b=slope,e=Ti,d=maxprot, extendInt="yes"))$root,
              silent=TRUE)
    if(class(Tm)=="try-error") {
      Tm <- NA
    }
    return(Tm)
  }
  
}


#' R2NLS function for calculating R2 for non-linear relationship
R2nls <- function(nls.obj,fitdata) {
  if (class(nls.obj) != "nls") {
    return(list(R2=-100))
  } else {
    # da <- eval(nls.obj$data)
    da <- fitdata
    resp.name <- all.vars(summary(nls.obj)$formula)[1]
    form <- paste(resp.name, "~1", sep = "")
    # y<-get("y",envir=parent.env(environment()))
    m0 <- stats::lm(form, da)
    an <- stats::anova(nls.obj, m0)
    sqn <- stats::deviance(nls.obj)
    sqe <- stats::deviance(m0)
    r2 <- 1 - (sqn/sqe)
    aov <- data.frame(fv = c("regression", "residuals"),
                      gl = c(-an$Df[2], an$Res.Df[1]),
                      sq = c(-an$Sum[2], an$Res.Sum[1]))
    aov$qm <- aov$sq/aov$gl
    aov$F <- c(aov$qm[1]/aov$qm[2], NA)
    aov$"Pr(>F)" <- c(1 - stats::pf(aov$F[1],
                                    df1 = aov$gl[1],
                                    df2 = aov$gl[2]),
                      NA)
    names(aov) <- c(" ", "Df", "Sum Sq", "Mean Sq",
                    "F value", "Pr(>F)")
    return(list(anova = aov, R2 = r2))
  }
}

#' ffit function for searching root
ffit <- function (y,a,b,d,e) {
  a*y+((d)/(1+exp(b*(log(y)-log(e)))))-d/2
}

#' replace_inf for replacing infinite values
replace.inf <- function(x) {
  ifelse(is.infinite(x),
         ifelse(x < 0, min(x[is.finite(x)]-1, na.rm = TRUE), max(x[is.finite(x)]+1, na.rm = TRUE)),
         x
         )
}

#' to calculate euclidean distance between two vectors
e.dist <- function(vect1, vect2, weights=1) {
  result <- sqrt(sum(weights*(vect1 - vect2)^2))
  return(result)
}

#' to convert vector of costs into cost matrix
Xconvert_costs <- function(costs) {
  cost_mats=list()
  if(length(costs)==1) {
    if(is.na(costs)) {
      cost_mats[[1]] <- matrix(c(0, 1, 1, 0), nrow = 2)
      rownames(cost_mats[[1]]) <- colnames(cost_mats[[1]]) <- c(1,0)
      costs=1
    } else {
      cost_mats[[1]] <- matrix(c(0, costs, 1, 0), nrow = 2)
      rownames(cost_mats[[1]]) <- colnames(cost_mats[[1]]) <- c(1,0)
    }
  } else {
    for (i in seq_along(costs)) {
      cost_mats[[i]] <- matrix(c(0, costs[i], 1, 0), nrow = 2)
      rownames(cost_mats[[i]]) <- colnames(cost_mats[[i]]) <- c(1,0)
    }
  }
  return(cost_mats[[1]])
}

Xdownsample <- function(data,factor) {
  downdata <- data %>%
    filter(complex==0) %>%
    slice_sample(n=as.integer(round(nrow(.)/factor)))
  newdata <- data %>%
    filter(complex==1) %>%
    bind_rows(downdata)
  return(newdata)
}

fit_sigmoid01 <- function(data, b_guess=-25) {
  if(length(data$x)!=length(data$y)) {
    stop("Vectors x and y do not have the same length.")
  }
  # guess the initial value of inflection point for fitting. 
  # The guess is made when difference between subsequent scaled abundances is >0.25
  
  # outpoints <- mean(data$y[(length(data$y)-2):length(data$y)])
  
  c_guess <- 0
  d_guess <- 1
  
  
  scaledata=data
  scaledata$x<-range.scale(scaledata$x)
  
  # to guess parameter e, fit lm model to limited data
  lmdata <- scaledata %>% filter(y>0.2 & y<0.8) %>% setNames(c("y","x"))
  lmmodel <- lm(y~x,lmdata)
  e_guess <- predict(lmmodel, newdata=data.frame(x=0.5)) %>% unname()
  
  # for (point in 1:(length(scaledata$y)-1)) {
  #   y_dif <- abs(scaledata$y[point+1] - scaledata$y[point])
  #   print(y_dif)
  #   if(abs(y_dif)>0.25) {
  #     e_guess <- (data$x[point+1]+data$x[point])/2
  #     e_vec <- point+1
  #     break
  #   } else {
  #     y_dif_saved <- y_dif
  #     e_guess <- median(data$x)
  #     e_vec <- 5
  #   }
  # }
  #guess slope
  # b_guess <- ((data$y[e_vec-1]-data$y[e_vec])*10)^3/7
  # if(c_guess-d_guess>0) {
  #   b_guess=b_guess*(-1)
  # }
  # form=as.formula(paste0(y,"~(c+a*",x,"+((d-c)/(1+exp(b*(log(",x,")-log(e))))))"))
  assign("my.fit.dat",
         try(minpack.lm::nlsLM(formula=y~(0+((1)/(1+exp(b*(log(x)-log(e)))))),
                               data=data,
                               start=list(b=b_guess, e=e_guess),
                               lower=c(-100,0.0005),
                               upper=c(0,1),
                               control=list(maxiter=100)),
             silent=TRUE)
  )
  return(my.fit.dat)
}
