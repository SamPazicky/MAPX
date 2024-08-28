#' X.crosstrain.tree
#'
#' Train classification tree model with a cross-training algorithm.
#' 
#' @param data Data frame with all data coming from an experiment. Each row corresponds to a unique pair of proteins and
#' other columns must include all predictors that are in the standard set.
#' @param standard.set Data frame with predictors for pairs of proteins and with labels.
#' @param train.cycles Integer: How many times should the training be repeated? Default is 25.
#' @param train.split Integer: To how many individual training sets should the data be split in each training cycle? Default is 4.
#' @param labels.col Character string: 'data' column name with labels (1 for complex-forming, 0 for others.) Default is NA, in which case
#' the last column will be considered as the column with labels.
#' @param tree.type Character string: What type of classification tree should be used for the model training? Options are "J48","CART", "PART", "C5.0","RF".
#' @param evaluate Logical: Should each model be evaluated? Default is FALSE. If TRUE, evaluation plots will be saved in the save.dir and
#' @param eval.metric Character string: What type of evaluation metric should be used for evaluation? Options as in X.evaluate
#' a table with evaluation metrics will be outputted.
#' @param plot Logical: If evaluate is TRUE, should the plots be saved?
#' @param train.pars Named list: Tuning parameters for classification tree training. Default is list(), in which case the default tuning parameters
#' as specified in the arguments of the function X.model.tree will be used.
#' @param save.dir Character string: Name of the folder to be used to save the models. Default is "models".
#' @return A list with three elements. $data for all prediction scores, $metric.data for all underlying data for evaluation
#' metric calculation, $metric.plots for metric plots.
#' 
#' @import tidyverse
#' @import caret
#' @examples 
#' cross.model <- X.crosstrain.tree(data=all.predictors[[tp]][[rep]],standard.set=GS_specific,train.pars=best.train.pars,
#'         evaluate=TRUE,plot=TRUE,
#'         train.split=3, train.cycles=10)
#' @export

X.crosstrain.tree = function(
  data=NULL, # all data coming from an experiment
  standard.set=NULL, # training dataset with labels and protein pairs
  train.cycles=25,
  train.split=4,
  labels.col=NA, # column in standard set with labels
  tree.type="RF",
  evaluate=FALSE,
  eval.metric="prc",
  plot=TRUE,
  train.pars=list(), # named list
  save.dir="models",
  save.output=c("models","plots")
) {
  
  if(is.null(data)) {
    stop("Please include data")
  } else {
    data <- as.data.frame(data)
    
  }
  
  if(is.null(standard.set)) {
    stop("Plase include standard set")
  } else {
    standard.set <- standard.set %>% as.data.frame()
    if(is.na(labels.col)) {
      labels.col <- names(standard.set)[ncol(standard.set)]
      cat("No labels column selected. Selecting the last column:",labels.col,"\n")
    }
    standard.set <- standard.set %>%
      rename(complex=!!sym(labels.col))
    if(length(setdiff(unique(standard.set$complex),c(0,1)))>0) {
      stop("The label must contain values 0 and 1 for non-interacting and interacting protein pairs, respectively.")
    } else {
      labeltable <- unname(table(standard.set$complex))
      cat(paste0("The model will be trained with ",labeltable[1], " non-interacting and ",labeltable[2]," interacting proteins.\n"))
    }
  }
  
  if(tree.type %in% c("J48","PART")) {
    options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx32g"))
  } else if(tree.type=="RF") {
    cat("Random forest does not accept NA values. All NA values will be replaced with an arbitrary value of -1E6. \n")
    data <- data %>%
      mutate(across(where(is.numeric), ~ ifelse(is.na(.x),-1000000,.x))) %>%
      mutate(across(where(is.character), ~ ifelse(is.na(.x), "NAVALUE",.x)))
    
  }
  
  functionname <- paste0("Xmodel.tree.",tree.type)
  
  dir.create(save.dir, showWarnings=FALSE)
  
  
  cat("Training a", tree.type, "model for prediction of protein complexes...\n")
  metrics <- list()
  prediction.data <- list()
  metric.plots <- list()
  
  counter=1
  pb <- txtProgressBar(min=0, max=train.cycles*train.split, style=3, initial="") # progress bar
  
  for(repmod in 1:train.cycles) {
    
    metrics[[repmod]] <- list()
    prediction.data[[repmod]] <- list()
    metric.plots[[repmod]] <- list()
    # split the data for cross-training
    traindata <- standard.set %>%
      mutate(complex=factor(complex,levels=c(0,1))) %>%
      group_by(complex) %>%
      mutate(Folds=createFolds(1:n(),k=train.split,list=FALSE)) %>%
      ungroup()
    
    for(k.rep in 1:train.split) {
      
      traindata.repmod <- traindata %>% 
        filter(Folds!=k.rep) %>%
        dplyr::select(!Folds)
      
      predictdata.repmod <- traindata.repmod %>% 
        dplyr::select(starts_with("protein"),complex) %>% 
        cross_join(data,vars=c("protein1","protein2"),mode="right") %>%
        filter(is.na(complex)) %>%
        dplyr::select(!complex)
      
      
      invisible(capture.output(
        modelresult <- do.call(functionname,list(data=traindata.repmod)%>%append(train.pars))
      ))
      
      if(length(modelresult[[1]])==1) {
        print(paste0("......Failed iteration ", repmod,", ", k.rep,"/",k))
        setTxtProgressBar(pb, counter)
        counter=counter+1
        next
      }
      saveRDS(modelresult$model,file=paste0(save.dir,"/model_",tree.type,"_repmod",repmod,"_k",k.rep,".RDS"))
      allpreds <- X.validate(predictdata.repmod, modelresult$model, modelresult$predict.type, labels.col="complex", eval.metric=eval.metric,
                         evaluate=FALSE)
      
      if(evaluate) {
        evaldata.repmod <- traindata %>%
          filter(Folds==k.rep) %>%
          dplyr::select(!Folds)
        
        eval <- X.validate(evaldata.repmod, modelresult$model, modelresult$predict.type, labels.col="complex", eval.metric=eval.metric,
                           evaluate=evaluate,plot=plot)
        
        if(plot) {
          suppressMessages(ggsave(paste0(save.dir,"/model_",tree.type,"_repmod",repmod,"_k",k.rep,"_",eval.metric,".png"),eval$plot))
          metric.plots[[repmod]][[k.rep]] <- eval$plot

        }
        metrics[[repmod]][[k.rep]] <- eval$eval.metric
      }
      
      prediction.data[[repmod]][[k.rep]] <- allpreds$data %>%
        dplyr::select(starts_with("protein"), score) %>%
        mutate(cycle=repmod) %>%
        mutate(split=k.rep)
      
      setTxtProgressBar(pb, counter)
      counter=counter+1
    }
  }
  close(pb)
  
  if(evaluate) {
    if(plot) {
      all.metric.plots <- lapply(metric.plots, function(x) x%>%setNames(paste0("split",1:train.split))) %>% setNames(paste0("cycle",1:train.cycles)) %>%
        unlist(recursive=FALSE) %>% setNames(str_replace(names(.),"\\.","_"))
    }
    
    all.metrics.data <- lapply(metrics, function(x) x%>%setNames(1:train.split)) %>% setNames(1:train.cycles) %>%
      unlist(recursive=FALSE) %>% stack() %>% setNames(c(eval.metric,"cyclesplit")) %>%
      separate_wider_delim(cols="cyclesplit",delim=".",names=c("cycle","split")) %>%
      dplyr::select(cycle,split,all_of(eval.metric))
  }
  
  all.prediction.data <- Reduce(bind_rows,prediction.data) %>%
    dplyr::select(starts_with("protein"),cycle,split,score)
  
  output <- list(
    data=all.prediction.data,
    metric.data=all.metrics.data,
    metric.plots=all.metric.plots
  )
  return(output)
}
