#' X.tune.tree
#'
#' Cross-validation and building of classification tree models for protein interaction prediction.
#' 
#' @param data Data frame with pairwise differential values. Must contain any number of columns with predictors for modelling and
#' a labels column with values 1 (for protein pairs that form a complex) and 0 (for protein pairs that do not form a complex).
#' @param labels.col Character string: A name of the data column that contains labels values (1 or 0). If not specified, last column will be considered as labels column.
#' @param predictors Vector of character strings: Names of the columns that contain predictors. If not specified, all columns except of labels column 
#' and columns with names starting with "protein" will be considered as predictors.
#' @param tree.type Character string: A type of tree model to build. Options are: CART, PART, J48, C5.0, C5.0Rules, RF.
#' @param cross.folds Numeric: To how many data pieces should the training data be split for cross-validation. Default is 5. The parameter is 
#' ignore when mode is set to 'm'.
#' @param eval.metric Character string: How should the model be evaluated in cross-validation? Defalt is "prc" for area under the precision-recall curve.
#' Other options are "roc" for area under the receiver-operator curve and "kappa" for Cohen's kappa.
#' @param CP Numeric vector: Complexity parameters to be cross-validated.
#' @param CF Numeric vector: Confidence factors to be cross-validated.
#' @param mTry A numeric vector: mtry parameter for random forest models. If not specified, five evenly spread mtry values between the number of predictors
#' and a fifth of that will be validated.
#' @param winnowing Vector of logicals: If winnowing should be validated, set to c(FALSE,TRUE). Default is FALSE.
#' @param boost Vector of integers: What boost iteratations should be used for cross-validation. Default is 1.
#' @param noGlobalPruning Vector of logicals: If global pruning should be validated, set to c(FALSE,TRUE). Default is FALSE.
#' @param costs Integer: A number specifying how much higher is the cost of falsely predicting non-interacting than interacting protein pairs and vice versa.
#' @param nodesize Vector of integers: Specifies minimums of terminal node sizes to be validated. Default is 10.
#' @param downsample Integer: How many times less of the non-interacting proteins should be used for the training? Default is 1. Applicable for C5.0 models.
#' @param rf.split Integer: To how many steps is the random forest training split? This is to avoid running out of RAM. Default is 5.
#' 
#' @import tidyverse
#' @import caret
#' @return A list with three elements. $data is a data frame with the tuning results, $eval.plot is the plotted data and $eval.plot.SDs are plotted standard deviations of the results.
#' @examples 
#' CV <- X.tune.tree(GS_specific,labels.col="complex",tree.type="RF",mTry=c(3,5),nodesize=c(10,20), downsample=3)
#' @export

X.tune.tree <- function(
  data=NULL,
  labels.col=NA,
  predictors=NA,
  tree.type="RF",
  cross.folds=5,
  eval.metric="prc",
  CP=0.01,
  CF=0.3,
  mTry=1,
  winnowing=FALSE,
  boost=1,
  noGlobalPruning=FALSE,
  costs=NA,
  nodesize=10,
  rf.split=5,
  downsample=1,
  ntree=151
) {
  
  # plot settings
  
  customPlot <- list(
    theme_bw(base_size = 12),
    # scale_fill_brewer(palette = "Set1"),
    # scale_colour_brewer(palette = "Set1"),
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())
  )
  
  Pal25 <- c(
    "dodgerblue2", "#E31A1C", # red
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "black", "gold1",
    "skyblue2", "#FB9A99", # lt pink
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
  )
  
  # Formatting data
  if(is.null(data)) {
    stop("Please define as data frame containing columns with predictors and label.")
  } else {
    data <- as.data.frame(data)
    if(is.na(labels.col)) {
      cat("Label column not selected. Selecting the last column.\n")
      labels.col <- names(data)[length(data)]
    }
    if(!labels.col %in% names(data)) {
      stop("Label must be a column name of the data.")
    } else {
      data <- data %>%
        rename(complex=!!sym(labels.col)) %>%
        mutate(complex=as.factor(complex))
    }
    if(length(setdiff(unique(data$complex),c(0,1)))>0) {
      stop("The label must contain values 0 and 1 for non-interacting and interacting protein pairs, respectively.")
    } else {
      labeltable <- unname(table(data$complex))
      cat(paste0("The model will be trained with ",labeltable[1], " non-interacting and ",labeltable[2]," interacting proteins.\n"))
    }
    if(any(grepl("protein",names(data)))) {
      data <- data %>%
        dplyr::select(!starts_with("protein"))
    }
    data <- data %>%
      dplyr::select(everything(),complex)
    if(is.na(predictors[1])) {
      predictors <- names(data)[-(length(data))]
      cat("Predictor columns not selected. Using following columns as predictors: ",paste(predictors,collapse=", "),"\n")
    } else {
      data <- data %>%
        dplyr::select(all_of(predictors),complex)
    }
  }
  
  

  if(tree.type %in% c("J48","PART")) {
    options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx32g"))
  } else if(tree.type=="RF") {
    cat("Random forest does not accept NA values. All NA values will be replaced with an arbitrary value of -1E6. \n")
    
  }
  
  ### This should be accessible all time anyways.
  tree.options <- list(
    CART=c("CP","costs"),
    J48=c("CF"),
    PART=c("CF"),
    C5.0=c("winnowing","CF","boost","noGlobalPruning","downsample"),
    C5.0Rules=c("winnowing","CF","boost","downsample"),
    RF=c("mTry","nodesize","rf.split","ntree","downsample")
  )
  ###
  
  cat("Preparing data...\n")
  # split the data for cross-validation
  traindata <- data %>%
    mutate(complex=factor(complex,levels=c(0,1))) %>%
    group_by(complex) %>%
    mutate(Folds=createFolds(1:n(),k=cross.folds,list=FALSE)) %>%
    ungroup()
  
  # create a loop table
  cur.options <- mget(tree.options[[tree.type]]) %>%
    append(list(CV=1:cross.folds)) %>%
    expand.grid() %>%
    as.data.frame() 
  
  functionname <- paste0("Xmodel.tree.",tree.type)
  
  cat(paste0("Tuning ", tree.type, " model parameters...\n"))
  pb <- txtProgressBar(min=0, max=nrow(cur.options), style=3, initial="") # progress bar
  
  metrics <- vector()
  for(i in 1:nrow(cur.options)) {
    cur.data <- cur.options[i,]
    cur.list <- cur.data %>% dplyr::select(!CV) %>% unlist() %>% as.list()
    cur.traindata <- traindata %>% filter(Folds!=cur.data$CV) %>% dplyr::select(!Folds)
    cur.evaldata <- traindata %>% filter(Folds==cur.data$CV) %>% dplyr::select(!Folds)
    
    modelresult <- do.call(functionname,list(data=cur.traindata) %>% append(cur.list))
    eval <- X.validate(cur.evaldata, modelresult$model, modelresult$predict.type, labels.col="complex", eval.metric=eval.metric)
    metrics[i] <- eval$eval.metric
    
    setTxtProgressBar(pb,i)
  }
  close(pb)
  
  result.table <- cur.options %>%
    mutate(metric=metrics) %>%
    mutate(success=ifelse(is.na(metric),0,1)) %>%
    mutate(metric=ifelse(is.na(metric),0,metric)) %>%
    group_by(across(!c(CV,metric,success))) %>%
    dplyr::summarise(metric.mean=mean(metric),metric.sd=sd(metric),success.rate=mean(success), .groups="keep") %>%
    ungroup() %>%
    dplyr::select(all_of(names(which(sapply(.,typeof)=="double"))),all_of(names(which(sapply(.,typeof)=="logical")))) %>%
    dplyr::relocate(starts_with("metric"), .after=last_col()) %>%
    dplyr::select(!success.rate)
  
  # plotting
  cat("Plotting\n")
  vars.toplot <- Filter(var,result.table) %>%
    dplyr::select(!starts_with("metric")) %>%
    names()
  
  if(length(vars.toplot)==1) {
    AUC_plot <- result.table %>%
      ggplot() + geom_point(mapping=aes(x=!!sym(vars.toplot[1]),y=metric.mean), size=3) +
      geom_line(mapping=aes(x=!!sym(vars.toplot[1]),y=metric.mean)) +
      geom_errorbar(mapping=aes(x=!!sym(vars.toplot[1]),ymin=metric.mean-metric.sd,ymax=metric.mean+metric.sd), 
                    width=(max(result.table[,1])-min(result.table[,1]))/nrow(unique(result.table[,1]))*0.15 ) +
      customPlot + ylim(c(0,1)) +
      scale_y_continuous(name=eval.metric)
    
    AUC_plot_SDs <- AUC_table %>%
      ggplot() + geom_point(mapping=aes(x=!!sym(current_tune_pars),y=sdAUC), size=3) +
      geom_line(mapping=aes(x=!!sym(current_tune_pars),y=sdAUC)) +
      customPlot + ylim(c(0,0.15)) +
      scale_y_continuous(name=paste0("SD(",eval.metric,")"))
  } 
  
  if(length(vars.toplot) >=2) {
    AUC_plot <- result.table %>%
      ggplot(mapping=aes(x=!!sym(vars.toplot[1]),y=metric.mean, color=as.character(!!sym(vars.toplot[2])))) + 
      geom_point(size=3,alpha=0.75) +
      geom_line(alpha=0.75) +
      geom_errorbar(mapping=aes(x=!!sym(vars.toplot[1]),ymin=metric.mean-metric.sd,ymax=metric.mean+metric.sd),
                    width=(max(result.table[,1])-min(result.table[,1]))/nrow(unique(result.table[,1]))*0.15, alpha=0.75) +
      customPlot +
      scale_color_manual(name=vars.toplot[2],values=Pal25) +
      scale_y_continuous(name=eval.metric, limits=c(0,1))
    
    AUC_plot_SDs <- result.table %>%
      ggplot(mapping=aes(x=!!sym(vars.toplot[1]),y=metric.sd, color=as.character(!!sym(vars.toplot[2])))) + 
      geom_point(size=3,alpha=0.75) +
      geom_line() +
      customPlot +
      scale_color_manual(name=vars.toplot[2],values=Pal25) +
      scale_y_continuous(name=paste0("SD(",eval.metric,")"), limits=c(0,0.2))
  } 

  if(length(vars.toplot)==3) {
    AUC_plot <- AUC_plot +
      facet_wrap(as.formula(paste0("~ ", vars.toplot[3]))) +
      ggtitle(vars.toplot[3]) +
      theme(plot.title=element_text(hjust=0.5))
    
    AUC_plot_SDs <- AUC_plot_SDs +
      facet_wrap(as.formula(paste0("~ ", vars.toplot[3]))) +
      ggtitle(vars.toplot[3]) +
      theme(plot.title=element_text(hjust=0.5))
  } 
  
  if(length(vars.toplot)==4) {
    AUC_plot <- AUC_plot +
      facet_wrap(as.formula(paste0(vars.toplot[3], " ~ ", vars.toplot[4]))) +
      ggtitle(paste0(vars.toplot[3],"\n",vars.toplot[4])) +
      theme(plot.title=element_text(hjust=0.5))
    
    AUC_plot_SDs <- AUC_plot_SDs +
      facet_wrap(as.formula(paste0(vars.toplot[3], " ~ ", vars.toplot[4]))) +
      ggtitle(paste0(vars.toplot[3],"\n",vars.toplot[4])) +
      theme(plot.title=element_text(hjust=0.5))
  }
  bestvalue <- max(result.table$metric.mean)
  cat(paste0("Done. The best model has a ",eval.metric," value of ",bestvalue,".\n"))
  output <- list(data=result.table,eval.plot=AUC_plot,eval.plot.SDs=AUC_plot_SDs)

  return(output)
}

