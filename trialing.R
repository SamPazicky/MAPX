library(data.table)
library(tidyverse)
library(MAPX)

reps <- c(1,2,3)
all.predictors <- list()
all.features <- list()
all.fitdata <- list()
all.fits <- list()
all.sdata <- list()
timepoints = paste0("tp",c(28,40))

for(tp in timepoints) {
  all.features[[tp]] <- list()
  all.predictors[[tp]] <- list()
  all.fitdata[[tp]] <- list()
  all.fits[[tp]] <- list()
  all.sdata[[tp]] <- list()
  
  for(rep in reps) {
    
    # 1. Loading the raw data.
    # First, load your data in R and configure it such that it is readable by MAPX commands.
    # The MAPX pipeline working with melting curve data expects the following column names:
    # - protein: column with protein IDs
    # - description: protein descriptions, such as in fasta headers
    # - length: length of the protein (number of amino acids)
    # - Ab1 to AbN where N is the number of temperatures measured: columns with raw abundances
    
    MCdata.raw <- meltdata[[tp]][[rep]] %>%
      # next two lines select the columns with protein IDs, their amino acid length and abundance values across the temperature gradient
      dplyr::select(Accession,Description,"# AAs",starts_with("Abundances (Grouped)")) %>% 
      dplyr::select(!contains("CV")) %>%
      # next row renames the columns to id, length and Ab1 to Ab10
      setNames(c("protein","description","length",paste0("Ab",1:10))) %>%
      # next row shortenst the IDs to remove alternative transcripts. This only simplifies the data and does not need to be done.
      mutate(protein=str_remove(protein,".[[:digit:]]-p[[:digit:]]$"))
    
    # 2. Data clean-up.
    # X.clean.meltcurves will remove contaminants, rows with NA values and average duplicate IDs.
    MCdata.clean <- X.clean.meltcurves(MCdata.raw)
    
    # The removed proteins can be viewed as a list:
    MCdata.clean$removed
    
    # The cleaned data is stored in a data frame:
    MCdata.clean$data %>% head()
    
    # The cleaned raw data is used for two things: Extraction of raw features for machine learning
    # and further melting curve processing (scaling and fitting).
    
    # 3a. Extraction of raw data features
    features.raw <- X.extract.raw.features(MCdata.clean$data, abundances=c(1,2), remains=10)
    head(features.raw)
    
    # 3b. Scaling of melting curves
    # The cleaned data are scaled from 0-1 and normalized into a sigmoid curve.
    MCdata.scaled <- X.scale.meltcurves(MCdata.clean$data)
    
    all.sdata[[tp]][[rep]] <- MCdata.scaled$data
    # The scaled data distribution before and after the normalization is plotted as a ggplot that
    # can be saved.
    # MCdata.scaled$plot
    # ggsave("MCdata_normalization.png",MCdata.scaled$plot)
    
    # normalized data are present in a data frame as:
    MCdata.scaled$data %>% head()
    
    # 4. Fitting of melting curves
    # For this function, the input data frame must contain the column 'protein' and one column for
    # each temperature starting with 'T' followed by a number.
    MCdata.fitted <- X.fit.meltcurves(MCdata.scaled$data)
    
    # The output gives a list of four elements. 
    # The first element contains fitted values. Some of these can be used as features for complex prediction.
    MCdata.fitted$data %>% head()
    
    # The fitdata will be important in later steps, let's save it in a list:
    all.fitdata[[tp]][[rep]] <- MCdata.fitted$data
    all.fits[[tp]][[rep]] <- MCdata.fitted$fits
    
    # The second element contains fitted values together with the original scaled and normalized abundance values
    MCdata.fitted$all.data %>% head()
    
    # The third element contains sigmoid fits for all melting curves. These are of class nls.
    MCdata.fitted$fits %>% head()
    
    # The fourth element contains plots that show distribution of the fitted and calculated parameters
    MCdata.fitted$plots %>% patchwork::wrap_plots()
    
    # The fitted melting curves can be visualized easily.
    # MCdata.curves <- X.plot.meltcurves(MCdata.scaled$data, MCdata.fitted$fits)
    
    # 5. Join raw features and fit features
    replace.inf <- function(x) {
      ifelse(is.infinite(x),
             ifelse(x < 0, min(x[is.finite(x)]-1, na.rm = TRUE), max(x[is.finite(x)]+1, na.rm = TRUE)),
             x
      )
    }
    
    features <- MCdata.fitted$data %>%
      full_join(features.raw) %>%
      mutate(logABL=log2(ABL)) %>%
      mutate(logAUC=log10(AUC)) %>%
      mutate(logLoss=log10(Loss)) %>%
      mutate(Penalty=ifelse(is.na(Penalty),2,Penalty)) %>%
      mutate(Penalty_trans=log(Penalty)) %>%
      mutate(Penalty_trans=ifelse(is.infinite(Penalty_trans)&Penalty_trans<0,min(Penalty_trans[which(is.finite(Penalty_trans))]),Penalty_trans)) %>%
      mutate(across(where(is.numeric), ~ replace.inf(.))) %>%
      select(protein,Ti,logAUC,logABL,logLoss,Penalty_trans)
    

    all.features[[tp]][[rep]] <- features
    # The features define each single proteins. Relating the proteins to each other for each protein
    # will the calculation of the predictors.
    
    # 6. Calculation of the predictors
    # In this function, 'data' is a data frame that has to contain a column protein and one column
    # for each element in the vector given in 'features'. 'funs' are functions that will be applied
    # to calculate for each feature and 'prefixes' will be used to names the reuslting columns.
    # For example, in the code below, the first feature is 'Ti' and so the function applied to 
    # 'Ti' between each pair of proteins will be absdif. Use ?X.calculate.predictors to see what
    # functions can be applied. The resulting column will be 'dTi' as the first prefix given
    # is 'd'.
    
    feature.predictors <- X.calculate.predictors(data=features, 
                                             features=c("Ti","logAUC","logABL","logLoss","Penalty_trans"), 
                                             funs=c(rep("absdif",4),"Xsum"),
                                             prefixes=c(rep("d",4),"sum")
    )
    
    # The previous function calculates 5 predictors from the features. An additional predictor is just a correlation
    # coefficient between data points of two proteins.
    
    cordata <- X.calculate.cors(MCdata.scaled$data)
    
    # The data frames feature.predictors and cordata should contain the same number of rows and rbind should be sufficient
    # to connect them. MAPX offers a generalized function to connect two data frames that contain two columns that are 
    # interchangeable: MAPX::cross_join.
    
    all.predictors[[tp]][[rep]] <- MAPX::cross_join(cordata,feature.predictors,vars=paste0("protein",c(1,2)),mode="full")
    
  }
}

X.plot.meltcurves(all.sdata %>% unlist(recursive=FALSE),all.fits %>% unlist(recursive=FALSE),
                  samples=rep(c("tp28","tp40"),each=3),replicates=c(1,2,3,1,2,3))

# At this point, the data has been reduced into predictors that can be fed into the machine learning model. 

# 7. Assemble gold standard dataset
# To train the model, gold standard dataset of interacting and non-interacting proteins has to be available.

# 7a. Assemble positive gold standard
# The positive gold standard is a pairwise list of protein-protein interactions that are already known. This,
# besides the data quality, is the most important model input, because the algorithm will search for protein 
# interactions with similar properties as those of protein-protein interactions included in the positive gold
# standard. A care should be taken when putting this together: presumably, MAP-X only detect stable and 
# stechiometrically defined interactions.
# Most databases list gold standard interactions with one column for complex name and one column for protein ID.
# The following command can take in this kind of data and convert it into a pair-wise interaction table.

GSpos <- X.assemble.posstandard(gold.standard, sep=";")

# 7b. Assemble negative gold standard
# The most comprehensive approach to assemble the gold standard of non-interacting proteins is using GO terms: the proteins
# that do not interact cannot share GO terms. A care should be taken when predicting interactome for organisms without
# well annotated proteome - in such cases, using GO terms from a similar organism basde on orthology would be a better
# strategy. The function expects a data frame with the first column containing protein IDs and the second column
# containing comma-separated GO terms.

# the following lines take a database from plasmoDB and prepares it for MAPX commands
GOterms <-  plasmoDB_data %>%
  dplyr::select(source_id,contains("GO", ignore.case=FALSE)) %>% # selects source_id column and all GO columns
  dplyr::select(contains("id")) %>% # selects all columns that contain id (gets rid of GO term descriptors)
  tidyr::unite("GO",contains("GO"),sep=";") %>% # connects all GO terms to semi-colon separate column "GO"
  #option 1: do not use any proteins that have N/A in at least on of the GO terms
  mutate(GO=ifelse(str_detect(GO,"N/A"),"",GO)) %>%
  #option 2: use everything. The algorithm will exclude those that have no GO term at all.
  # mutate(GO=ifelse(GO=="N/A;N/A;N/A","N/A",GO)) %>% # converts rows without GO to N/A
  # mutate(GO=str_remove_all(GO,"N/A;")) %>% # removes all N/As followed by semi-colon
  # mutate(GO=str_remove_all(GO,";N/A")) %>% # removes all N/As at the and of the row
  mutate(GO=ifelse(GO=="N/A","",GO)) %>%
  rename(protein=source_id) %>%
  mutate(protein=str_remove(protein,".[[:digit:]]$")) # removes '.1' at the end of the IDs

# mind that the function to assemble the negative pair-wise table is slow.
GSneg <- X.assemble.negstandard(data=GOterms, sep=";") %>%
  mutate(complex=0)

# finally, the two can be connected together
GS <- rbind(GSpos,GSneg)

# Let's build a global protein-protein interaction network for tp28.
comp.models <- list()
tp="tp28"
for(rep in reps) {
  # 8. Model tuning
  # Machine learning models can be calculated with different model parameters. The best choice of the parameters
  # depends on the nature of the dataset. In this step, the best parameters for particular model can be chosen.
  # MAP-X comes with commands for tree-type models. The best-working model in my hands was random forest,
  # but any type of model could theoretically be integrated into the pipeline.
  
  # First, extract from gold standard (GS) those pairs of proteins that you have data from:
  GS_specific <- X.assemble.traindata(all.predictors[[tp]][[rep]],GS)
  
  # This specific standard will be used to tune the parameters. In the example below, I chose random forest tree type (
  # tree.type="RF") and I optimize the parameters mTry with values 3, 5, 10 and 20 and nodesizes with 10 and 20. This
  # gives a total of 8 combinations and thus there will be 8 tuning steps.
  
  CV <- X.tune.tree(GS_specific,labels.col="complex",tree.type="RF",mTry=c(3,5),nodesize=c(10,20), downsample=3)
  
  # Let's look at the resulting data frame and choose the best model based on the chosen metric (in our case area under the precision-recall curve).
  # We can take the row with the best metric and extract the parameters as a named list - this named list is an argument for final model training.
  best.train.pars <- CV$data %>% arrange(desc(metric.mean)) %>% dplyr::select(!starts_with("metric")) %>% slice(1) %>% as.vector()
  
  # 9. X.crosstrain.tree is used for final model training. The function goes through a specified amount of train cycles. In each train cycle,
  # the data training data is randomly split to train.split number of subdata and one model is trained for each split. This results
  # in train.cycle x train.split number of models. The example below goes through 2 train cycles split to 3, resulting in 6 models.
  
  cross.model <- X.crosstrain.tree(data=all.predictors[[tp]][[rep]],standard.set=GS_specific,train.pars=best.train.pars,
                                   evaluate=TRUE,plot=TRUE,
                                   train.split=3, train.cycles=10)
  
  # If evaluate=TRUE, each submodel will be evaluated against the unused portion of data. If plot=TRUE, the evaluation metric will also be plotted.
  # The evaluation metric data as well as plot are outputted and can be saved:
  fwrite(cross.model$metric.data,"model_evals.csv")
  for(i in names(cross.model$metric.plots)) {
    plot <- cross.model$metric.plots[[i]]
    ggsave(paste0("metric.plot_",i,".png"),plot)
  }
  
  # 10. The submodels can be averaged into a composite model using X.average.crossmodels
  comp.models[[rep]] <- X.average.crossmodels(data=cross.model$data, eval.metric="prc", evaluate=TRUE, standard.set=GS_specific)
  
}

# 11. Finally, the final model is the result of averaging of composite models of multiple replicates.

comp.models.data <- lapply(comp.models, function(x) x$data)
final.model <- X.average.reps(data=comp.models.data, evaluate=TRUE,standard.set=GS)

# Next, let's build  a global protein interaction network. For this, a cut-off value for the data needs to be decided.
# To do this, let's loop through multiple precision values

statsplots <- list()
statsdata <- list()
final_networks <- list()

precision_scan <- seq(0.2,0.9,0.02)
for(ps in precision_scan) {
  
  # lets find the closest score value to the current precision value
  cutoff=final.model$eval.data %>% 
    slice_min(abs(Precision-ps)) %>%
    pull(Probability) %>% unique() %>% .[1]
  # based on this value, the data will be cut
  cut_data <- final.model$data %>% filter(score>=cutoff)
  
  # 12. The cut_data is prepared to be used for global network building. Let's first build the initial network:
  network.built <- X.build.complexes(data=cut_data, algo="SP",scores.col="score",init.stats=TRUE,final.stats=TRUE,
                                       standard.set=GS,labels.col="complex",  SP.finalpreds="original", rep.steps=20)
  
  # 13. The initial network is then refined.
  network.refined <- X.refine.complexes(data=network.built$data, algo="SCBP",final.stats=TRUE,
                                        standard.set=GS, rep.steps=20)
  
  # 14. Finally, the initial network is post-processed.
  network.post.split <- X.postprocess(data=network.refined$data, mode="split", scores.col="score", final.stats=TRUE, 
                                   standard.set=GS,labels="complex", weighted=FALSE)
  
  network.post.trim <- X.postprocess(data=network.post.split$data, mode="trim", scores.col="score", final.stats=TRUE, 
                                     standard.set=GS,labels="complex", weighted=TRUE)
  
  final_networks[[paste0("ps",ps)]] <- network.post.trim$data
  
  # Let's save the resulting plots and the network statistics into lists.
  
  statsplots[[paste0("ps",ps)]] <- 
    network.built$stats_initial$plot + ggtitle("Averaged") +
    network.built$stats_final$plot + ggtitle("Built") +
    network.refined$stats_final$plot + ggtitle("Refined") +
    network.post.split$stats_final$plot + ggtitle("Split") +
    network.post.trim$stats_final$plot + ggtitle("Trimmed") +
    patchwork::plot_layout(ncol=5,nrow=1)
  
  statsdata[[paste0("ps",ps)]] <- network.post.trim$stats_final$data %>% bind_cols(network.post.trim$stats_final$network.stats) %>% mutate(PS=ps) %>%
    mutate()
}

# 15. Now let's calculate a final network score to choose the final cutoff

choicedata <- statsdata %>%
  purrr::reduce(bind_rows) %>%
  mutate(network.score=stringency + fragmentation + purity + retrieval + recall/2)
chosen_precision <- choicedata %>% slice_max(network.score) %>% pull(PS)


choiceplot <- choicedata %>%
  ggplot(aes(x=PS,y=network.score)) +
  geom_point() + geom_line() + theme_bw()
choiceplot

# Using the final cutoff, we can plot the data:

# Save cut data as final_data
final_data <- final_networks[[paste0("ps",chosen_precision)]]
X.network.stats(final_data,GS,complex.stats=TRUE)

# Set up design parameters for the cytoscape network
design=list('setNodeShapeDefault("ELLIPSE")',
            'setNodeColorDefault("#5A5A5A")',
            'lockNodeDimensions(TRUE)',
            'setNodeLabelOpacityDefault(1)',
            'setEdgeColorMapping(table.column="complex",list("1","0"), list("#00FF00","#FF0000"),mapping.type="d")',
            'setEdgeLineWidthMapping(mapping.type="c",table.column="score", widths=c(1,20))'
)

# 16. Plot the cytoscape network using X.plot.network
X.plot.network(final_data,scores.col="score",standard.set=GS,labels.col="complex",annotation=plasmoDB_data%>%rename(protein=Accession), 
               design.params=design, network.name="MAPX_network", cytoscape.path="C://Program Files/Cytoscape_v3.9.1/Cytoscape.exe")


# Besides plotting of global protein-protein interaction maps, the MAP-X data can be used to plot local protein networks.
# In this approach, there is a specific protein of interests and our research question is what the interaction partners of this
# protein are. This approach is an antipole to the global approach where we do not ask about a specific protein but want to profile
# a global map of protein-protein interactions.

# 17. In local approaches, we work with probabilities (i.e., what is the probability that the protein of interest interacts with this other protein?)
# The model scores of some models (like classification trees) do not approximate the probabilities well and need to be calibrate.

calibrated.model <- X.calibrate.scores(final.model$data, GS)

# 18. Using the calibrated model, we can now find interactors of a specific protein of interest. Let's say, we are interested in a
# conserved protein of unknown function PF3D7_0412200.

local.subnetwork <- X.local.network(calibrated.model$data%>%dplyr::select(!score),"PF3D7_0412200",
                                    plot="cytoscape",
                                    plot.annotation=plasmoDB_data%>%rename(protein=Accession),
                                    plot.cytoscape.path="C://Program Files/Cytoscape_v3.9.1/Cytoscape.exe"
                                    )
local.subnetwork <- X.local.network(calibrated.model$data%>%dplyr::select(!score),"PF3D7_1022000",0.2,min.conditions=1,
                                    plot=NULL,
                                    plot.annotation=plasmoDB_data%>%rename(protein=Accession),
                                    plot.cytoscape.path="C://Program Files/Cytoscape_v3.9.1/Cytoscape.exe"
)

# 19. For several applications, it does not make sense to spend time and resource to train a model for each sample or replicate: for comparison 
# of assemblies of known protein complexes or for prediction of moonlighting proteins. In these cases, we can use the features extracted
# for each protein and reduce them into a PCA plot:

all.features.data <- data.frame()
for(tp in timepoints) {
  for(rep in reps) {
    all.features.data <- bind_rows(all.features.data,
                                   all.features[[tp]][[rep]] %>% mutate(replicate=as.character(rep)) %>% mutate(condition=tp)
    )
  }
}
data.PCAs <- X.PCA.reduction(all.features.data,plot=TRUE,plot.complexes=complexes)

# 20. The spread of the subunits on PCA plots nicely visualises the assembly state of the protein complexes. To put a number on it, we can
# calculate Euclidean distance from the centroid of the points of a protein complex and compare it to Euclidean distances of randomly selected proteins.
# Then, the program compares the random and complex-related Euclidean distances using Z-statistics to calculate the assembly factor.

# If multiple replicates are available, their averaging can be weighted by the data quality. In this case, we use the R2 goodness of fit
# that is calculated during curve fitting. The quality.data has two columns, protein and quality.

all.fitdata.frame <- data.frame()
for(tp in timepoints) {
  for(rep in reps) {
    all.fitdata.frame <- bind_rows(all.fitdata.frame,
                                   all.fitdata[[tp]][[rep]] %>% mutate(replicate=as.character(rep)) %>% mutate(condition=tp)
    )
  }
}
quality.data <- all.fitdata.frame %>%
  dplyr::select(protein,condition,replicate,R2) %>%
  rename(quality=R2)

# The function X.AI will calculate the assembly factors
data.AIs <- X.AI(data.PCAs$data,complexes, quality.data)

# 21. PCA data can also be used to predict moon-lighting subunits of known protein complexes.

data.moonlighting <- X.predict.moonlighters(data=data.PCAs$data,complexes=table1_complexes,min.reps=2,min.conditions=2,weights=data.PCAs$var.explained,
                                            p.adj.method="BH")

trialML <- X.predict.moonlighters(data=pcatrial$data,complexes=list(complex=c("PF3D7_1312600","PF3D7_0504600","PF3D7_0303700","PF3D7_1232200")),
                                  min.reps=1, min.conditions=2, weights=pcatrial$var.explained)

trialPCA <- X.PCA.reduction(all.features.data,plot=FALSE, average=TRUE)
trialAF <- X.AI(trialPCA$data,complexes, quality.data%>%rename(q=quality)%>%group_by(protein,condition)%>%dplyr::summarise(quality=mean(q))%>%
                  mutate(replicate="averaged"), animate=TRUE)

animation_to_save <- trialPCA$plots$complex +
  transition_states(condition) +
  labs(title = "{closest_state}hpi")

anim_save(paste0("animated_complex.gif"), animation = animation_to_save, height=9, width=9, units="cm",res=300,
          renderer = gifski_renderer())

animate(animation_to_save, fps = 60, renderer = av_renderer("animated_comple.mkv"))

anim_save(paste0("animated_complex.mvk"), animation = animation_to_save, height=9, width=9, units="cm",res=300,
          renderer = av_renderer())
