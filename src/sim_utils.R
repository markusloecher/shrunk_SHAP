#'Generate indexing matrix for inbag data of a ranger object
#'
#'@description 
#'get_inbag takes the ranger.inbag.counts object and transforms it into an indexing matrix
#'
#'@details 
#'This is a helper function. inbag.counts in the ranger package is giving the frequency of 
#'an observation as an element, so on position (1,1) we will get the frequency of the first 
#'observation in the first decision tree of the forest.
#'This function will transform this to get a matrix with which we can index the dataframe to retrieve
#'the inbag observations with their respective frequency.
#'
#'@param inbag this a ranger.inbag.counts object so a list of length num.trees where each element is a vector
#'of length nrow(data). The n'th entry will give the frequency of the n'th observation in the bootstrap of the
#'respective decision tree. Alternatively one can use an already transformed inbag.counts object where
#'the tranformation in the first if statement is already done
#'
#'@returns matrix of indices with dimensions nrow(data) times num.trees. Each column represents 
#'one decision tree and indexing the data with a column 
#'returns the bootstrap used for training the respective decision tree
get_inbag <- function(inbag){
  if(class(inbag) == "list"){
    inbag <- as.data.frame(t(as.matrix(as.data.frame(inbag))))
  }
  names(inbag) = 1:ncol(inbag)
  new_inbag = sapply(1:nrow(inbag),FUN = function(x){
    as.numeric(rep(names(inbag),inbag[x,]))})
  return(as.matrix(new_inbag))
}

#'Build a treeshap object for each inbag and oob data
#'
#'@description 
#'Built on the treeshap library and its c++ code for shap value calculation. This function generates shap
#'values for each individual decision tree for its inbag and out-of-bag observations. 
#'
#'@details 
#'The treeshap::treeshap function does not support treewise shap value generation. This function fills that gap.
#'Since the supporting c++ code of the treeshap function recurses over the forest roots, we can loop over each root
#'individually. This is what this function does. There is the option via returnTrees to either retrieve 
#'the shap values for each observation and tree, otherwise the shap values are summed over all trees. This
#'would be the equivalent of the original treeshap function but split into inbag and out-of-bag SHAP values.
#'
#'@param unified_model A helper object generated from treeshap::ranger_unify which builds a unified_model object 
#'from a ranger model
#'
#'@param x the data.frame which the ranger model used as an input. Treeshap does not at the time of building
#'this function support categorical data, therefore the data.frame needs to contain numeric data only.
#'
#'@param rf_inbag data which contains the information of how the bootstraps for each decision tree looked like
#'This can and should be the inbag.counts object from the ranger model. If None is given, the function just
#'calculates normal treeshap SHAP values. Default is `None`
#'
#'@param returnTrees boolean parameter. If true, returns SHAP values for each observation for each tree.
#'Otherwise sums up the SHAP values over all trees for each observation. Default is `FALSE`
#'
#'@returns if rf_inbag is None, returns a treeshap object, if rf_inbag is given, returns a list of 
#'two treeshap objects. One for inbag and one for out-of-bag SHAP values. The treeshap object is identical
#'to the one returned by the treeshap::treeshap function.
treewise_shap <- function(unified_model,x,rf_inbag = NULL, returnTrees=FALSE){
  if (is.null(rf_inbag)){
    return(treeshap(unified_model,x))
  }
  else{
    inbag_cnt <<- t(as.matrix(as.data.frame(rf_inbag)))
    #Raphael: make matrix of oob observations. Each column is a tree, rows are observations
    oob = t(inbag_cnt) == 0
    inbag = get_inbag(as.data.frame(inbag_cnt))
    model <- unified_model$model
    #Rewriting the unified_model columns such that c++ can use them, copied from treeshap function
    roots <- which(model$Node == 0) - 1
    nTrees=length(roots)
    yes <- model$Yes - 1
    no <- model$No - 1
    missing <- model$Missing - 1
    feature <- match(model$Feature, colnames(x)) - 1
    split <- model$Split
    decision_type <- unclass(model$Decision.type)
    is_leaf <- is.na(model$Feature)
    value <- model$Prediction
    verbose <- F
    #make matrix of correct dimensions to fill in shap values
    #if (returnTrees) {
      shap_vals_in = shap_vals_oob = array(NA,dim=c(dim(x),nTrees))
      dimnames(shap_vals_in) <- list(rownames(x), colnames(x),1:nTrees)
      dimnames(shap_vals_oob) <- list(rownames(x), colnames(x),1:nTrees)
    # } else {
    #   shap_vals_in = shap_vals_oob = matrix(0,nrow = dim(x)[1],ncol = dim(x)[2])
    #   dimnames(shap_vals_in) <- list(rownames(x), colnames(x))
    #   dimnames(shap_vals_oob) <- list(rownames(x), colnames(x))
    # }
    #calculate shap value for each individual tree and add them up
    for (i in 1:nTrees){
      #k is either 1 or 2 and switches between the inbag and oob case
      k = 0
      for(index in list(oob[,i],inbag[,i])){
        k = k + 1
        #x2 is the data subsetted by one tree in our indexing matrix
        x2 <- as.data.frame(t(as.matrix(x[index,])))
        #is_na and cover are relevant for the c++ code and depend on x2 so they are only calculated now
        is_na <- is.na(x2)
        cover<-  set_reference_dataset(unified_model,x[index,])$model$Cover
        #Calls c++ code but only for one root and as it is recursive, this works
        shap_vals_step = .Call('_treeshap_treeshap_cpp', PACKAGE = 'treeshap', x2, is_na, roots[i], yes, no, missing, feature, split, decision_type, is_leaf, value, cover, verbose)
        #add result for the correct observations to shap vals (addition not mean)
        if (k == 1){
          oob_index <- which(index)
          #if (returnTrees) 
          shap_vals_oob[oob_index,,i] <- shap_vals_step
          #else shap_vals_oob[oob_index,] <- shap_vals_oob[oob_index,] + shap_vals_step
        }
        else{
          #if (returnTrees) 
          shap_vals_in[index,,i] <- shap_vals_step
          #else shap_vals_in[index,] <- shap_vals_in[index,] + shap_vals_step
        }
      }
    }
    #only bookkeeping, so the output format is the same as for original treeshap
    interactions_array <- NULL
    if (returnTrees) {
      treeshap_obj_in <- list(shaps = shap_vals_in, interactions = interactions_array, 
                              unified_model = unified_model, observations = x)
      treeshap_obj_oob <- list(shaps = shap_vals_oob, interactions = interactions_array, 
                               unified_model = unified_model, observations = x)
    } else {
      #browser()
      #dim(shap_vals_in)
      shap_vals_in=apply(shap_vals_in, c(1,2), mean, na.rm=TRUE)*nTrees
      shap_vals_oob=apply(shap_vals_oob, c(1,2), mean, na.rm=TRUE)*nTrees
      treeshap_obj_in <- list(shaps = as.data.frame(shap_vals_in), interactions = interactions_array, 
                              unified_model = unified_model, observations = x)
      treeshap_obj_oob <- list(shaps = as.data.frame(shap_vals_oob), interactions = interactions_array, 
                               unified_model = unified_model, observations = x)
    }
    
    class(treeshap_obj_in) <- "treeshap"
    
    class(treeshap_obj_oob) <- "treeshap"
    #instead of one treeshap object, we get a list of two. One for the inbag data and one for the oob data
    return(list(inbag=treeshap_obj_in,oob=treeshap_obj_oob))
  }
}

SHAPviolin= function(res_wide, r=0.15, main = "test/test", 
                     ylabel = "SHAP",
                     YLIM=NULL,
                     my_geom = geom_violin(),
                     yAxis = "I") {
  main=paste0(main, "(r=",r,")" )
  
  colnames(res_wide) = paste0("x",1:5)
  
  res_long <- tidyr::pivot_longer(res_wide,
                                  x1:x5,
                                  names_to = "predictor",
                                  values_to = "SHAP"
  )
  p1 = ggplot(res_long, aes(x = predictor, y = SHAP)) 
  if (yAxis == "sqrt") p1 = p1 + scale_y_continuous(trans="S_sqrt")
  if (!is.null(YLIM)) p1 = p1 + coord_cartesian(ylim=YLIM) #ylim(YLIM[1],YLIM[2])
  p1 = p1 + my_geom + ggtitle(main) + ylab(ylabel)
  p1=p1+
    stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = .75,col=3)
  return(p1)
  #par(mfrow=c(2,2))
  #for (j in 1:5) hist(res$t1[,j], main = paste("X",j),xlab="t-score")
}

# j=1
# summary(lm(shap_vals_oob[,j] ~ shap_vals_in[,j]-1))$coefficients

slopes_generation_M <- function(shap_vals_in,shap_vals_oob,M = 1){
  n = nrow(shap_vals_in)/M
  slopes = matrix(0,nrow = M,ncol = ncol(shap_vals_oob))
  for (j in 1:ncol(slopes)){
    fit = sapply(1:M,FUN = function(x) summary(lm(shap_vals_oob[((n*(x-1))+1):(n*x),j] ~ shap_vals_in[((n*(x-1))+1):(n*x),j]-1))$coefficients[,3])
    slopes[,j] = fit
  }
  return(slopes)
}

slopes_generation <- function(shap_vals_in,
                  shap_vals_oob,
                  retSmooth = FALSE){
  n = nrow(shap_vals_in)
  p = ncol(shap_vals_oob)
  slopes = rep(0, p)
  if (retSmooth) shap_vals_oob_smooth=shap_vals_oob#silly initialization
  
  for (j in 1:p){
    fit = lm(shap_vals_oob[,j] ~ shap_vals_in[,j]-1)
    slopes[j] = as.numeric(fit$coefficients)
    if (is.na(slopes[j])) {
      #browser()
      slopes[j] = 1
    }
    #cat("j=", j, sum(is.na(shap_vals_in[,j])), sum(is.na(shap_vals_oob[,j])), "\n")
    if (retSmooth) shap_vals_oob_smooth[,j] = slopes[j]*shap_vals_in[,j]
  }
  if (retSmooth) return(shap_vals_oob_smooth)
  return(slopes)
}
