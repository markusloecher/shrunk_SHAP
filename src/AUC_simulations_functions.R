#The data has 1000 samples with 50 features. All
#features are discrete, with the j th feature containing j + 1 distinct values 0,1,...,j . We randomly
#select a set S of 5 features from the first ten as relevant features. The remaining features are noisy
#feature
library(foreach)
library(doParallel)
#library(rfVarImpOOB)
library(ranger)
library(treeshap)


NoisyFeatureSim = function(N=500, p=50, ntree=100, onlyX = FALSE, verbose=0){
  x=matrix(0,nrow=N,ncol=p)
  for (k in 1:p){
    x[,k] = sample(0:k,N,replace=TRUE)
  }
  
  if (onlyX) return(x)
  #rlvFtrs = sample(10,5) #relevant Features
  rlvFtrs = rep(0,p)
  rlvFtrs[sample(10,5)] =  1
  y=Xrlv=rep(0,N)
  for (k in which(rlvFtrs==1)) {
    Xrlv=Xrlv+x[,k]/k
  }
  #hist(2*Xrlv/5-1) # nicely symmetric around 0 !
  y = factor(rbinom(N,1,plogis(2*Xrlv/5 - 1)))
  data=cbind.data.frame(x,y)
  colnames(data)[1:p] = paste0("x",1:p)
  #ranger.unify does not take categorical values so we change everything to numeric
  data <- sapply(data, as.numeric)
  rfobj <- ranger(y ~ ., data = data, keep.inbag=TRUE,  importance = c("impurity"), num.trees =ntree, classification = T)
  if (verbose>1) browser()
  imp_MDI=rfobj$variable.importance
  #try({impRf = importance(rf)})
  rfobj <- ranger(y ~ ., data = data, keep.inbag=TRUE,  importance = c("impurity_corrected"), num.trees =ntree, classification =T)
  if (verbose>1) browser()
  imp_AIR=rfobj$variable.importance

  #non-smoothed shap values
  rf_unif <- ranger.unify(rfobj,data)
  treeshap_obj = treewise_shap(rf_unif,data,rf_inbag = rfobj$inbag.counts)
  shap_total = treeshap_obj[[2]]$shaps[,-(p+1)]
  shap = colMeans(abs(shap_total))
  
  #train/test splitting for smoothed shaps
  train_data = data[1:round(N/2),]
  test_data = data[(round(N/2)+1):N,]
  #slopes generation
  rf_train = ranger(y ~ ., data = train_data, num.trees =ntree, mtry = 3)
  rf_test = ranger(y ~ ., data = test_data, num.trees =ntree, mtry = 3)
  rf_train_unif = ranger.unify(rf_train,train_data)
  rf_test_unif = ranger.unify(rf_test,test_data)
  zz <- file("treeshap.log", open = "wt")
  sink(zz)
  #forward train/test split
  shap_in_f = treeshap(rf_test_unif,test_data)
  shap_out_f = treeshap(rf_train_unif,test_data)
  #backward train/test split, aka. train and test data swapped
  shap_in_b = treeshap(rf_train_unif,train_data)
  shap_out_b = treeshap(rf_test_unif,train_data)
  sink()
  #generating forward and backward slopes
  slopes_f <- slopes_generation(shap_in_f$shaps[,-(p+1)],shap_out_f$shaps[,-(p+1)])
  slopes_b <- slopes_generation(shap_in_b$shaps[,-(p+1)],shap_out_b$shaps[,-(p+1)])
  #mean of slopes
  slopes = (slopes_f+slopes_b)/2
  #smoothing of the "true" shap values
  shap_smooth = c(slopes)*t(shap_total)
  #browser()
  shap_smooth = rowMeans(abs(shap_smooth))

  return(cbind(MDI = imp_MDI, AIR = imp_AIR, shap = shap,shap_smooth = shap_smooth,rlvFtrs = rlvFtrs))
}


