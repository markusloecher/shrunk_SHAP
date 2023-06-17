#The data has 1000 samples with 50 features. All
#features are discrete, with the j th feature containing j + 1 distinct values 0,1,...,j . We randomly
#select a set S of 5 features from the first ten as relevant features. The remaining features are noisy
#feature
library(foreach)
library(doParallel)
library(ranger)

N=2000
p=50
ntree=500

Nsims =  100
ncores = 8
set.seed(NULL)
verbose=FALSE


source("src/AUC_simulations_functions.R")
source("src/sim_utils.R")

seriousRun=TRUE

if (seriousRun){
  t = Sys.time()
  if (ncores>1){
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
     parSims = foreach(i=1:Nsims, .combine=rbind, .packages=c("ranger","treeshap")) %dopar% {
      NoisyFeatureSim(N,p,ntree)
    }
    stopCluster(cl)
  } else {
    parSims = list()
    for (i in 1:Nsims) {
      print(i)
      parSims[[i]] <-  NoisyFeatureSim(N,p,ntree)
    }
    parSims = do.call("rbind", parSims)
  }
  print(Sys.time()-t)
} else {
  #single run:
  res = NoisyFeatureSim(N=50, p=50, ntree=10)
}
closeAllConnections()
save(parSims, file = paste0("AUCsims_",Sys.Date(),".rda"))

if (1){
  library(AUC)
  
  M=ncol(parSims)
  aucSims = matrix(0,nrow=1,ncol=M-1);k=1
  colnames(aucSims) = colnames(parSims)[1:(M-1)] #, "MDA", "MDI","Mod_GOOB_node1","Mod_GOOB_node5")
  options(digits=3)
  for (j in 1:(M-1)){
    aucSims[1,k] = auc(roc(parSims[,j],factor(parSims[,"rlvFtrs"])))
    cat(colnames(parSims)[j],aucSims[1,k], "\n")
    k=k+1
  }
  library(xtable)
  #print(xtable(aucSims),comment=FALSE)
}
