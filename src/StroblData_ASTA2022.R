
library(xgboost)
library(ranger)
library(treeshap)

library(future.apply)
#plan(multisession, workers = 1)#still waiting for access to our server, sigh
library(ggplot2)
library(ggpubr)
#define "signed" log/sqrt axes:
library(scales)
S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
IS_sqrt <- function(x){x^2*sign(x)}
S_sqrt_trans <- function() trans_new("S_sqrt",S_sqrt,IS_sqrt)


SHAPcorrs_lm = function(r=0.15, M=100,n=400,verbose=1){
  
  beta1=beta2=cor1=cor2=matrix(0,nrow=M,ncol=5)
  t1=t2=matrix(0,nrow=M,ncol=5)
  SHAPimp_shrunk=SHAPimp=matrix(0,nrow=M,ncol=5)
  
  for (i in 1:M){
    train= StroblData(n=n,relevance=r, yFactor=FALSE)
    test = StroblData(n=n,relevance=r, yFactor=FALSE)
    
    lm_train=lm(y ~ .,data=train)
    lm_test=lm(y ~ ., data=test)
    
    #=predict()
   }
}    
    
SHAPcorrs = function(r=0.15, M=100,n=400,
       alg = c("rf", "xgb")[1], 
       genData = c("strobl", "MARS", "MARSadd")[1],
       p=5,#<< num of features
       k=0,#<< if k>0 forced cardinality of informative features
       ntree = 500, 
       max_depth = 50,
       verbose=1){
  
  beta1=beta1_s=beta2=cor1=cor2=matrix(0,nrow=M,ncol=p)
  t1=t1_s=t2=matrix(0,nrow=M,ncol=p)
  SHAPimp_shrunk=SHAPimp_shrunk_s=SHAPimp=matrix(0,nrow=M,ncol=p)
  
  for (i in 1:M){
    if (genData == "strobl"){
      train= StroblData(n=n,relevance=r, yFactor=FALSE)
      test = StroblData(n=n,relevance=r, yFactor=FALSE)
    } else if (genData == "MARS"){
      #in this case we treat r = SNR !
      snr=r
      dataList =sim.mars1(n,nval=n,sigma=sqrt(6.29620/snr),p=p,k=k)
      
      train= dataList$train
      test = dataList$test
    }
    
    X_train = as.matrix(train[,-which(colnames(train)=="y")])
    y_train = train$y
    X_test = as.matrix(test[,-which(colnames(test)=="y")])
    y_test = test$y
    
    if (alg == "rf") {
      rf_train=ranger(y ~ .,data=train,num.trees = ntree,
                max.depth = max_depth, classification = T)
      rf_test=ranger(y ~ .,data=test,num.trees = ntree,
                  max.depth = max_depth, classification = T)
      
      urf_train <- ranger.unify(rf_train, train)
      urf_test <- ranger.unify(rf_test, test)
    } else if (alg == "xgb") {
      param <- list(objective = "reg:squarederror", max_depth = 6)
      
      xgb_train <- xgboost::xgboost(X_train, label = y_train, params = param, nrounds = 200, verbose = 0)
      xgb_test <- xgboost::xgboost(X_test, label = y_test, params = param, nrounds = 200, verbose = 0)
      
      urf_train <- xgboost.unify(xgb_train, X_train)
      urf_test <- xgboost.unify(xgb_test, X_test)
      
    } else stopifnot(FALSE)
    
    ts_train <- treeshap(urf_train, X_train, verbose = FALSE)
    #rf_ts <- colMeans(abs(ts_train$shaps))
    
    #trained on test, predict on test:
    ts_test_in <- treeshap(urf_test, X_test, verbose = FALSE)
    #trained on train, predict on test:
    ts_test_out <- treeshap(urf_train, X_test, verbose = FALSE)
    
    ts_test_in_shrunk = ts_test_in_shrunk_s = ts_test_in$shaps
    out_shap_s = as.data.frame(scale(ts_test_out$shaps))
    in_shap_s = as.data.frame(scale(ts_test_in$shaps))
    
    for (j in 1:p) {
      fit1   = lm(ts_test_out$shaps[,j] ~ ts_test_in$shaps[,j] -1)
      fit1_s = lm(out_shap_s[,j] ~ in_shap_s[,j] - 1)
      
      if (verbose>2){
        plot(ts_test_out$shaps[,j] ~ ts_test_in$shaps[,j], pch=20,cex=0.75, 
             col=rgb(0,0,1,0.5),xlab="out", ylab="in", main =paste0("X",j));grid()
        abline(fit1,col="darkgreen")
        fit1a = lm(ts_test_out$shaps[,j] ~ ts_test_in$shaps[,j] )
        abline(fit1a,col="brown")
      }
      beta1[i,j] = as.numeric(fit1$coefficients)
      beta1_s[i,j] = as.numeric(fit1_s$coefficients)
      
      fit1Sum   = summary(fit1)
      fit1Sum_s = summary(fit1_s)
      
      t1[i,j]   = fit1Sum$coefficients[1,"t value"]
      t1_s[i,j] = fit1Sum_s$coefficients[1,"t value"]
      
      cor1[i,j] = cor(ts_test_in$shaps[,j],ts_test_out$shaps[,j])
      ts_test_in_shrunk[,j]   = fit1$fitted.values
      ts_test_in_shrunk_s[,j] = fit1_s$fitted.values
      #ts_test_in_shrunk[,j] = beta1[i,j] * ts_test_in$shaps[,j]
      #ts_test_in_shrunk_s[,j] = beta1_s[i,j] * in_shap_s[,j]
      
      fit2 = lm(ts_train$shaps[,j] ~ ts_test_in$shaps[,j] -1)
      beta2[i,j] = as.numeric(fit2$coefficients)
      fit2Sum=summary(fit2)
      t2[i,j] = fit2Sum$coefficients[1,"t value"]
      cor2[i,j] = cor(ts_train$shaps[,j],ts_test_in$shaps[,j])
    }
    
    #scaling factor:
    if (0){
      g = rowSums(ts_test_in$shaps)/rowSums(ts_test_in_shrunk)
      g_s = rowSums(ts_test_in$shaps)/rowSums(ts_test_in_shrunk_s)
      #yes, I know I should use sweep instead of a for loop!
      for (kk in 1:nrow(ts_test_in_shrunk)){
        ts_test_in_shrunk[kk,] = g[kk]*ts_test_in_shrunk[kk,]
        ts_test_in_shrunk_s[kk,] = g_s[kk]*ts_test_in_shrunk_s[kk,]
      }
    }
    #g = rowSums(ts_test_in$shaps)/rowSums(ts_test_in_shrunk)
    
    SHAPimp[i,] = colMeans(abs(ts_test_in$shaps))
    SHAPimp_shrunk[i,] = colMeans(abs(ts_test_in_shrunk))
    SHAPimp_shrunk_s[i,] = colMeans(abs(ts_test_in_shrunk_s))
    
    if (verbose) if (i %% 20 == 0) cat("done with iteration", i, "\n")
  }
  res = list(cor1, cor2, beta1, beta1_s,beta2, SHAPimp,SHAPimp_shrunk, SHAPimp_shrunk_s, t1,t1_s, t2, ts_test_in, ts_test_out)
  names(res) = c("cor1", "cor2", "beta1", "beta1_s", "beta2", "SHAPimp","SHAPimp_shrunk","SHAPimp_shrunk_s","t1","t1_s", "t2", "ts_test_in", "ts_test_out")
  return(res)
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
pltSHAPcorrs = function(in_out, r=0.15, type=1){
  # type==1 boxplots of the slopes beta1/beta2 for train/test
  # type=="beta" violin plots of the slopes beta1/beta2 for train/test
  # type==2 boxplots of the correlations cor1/cor2 for train/test
  # type==3 boxplots of the raw/shrunk SHAP importances 
  # type==3.5 boxplots of the raw/shrunk/conf.shrunk SHAP importances
  # type=="ScaledBeta" boxplots of just the shrunk (scaled regression) SHAP importances
  # type==4 violin plots of the t scores
  # type=="beta-SE" strange pots of standard errors
  # type=="CI-beta" violin plots of the CI95 limits of the slopes
  # type=="SE" violin plots of standard errors
# NOTE: 5 and 6 are special since they use 
# local SHAP scores from just one simulation!
  # type==5 base R boxplots of the raw/shrunk SHAP importances 
  # type==6 violin plots of the SHAP scores for train/test
  # type==7 violin plots of the SHAP importance scores for train/test  
  library(ggplot2)
  #define "signed" log/sqrt axes:
  library(scales)
  S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
  IS_sqrt <- function(x){x^2*sign(x)}
  S_sqrt_trans <- function() trans_new("S_sqrt",S_sqrt,IS_sqrt)
  
  
  res = in_out[[as.character(r)]]
  main=paste0(" (r=",r,")" )
  tscore1_wide = as.data.frame(res$t1)
  colnames(tscore1_wide) = paste0("x",1:5)
  #compute standard errors from t-scores
  SE1_wide=tscore1_wide
  for (j in 1:ncol(SE1_wide)) SE1_wide[,j] = res$beta1[,j]/tscore1_wide[,j]
  
  SE1_long <- tidyr::pivot_longer(SE1_wide,
                                  x1:x5,
                                  names_to = "predictor",
                                  values_to = "StdErr"
  )
  beta1_wide = as.data.frame(res$beta1)
  colnames(beta1_wide) = paste0("x",1:5)
  beta1_long <- tidyr::pivot_longer(beta1_wide,
                                    x1:x5,
                                    names_to = "predictor",
                                    values_to = "beta1"
  )
  
  
  #browser()
  if (type %in% c(1:3,5)) par(mfrow=c(1,2))
  
  if (type==1){
    boxplot(res$beta1, main =paste0("test/test (r=",r,")" ), ylab ="betas")
    boxplot(res$beta2, main =paste0("train/test (r=",r,")" ), ylab ="betas")
    print(colMeans(res$beta1))
    print(colMeans(res$beta2))
  }
  if (type==2){
    #boxplot(res$cor1, main =paste0("test/test (r=",r,")" ), ylab ="corrs")
    #boxplot(res$cor2, main =paste0("train/test (r=",r,")" ), ylab ="corrs")
    #print(colMeans(res$cor1))
    res_wide = as.data.frame(res$cor1)
    p1 = SHAPviolin(res_wide, r, main="test/test", ylabel="corr")
    res_wide = as.data.frame(res$cor2)
    p2 = SHAPviolin(res_wide, r, main="train/test", ylabel="corr")
    
    p12 = ggarrange(p1, p2, nrow=1)
    return(p12)
  }
  if (type=="beta"){
    #boxplot(res$cor1, main =paste0("test/test (r=",r,")" ), ylab ="corrs")
    #boxplot(res$cor2, main =paste0("train/test (r=",r,")" ), ylab ="corrs")
    #print(colMeans(res$cor1))
    res_wide = as.data.frame(res$beta1)
    p1 = SHAPviolin(res_wide, r, main="test/test", ylabel="beta1")
    res_wide = as.data.frame(res$beta2)
    p2 = SHAPviolin(res_wide, r, main="train/test", ylabel="beta2")
    
    p12 = ggarrange(p1, p2, nrow=1)
    return(p12)
  }
  if (type==3){
    #boxplot(res$SHAPimp_shrunk, main =paste0("test/test (r=",r,")" ), ylab ="SHAPimp_shrunk")
    #boxplot(res$SHAPimp, main =paste0("train/test (r=",r,")" ), ylab ="SHAPimp")
    #print(colMeans(res$cor1))
    res_wide1 = as.data.frame(res$SHAPimp)
    res_wide2 = as.data.frame(res$SHAPimp_shrunk)
    YLIM1=range(res_wide1);YLIM2=range(res_wide2)
    YL=range(c(YLIM1,YLIM2))
    #cat("YLIM=", YL, "\n")
    p1 = SHAPviolin(res_wide1, r, my_geom= geom_boxplot(outlier.size = 0.3), main="raw ", ylabel="SHAP",YLIM=YL)
    p2 = SHAPviolin(res_wide2, r, my_geom= geom_boxplot(outlier.size = 0.3), main="shrunk ", ylabel="SHAP",YLIM=YL)
    
    p12=ggarrange(p1, p2, nrow=1)
    return(p12)
  }
  if (type==3.5){
    #boxplot(res$SHAPimp_shrunk, main =paste0("test/test (r=",r,")" ), ylab ="SHAPimp_shrunk")
    #boxplot(res$SHAPimp, main =paste0("train/test (r=",r,")" ), ylab ="SHAPimp")
    #print(colMeans(res$cor1))
    res_wide1 = as.data.frame(res$SHAPimp)
    res_wide2 = as.data.frame(res$SHAPimp_shrunk)
    YLIM1=range(res_wide1);YLIM2=range(res_wide2)
    YL=range(c(YLIM1,YLIM2))
    #cat("YLIM=", YL, "\n")
    p1 = SHAPviolin(res_wide1, r, my_geom= geom_boxplot(outlier.size = 0.3), main="raw ", ylabel="SHAP",YLIM=YL)
    p2 = SHAPviolin(res_wide2, r, my_geom= geom_boxplot(outlier.size = 0.3), main="shrunk ", ylabel="SHAP",YLIM=YL)
    
    for (j in 1:5){
      m = as.numeric(abs(res$t1[,j])>6)
      res_wide2[,j] = m* res_wide2[,j]
    }
    p3 = SHAPviolin(res_wide2, r, my_geom= geom_boxplot(outlier.size = 0.3), main="conf.shrunk ", ylabel="SHAP",YLIM=YL)
    
    
    p13=ggarrange(p1, p2, p3, nrow=1)
    return(p13)
  }
  if (type=="ScaledBeta"){#standardized regression
    
    #res_wide1 = as.data.frame(res$SHAPimp)
    res_wide2 = as.data.frame(res$SHAPimp_shrunk_s)
    #YLIM1=range(res_wide1);YLIM2=range(res_wide2)
    YL=range(res_wide2)#range(c(YLIM1,YLIM2))
    #cat("YLIM=", YL, "\n")
    #p1 = SHAPviolin(res_wide1, r, my_geom= geom_boxplot(outlier.size = 0.3), main="raw ", ylabel="SHAP",YLIM=YL)
    p2 = SHAPviolin(res_wide2, r, my_geom= geom_boxplot(outlier.size = 0.3), main="shrunk ", ylabel="SHAP",YLIM=YL)
    
    #p12=ggarrange(p1, p2, nrow=1)
    return(p2)
  }
  if (type==4){
    res_wide = as.data.frame(res$t1)
    colnames(res_wide) = paste0("x",1:5)
    #b1=boxplot(res_wide, main =paste0("test/test (r=",r,")" ), ylab = "t-score")
    main=paste0(" (r=",r,")" )
    
    res_long <- tidyr::pivot_longer(res_wide,
      x1:x5,
      names_to = "predictor",
      values_to = "tScore"
    )
    p1 = ggplot(res_long, aes(x = predictor, y = tScore)) + #scale_y_sqrt()
      scale_y_continuous(trans="S_sqrt")#,breaks=seq(-0.1,0.5,0.05
    #p1 = p1 +geom_boxplot(outlier.size = 0.3) 
    p1 = p1 + geom_violin() + ggtitle(main)
    p1 = p1 + geom_hline(yintercept = c(-1,1)*sqrt(2.5), linetype="dashed", color = "red", size=.5) 
    p1=p1+
      stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = .75,col=3)
    return(p1)
    #par(mfrow=c(2,2))
    #for (j in 1:5) hist(res$t1[,j], main = paste("X",j),xlab="t-score")
  }
  if (type=="beta-SE"){
    main=paste0(" (r=",r,")" )
    
    beta1_SE1_long <- cbind.data.frame(beta1_long, SE1_long)[,-3] #merge(beta1_long, SE1_long, by.x="predictor")
    
    p1 = ggplot(beta1_SE1_long, aes(x = beta1, y = StdErr)) + geom_point() #scale_y_sqrt()
    p1 = p1 + facet_wrap(~predictor,ncol=1, scales = "free_y") + ggtitle(main)
    
    return(p1)
  }
  
  if (type=="CI-beta"){
    
    CIl95_beta1_wide = beta1_wide
    for (j in 1:5){
      CIupper = beta1_wide[,j] + 1.96*SE1_wide[,j]
      CIlower = beta1_wide[,j] - 1.96*SE1_wide[,j]
      
      mixedSign = CIupper*CIlower<0
      CIl95_beta1_wide[mixedSign,j] = 0 #mixed sign means that 0 is contained
      allNegative = CIupper<0 & CIlower <0
      CIl95_beta1_wide[allNegative,j] = CIupper[allNegative] #take the upper 
      allPositive = CIupper>0 & CIlower >0
      CIl95_beta1_wide[allPositive,j] = CIlower[allPositive] #take the lower
      if (!all(mixedSign | allNegative | allPositive)) print("not a partition!")
    }
    CIl95_beta1_long <- tidyr::pivot_longer(CIl95_beta1_wide,
                                      x1:x5,
                                      names_to = "predictor",
                                      values_to = "beta1l95"
    )
    p1 = ggplot(CIl95_beta1_long, aes(x = predictor, y = beta1l95)) 
    p1 = p1 + geom_violin() + ggtitle(main)
    
    p1=p1+
      stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = .75,col=3)
    return(p1)
    
  }
  if (type=="SE"){
    main=paste0(" (r=",r,")" )
    SE1_long <- tidyr::pivot_longer(SE1_wide,
                                    x1:x5,
                                    names_to = "predictor",
                                    values_to = "StdErr"
    )
    p1 = ggplot(SE1_long, aes(x = predictor, y = StdErr)) + #scale_y_sqrt()
      scale_y_continuous(trans="S_sqrt")#,breaks=seq(-0.1,0.5,0.05
    #p1 = p1 +geom_boxplot(outlier.size = 0.3) 
    p1 = p1 + geom_violin() + ggtitle(main)
    #p1 = p1 + geom_hline(yintercept = c(-1,1)*sqrt(2.5), linetype="dashed", color = "red", size=.5) 
    p1=p1+
      stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = .75,col=3)
    return(p1)
    #par(mfrow=c(2,2))
    #for (j in 1:5) hist(res$t1[,j], main = paste("X",j),xlab="t-score")
  }
  if (type==5){
    boxplot(res$ts_test_in$shaps, main =paste0("test/test (r=",r,")" ), ylab = "SHAP")
    boxplot(res$ts_test_out$shaps, main =paste0("train/test (r=",r,")" ), ylab = "SHAP")
  }
  
  if (type==6){
    res_wide = as.data.frame(res$ts_test_in$shaps)
    p1 = SHAPviolin(res_wide, r, "test/test")
    res_wide = as.data.frame(res$ts_test_out$shaps)
    p2 = SHAPviolin(res_wide, r, "train/test")
    
    p12=ggarrange(p1, p2, nrow=1)
    return(p12)
  }
  if (type==7){
    res_wide = as.data.frame(res$SHAPimp)
    p1 = SHAPviolin(res_wide, r, "raw ", ylabel="SHAP")
    res_wide = as.data.frame(res$SHAPimp_shrunk)
    p2 = SHAPviolin(res_wide, r, main="shrunk ", ylabel="SHAP")
    
    p12=ggarrange(p1, p2, nrow=1)
    return(p12)
  }
}

StroblData = function(
  n=120, # number of rows in data
  # number of simulations
  #nCores = M, # number of cores to use; set to 1 on Windows!
  relevance = 0.15, # signal strength (0 for NULL)
  data_type = 1, #1 for classification, 2 for regression
  yFactor = TRUE, #y is factor ?
  #ntree = 50, #number of trees in forest
  #correctBias = c(inbag=TRUE,outbag=TRUE),
  verbose=0
){
 
    x1 = rnorm(n)
    x2 = base::sample(1:2,n,replace=TRUE)
    x3 = base::sample(1:4,n,replace=TRUE)#rmultinom(n, 1, prob=rep(1/k,k))
    x4 = base::sample(1:10,n,replace=TRUE)
    x5 = base::sample(1:20,n,replace=TRUE)
    
    if (data_type==1){
      y= rbinom(n,1,p=0.5 + relevance*c(-1,1)[as.numeric(x2)])#NULL case
      if (yFactor) y=factor(y)
    #test: rbinom(8,1,p=0.5 + rep(c(-1,1),4)*0.4)
    } else if (data_type==2){
      #We call relevance SNR
      #SNR is defined as Var(F(x))/Var(noise)
      #Var(F(x)) = Var(beta2*x2) = beta2^2*0.25, Var(noise)=1
      #So relevance = beta2^2*0.25 -> beta2=sqrt(relevance)/0.5
      beta2 = 2*sqrt(relevance)
      if (verbose) cat("SNR=", relevance,  ", so beta2=", beta2, "\n")
      y = beta2*as.numeric(x2) + rnorm(n)
      
      if (verbose>2) boxplot(y ~ x2)
        #plot(x2,y,col=rgb(0,0,1,0.5),pch=20)
    }  
    
    data = cbind.data.frame(x1,x2,x3,x4,x5,y)
    colnames(data) = c("x1", "x2", "x3", "x4", "x5", "y")
      
    return(data)
}

TuneRF = function(SNR=1, n=200, mtry=c(2:5), 
                  min.node.size = 4*(1:12), max.depth = NULL, Navs = 10, verbose=0){
  #xy = StroblData(data_type=2, relevance=SNR)
  #rf=ranger(y ~ .,data=xy)
  #prediction.error	=Overall out of bag prediction error=MSE
  
  pars = expand.grid(mtry=mtry, min.node.size = min.node.size, 
                     max.depth = max.depth, n=n, SNR = SNR)
  K = nrow(pars)
  pars$OOBerr = NA
  
  for (i in 1:K){
    oob=0
    for (j in 1:Navs){
      xy = StroblData(data_type=2, relevance=pars[i,"SNR"],n = pars[i,"n"])
      rf=ranger(y ~ .,data=xy, min.node.size = pars[i,"min.node.size"], 
                mtry = pars[i,"mtry"], max.depth=pars[i,"max.depth"])
      oob=oob+rf$prediction.error
    }
    pars$OOBerr[i] = oob/Navs
    if (verbose) cat("done with", i, "\n")
  }
  #attr(pars,"min.node.size")=min.node.size
  #attr(pars,"mtry")=mtry
  return(pars)
}

StroblSHAP = function(SNR=1, n=200, mtry=2, max.depth = 10,
                  min.node.size = 5, onlyGini=FALSE, verbose=0){
  library(ranger)
  library(treeshap)
  
  
  xy = StroblData(data_type=2, relevance=SNR,n = n)
  if (verbose>1) {
    #tmp=rnorm(5)+ (1:5)
    #names(tmp) = paste0("x",1:5)
    #tmp <- c(max.depth=max.depth,mtry = mtry, tmp)
    #return(tmp)#just a dummy !
    mlm = lm(y ~ .,data=xy)
    mlm_sum = summary(mlm)
    tVals = mlm_sum$coefficients[-1,"t value"]
    
    return(c(SNR=SNR,tVals))
  }
  
  if (onlyGini){
    rf=ranger(y ~ .,data=xy, min.node.size = min.node.size, 
              mtry = mtry, max.depth=max.depth, importance="impurity")
    rf_ts = rf$variable.importance
  } else {
    rf=ranger(y ~ .,data=xy, min.node.size = min.node.size, 
              mtry = mtry, max.depth=max.depth)
    urf <- ranger.unify(rf, xy)
    ts <- treeshap(urf, xy[, -which(colnames(xy)=="y")], verbose = FALSE)
    rf_ts <- colMeans(abs(ts$shaps))
  }
  
  
  rf_ts <- c(max.depth=max.depth,mtry = mtry, rf_ts)
  return(rf_ts)
}

StroblSHAP_replicate = function(#compute SHAP values on Strobl data sweeping over various parameters
  SNR=1,
  max.depth.range = c(1:3,5,10,30),
  mtry.range = 1:5,
  n=200,
  Nreps=20,
  onlyGini=FALSE,
  verbose=1
){
  on.exit(save(res_SNR, file = "res_SNR_tmp.rda"))
  
  res_SNR=list();i=1
  for (mtry in mtry.range){
    for (max.depth in max.depth.range){
      res_SNR[[i]] <- t(future_replicate(Nreps, StroblSHAP(SNR=SNR, n=n, mtry=mtry, max.depth = max.depth,onlyGini=onlyGini, verbose=0), future.seed=TRUE))
      i=i+1
    }
    if(verbose) cat("done with mtry", mtry, "\n")
  }
  if(verbose>1) save(res_SNR, file = "res_SNR_tmp.rda")#to be safe
  
  res_SNR_wide = do.call("rbind.data.frame",res_SNR)
  
  res_SNR_long <- tidyr::pivot_longer(res_SNR_wide,
    x1:x5,
    names_to = "predictor",
    values_to = "SHAP"
  )
  res_SNR_long$max.depth = factor(res_SNR_long$max.depth)
  res_SNR_long$mtry = factor(res_SNR_long$mtry)  
  
  return(list(res_SNR=res_SNR,res_SNR_long=res_SNR_long))
}

plot_res_SNR= function(res_SNR,SNR,sy =scale_y_sqrt(), os = 0.2){
  library(ggplot2)
  res_SNR_long=res_SNR$res_SNR_long
  p1 = ggplot(res_SNR_long, aes(x = max.depth, y = SHAP, 
           fill = predictor)) + geom_boxplot(outlier.size = os) +
           facet_wrap(~ mtry) + geom_hline(yintercept = sqrt(SNR), linetype="dashed", color = "red", size=.5) 
  p1 + sy
}

plot_OOBerr = function(pars, xPar=c("min.node.size", "max.depth")[1],
                       BayesErr=1,l="",showLeg=FALSE, main = ""){
  mtry=sort(unique(pars$mtry))
  K=length(mtry)
  x=sort(unique(pars[,xPar]))
  mm=max(x)
  
  OOBerr=matrix(pars$OOBerr, ncol=K,byrow = TRUE)
  #browser()
  colnames(OOBerr) = paste("mtry",mtry,sep="")
  rownames(OOBerr) = paste(xPar,x,sep="")
  
  ylim=range(c(BayesErr,range(pars$OOBerr)))
  if (showLeg) xlim=c(1,mm*1.2) else xlim=c(1,mm)
  if (l=="x") xlim = log(xlim)
  
  matplot(x,OOBerr, type="l", lwd=2, xlab="",lty=1,
           ylim=ylim, log=l);grid()
  title(main);grid()
  abline(h=BayesErr, lty=2, col = "darkgreen",lwd=2)
  mtext(xPar, side=1,at=mm*0.9)
  if (showLeg) legend("topright", col = 1:K, legend=paste(c(mtry)), title="mtry",lwd=1.5)
}
if (0){
  library(dplyr)
  xy = StroblData(data_type=2,verbose=3)
  res=list()
  for (max.depth in c(1:3,5,7,10,30)){
    res[[as.character(max.depth)]] <- t(future_replicate(20, StroblSHAP(SNR=1, n=200, mtry=2, max.depth = max.depth, verbose=0), future.seed=TRUE))
    #xm=colMeans(res[[as.character(max.depth)]][,paste0("X",1:5)])
    #xs=apply(res[[as.character(max.depth)]][,paste0("X",1:5)],2,IQR)
    
  }
  
  res_wide = do.call("rbind.data.frame",res)
  #res_long <- tidyr::gather(res_wide, SHAP, x1:x5, factor_key=TRUE)
  
  res_long <- res_wide %>% tidyr::pivot_longer(
    x1:x5,
    names_to = "predictor",
    values_to = "SHAP"
  )
  res_long$max.depth = factor(res_long$max.depth)
  
  library(ggplot2)
  p1 = ggplot(res_long, aes(x = max.depth, y = SHAP, fill = predictor)) + 
    geom_boxplot(outlier.size = 0.2)
  p1
  #res <- future_sapply(rep(0.2, 5), StroblSHAP,future.seed=TRUE)
  #res <- future_lapply(rep(0.2, 5), StroblSHAP,n=200, mtry=2, max.depth = 10,verbose=3,future.seed=TRUE)
  #res=do.call("rbind.data.frame",res)
  
  #test:
  res_SNR004 = StroblSHAP_replicate(0.04,max.depth.range = c(1:3,5,10,30)[1:2],
                                    mtry.range = 1:2,
                                    n=200,
                                    Nreps=2,
                                    verbose=2)
}