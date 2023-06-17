AppendSims = function(in_out, in_out_new,verbose=1){
  stopifnot(all(names(in_out) == names(in_out_new)))
  
  for (r_str in names(in_out)){
    x_old = in_out[[r_str]]
    x_new = in_out_new[[r_str]]
    stopifnot(all(names(x_old) == names(x_new)))
    Metrics = names(x_old)
    Metrics = Metrics[!Metrics %in% c("ts_test_in", "ts_test_out") ]
    for (x_str in Metrics){
      if (verbose){
        cat("xold", x_str, dim(x_old[[x_str]]), "\n")
        cat("x_new", x_str, dim(x_new[[x_str]]), "\n")
      }
      in_out[[r_str]][[x_str]] = rbind.data.frame(x_old[[x_str]], x_new[[x_str]])
    }
  }
  return(in_out)
}

#test:
if (0){
  tmp = AppendSims(in_out, in_out)
}

cont_to_cat <- function(x, k, ranFun =c("norm","unif")[1] , ...) {
  if (k > 0) {
    if (ranFun == "norm") 
      breaks = qnorm(p = seq(0, 1, length.out = k+1), ...)
    if (ranFun == "unif") 
      breaks = qunif(p = seq(0, 1, length.out = k+1), ...)
    res <- cut(x, 
               breaks = breaks, 
               labels = LETTERS[1:k], 
               include.lowest = TRUE) 
    (as.numeric(res) - 1)/(k-1) #- .5
  } else {
    x
  }
  
}
#test:
if (0){
  tmp = cont_to_cat(rnorm(20), 2)
  table(tmp)
}
################Mentch Zhou paper:####################
enlist <- function (...) 
{
  result <- list(...)
  if ((nargs() == 1) & is.character(n <- result[[1]])) {
    result <- as.list(seq(n))
    names(result) <- n
    for (i in n) result[[i]] <- get(i)
  }
  else {
    n <- sys.call()
    n <- as.character(n)[-1]
    if (!is.null(n2 <- names(result))) {
      which <- n2 != ""
      n[which] <- n2[which]
    }
    names(result) <- n
  }
  result
}
# Simulating data from the MARS model, which contains 
# interaction terms.
sim.mars1 <- function(n, #<< num rows in the train data
                      nval,#<< num rows in the test data
                      sigma, #<< noise stdev
                      p = 10,#<< num of features
                      k = 0 #<< if k>0 forced cardinality of informative features
                     # s=5 #<< s <= p is the number of features with a nonzero coefficient thus considered signal
){
  stopifnot(p>=5)
  
  if (length(k)==1) k = rep(k,p)
  
  x <- matrix(runif(p*n,0,1),n,p)
  xval <- matrix(runif(p*nval,0,1),nval,p)
  for (j in 1:p){
    if (k[j]>0){
      x[,j] = cont_to_cat(x[,j],k[j], ranFun = "unif")
      xval[,j] = cont_to_cat(xval[,j],k[j], ranFun = "unif")
    }
  }
  e <- matrix(rnorm(n,0,sigma),n,1)
  eval <- matrix(rnorm(nval,0,sigma),nval,1)
  mu <- as.vector(10*sin(pi*x[,1]*x[,2])+20*(x[,3]-0.05)^2 + 10*x[,4]+5*x[,5])
  muval <- as.vector(10*sin(pi*xval[,1]*xval[,2])+20*(xval[,3]-0.05)^2 + 10*xval[,4]+5*xval[,5])
  y <- as.vector(10*sin(pi*x[,1]*x[,2])+20*(x[,3]-0.05)^2 + 10*x[,4]+5*x[,5]+e)
  yval <- as.vector(10*sin(pi*xval[,1]*xval[,2])+20*(xval[,3]-0.05)^2 + 10*xval[,4]+5*xval[,5]+eval)
  
  train = cbind.data.frame(x,y);colnames(train) = c(paste0("x",1:p),"y")
  test = cbind.data.frame(xval,yval);colnames(test) = c(paste0("x",1:p),"y")
  enlist(train,test)
  #enlist(x,y,xval,yval,mu,muval)
}
#test:
if (0){
  tmp = sim.mars1(n=200,nval=200,sigma=sqrt(6.29620))
}
# Simulating data from the MARSadd model.
sim.mars2 <- function(n,nval,sigma){
  x <- matrix(runif(10*n,0,1),n,10)
  xval <- matrix(runif(10*nval,0,1),nval,10)
  e <- matrix(rnorm(n,0,sigma),n,1)
  eval <- matrix(rnorm(nval,0,sigma),nval,1)
  mu <- as.vector(0.1*exp(4*x[,1])+4/(1+exp(-20*(x[,2]-0.5)))+
                    3*x[,3]+2*x[,4]+x[,5]        )
  muval <- as.vector(0.1*exp(4*xval[,1])+4/(1+exp(-20*(xval[,2]-0.5)))+
                       3*xval[,3]+2*xval[,4]+xval[,5]        )
  y <- as.vector(0.1*exp(4*x[,1])+4/(1+exp(-20*(x[,2]-0.5)))+
                   3*x[,3]+2*x[,4]+x[,5] + e          )
  yval <- as.vector(0.1*exp(4*xval[,1])+4/(1+exp(-20*(xval[,2]-0.5)))+
                      3*xval[,3]+2*xval[,4]+xval[,5] + eval         )
  
  enlist(x,y,xval,yval,mu,muval)
}

RunSim = function(
  n=120,#<< sample size
  M=100, #<< number of replication
  p=5,#<< num of features
  k=0,#<< if k>0 forced cardinality of informative features
  rVals = c(0.001, 0.005, 0.01, 0.02, 0.05,0.15),
  alg = c("rf", "xgb")[1], 
  genData = c("strobl", "MARS", "MARSadd")[1],
  height=5, 
  width=5,
  main ="MARS, k=2 for p=1:3",
  rerun = c(FALSE,TRUE)[2],
  addSims = c(FALSE,TRUE)[1],
  verbose 
  = 1
){
  kAll = paste0(k,collapse="")
  fileInfo = paste0("in_out_",alg,"_",genData,"_k",kAll,"_n",n,"_p",p)
  rdaFile = paste0("data/",fileInfo,".rda")
  
  if (verbose > 2) {
    print(fileInfo)
    return()
  }
  pars = list(n=n,M=M, alg = alg, genData = genData, k=k,p=p, desc="")
  
  if (rerun | addSims){
    #snr in 0.05, 0.09, 0.14, ..., 6.00
    in_out_new=list()
    for (r in rVals){
      cat("r=", r, "\n")
      in_out_new[[as.character(r)]] = SHAPcorrs(r,n=n, M=M, alg = alg, genData = genData, p=p, k=k)
    }
    
    if (addSims) {
      rdaFile_new = paste0("data/",fileInfo,"new.rda")
      save(in_out_new, file = rdaFile)
      load(rdaFile)
      in_out = AppendSims(in_out, in_out_new)
    } else {
      in_out = in_out_new
    }
    #in_out_xgb = in_out
    save(in_out, pars, file = rdaFile)
  } else {
    if (verbose) cat("opening", rdaFile, "\n")
    load(rdaFile)
  }
  ### plot SHAP values
  
  plotShrunkSHAP(in_out, main = main, log_r=TRUE)
  
  #fname1 = paste0("figures/", fileInfo,Sys.Date(),".pdf")
  fname2 = paste0("figures/", fileInfo,".pdf")
  ggsave(fname2, height=height, width=width)
  #system(paste("cp ", fname1, fname2))
  
  invisible(in_out)
}
#test:
if (0){
  tmp = RunSim(n=20,M=3,p=5,k=0,alg = "xgb",genData = "MARS")
}
