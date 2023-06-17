library(ranger)
library(treeshap)
library(dplyr)
#library(doParallel)


# StroblData = function(n = 120, 
#                       relevance = 0.15){
#   x2 = sample(2,n,replace=T)
#   xy = data.frame(x1 = rnorm(n),
#                     x2 = x2,
#                     x3 = sample(4,n,replace=T),
#                     x4 = sample(10,n,replace=T),
#                     x5 = sample(20,n,replace=T),
#                   y = rbinom(n,1,0.5+(c(-1,1)[x2]*relevance)))
#   return(xy)
# }

#wrapper function
simulation <- function(n = 120, M = 100, 
        relevance = 0.15, ntree = 100, max_depth = 50){
  
}
  
simulation_in_oob <- function(relevance = 0.15, n = 120, M = 100, 
           ntree = 500, max_depth = 50){
  time = Sys.time()
  shap_vals = shap_avs = shap_avs_oob = shap_avs_in= shap_vals_oob = shap_vals_in = data.frame(x1 = NA,x2 = NA,x3 =NA,x4=NA,x5=NA)
  slopes = R2 = matrix(0,nrow = M,ncol = ncol(shap_vals))
  for (i in 1:M){
    data = StroblData(n, relevance, yFactor = FALSE)
    suppressWarnings(rf_model <- ranger(
      formula = y ~ ., 
      data = data,
      num.trees = ntree,
      importance = "impurity",
      #mtry = 1,
      keep.inbag = T,
      max.depth = max_depth,
      classification = T
    ))
    #MDI as reference metric
    impurity_sim = rf_model$variable.importance
    unified = ranger.unify(rf_model,data)
    #inbag/oob treeshap calculation for smoothing
    zz <- file("treeshap.log", open = "wt")
    sink(zz)
    treeshap_obj = treewise_shap(unified,data,rf_inbag = rf_model$inbag.counts)
    shap_vals_in_sim = treeshap_obj[[1]]$shaps[,-6]
    shap_vals_oob_sim = treeshap_obj[[2]]$shaps[,-6]
    shap_avs_oob_sim <- shap_vals_oob_sim %>% abs() %>% colMeans() 
    shap_avs_in_sim <- shap_vals_in_sim %>% abs() %>% colMeans() 
    ###Raw treeshap calculation with original treeshap::treeshap function
    shap_obj = treeshap(unified,data)
    sink()
    shap_vals_sim = shap_obj$shaps[,-6]
    shap_avs_sim = shap_vals_sim %>% abs() %>% colMeans()
    shap_vals = rbind(shap_vals,shap_vals_sim)
    shap_avs = rbind(shap_avs,shap_avs_sim)
    for (j in 1:ncol(slopes)){
      fit = lm(shap_vals_oob_sim[,j] ~ shap_vals_in_sim[,j]-1)
      slopes[i,j] = fit$coefficients
      R2[i,j] =summary(fit)$r.squared
      
    }
    shap_vals_oob <- rbind(shap_vals_oob,shap_vals_oob_sim)
    shap_vals_in <- rbind(shap_vals_in,shap_vals_in_sim)
    shap_avs_oob <- rbind(shap_avs_oob,shap_avs_oob_sim)
    shap_avs_in <- rbind(shap_avs_in,shap_avs_in_sim)
  }
  cat("in simulation_in_oob, relevance=", relevance)
  print(Sys.time()-time)
  return(list(shap_vals_oob=shap_vals_oob[-1,],
              shap_vals_in=shap_vals_in[-1,],
              shap_avs_oob=shap_avs_oob[-1,],
              shap_avs_in=shap_avs_in[-1,],
              shap_vals=shap_vals[-1,],
              shap_avs=shap_avs[-1,],
              slopes=slopes, 
              R2=R2))
}


simulation_traintest <- function(relevance = 0.15, n = 120, M = 100, 
           ntree = 500, max_depth = 50){
  time = Sys.time()
  shaps = shaps_out = shaps_in = shap_avs = shap_avs_out = shap_avs_in = data.frame(x1 = NA,x2 = NA,x3 =NA,x4=NA,x5=NA)
  slopes = R2 = matrix(0,nrow=M,ncol=ncol(shaps))
  for (i in 1:M){
    data = StroblData(n, relevance, yFactor = FALSE)
    #split_temp = sample(1:n,round(0.5*n))
    #split = 1:n %in% split_temp
    split = 1:round(n/2)
    split2 = (round(n/2)+1):n
    #print(sum(split))
    train = data[split,]
    test = data[split2,]
    rf_train <- ranger(
      formula = y ~ ., 
      data = train,
      num.trees = ntree,
      max.depth = max_depth,
      classification = T
    )
    rf_test <- ranger(
      formula = y ~ ., 
      data = test,
      num.trees = ntree,
      max.depth = max_depth,
      classification = T
    )
    rf_model <- ranger(
      formula = y ~ ., 
      data = data,
      num.trees = ntree,
      max.depth = max_depth,
      classification = T
    )
    unif_train = ranger.unify(rf_train,train)
    unif_test = ranger.unify(rf_test,test)
    unif_model = ranger.unify(rf_model,data)
    zz <- file("treeshap.log", open = "wt")
    sink(zz)
    shap_vals_in_sim = treewise_shap(unif_test,test)$shaps[,-6]
    shap_vals_out_sim = treewise_shap(unif_train,test)$shaps[,-6]
    shap_vals_sim = treewise_shap(unif_model,data)$shaps[,-6]
    sink()
    
    for (j in 1:5) {
      fit = lm(shap_vals_out_sim[,j] ~ shap_vals_in_sim[,j]-1)
      slopes[i,j] = fit$coefficients
      R2[i,j] =summary(fit)$r.squared
    }
    shap_avs_in_sim <- shap_vals_in_sim %>% abs() %>% colMeans() 
    shap_avs_out_sim <- shap_vals_out_sim %>% abs() %>% colMeans()
    shap_avs_sim <- shap_vals_sim %>% abs() %>% colMeans()
    shaps_in = rbind(shaps_in,shap_vals_in_sim)
    shaps_out = rbind(shaps_out,shap_vals_out_sim)
    shaps = rbind(shaps,shap_vals_sim)
    shap_avs_in <- rbind(shap_avs_in,shap_avs_in_sim)
    shap_avs_out <- rbind(shap_avs_out,shap_avs_out_sim)
    shap_avs <- rbind(shap_avs,shap_avs_sim)
  }
  cat("in simulation_traintest, relevance=", relevance)
  print(Sys.time()-time)
  return(list(shaps_out=shaps_out[-1,],
              shaps_in=shaps_in[-1,],
              shap_avs_out=shap_avs_out[-1,],
              shap_avs_in=shap_avs_in[-1,],
              shaps=shaps[-1,],
              shap_avs=shap_avs[-1,],
              slopes=slopes, 
              R2=R2))
  #return(list(shaps_out[-1,],shaps_in[-1,],shap_avs_out[-1,],shap_avs_in[-1,],
   #           shaps[-1,],shap_avs[-1,],slopes, R2))
}

plot_figures <- function(res_list,labels,title){
  shap_avs_out <- res_list[[3]]
  shap_avs_in <- res_list[[4]]
  shap_avs <- res_list[[6]]
  slopes <- res_list[[7]]
  shap_avs_out_hat = shap_avs_out * slopes
  par(mar=c(2,2,8,2)+0.2,mfrow=c(1,4))
  boxplot(shap_avs, col = "darkorange", main = labels[1],ylab = "SHAP");grid()
  boxplot(shap_avs_in, col = "darkred",main = labels[2]);grid()
  boxplot(shap_avs_out, col = "darkblue",main = labels[3]);grid()
  boxplot(shap_avs_out_hat, col = "darkgreen",main = labels[4]);grid()
  mtext(title, outer = T,cex = 1.5, line = -2)
  return(recordPlot())
}

all_sims <- function(n,M,seed=NULL){
  set.seed(seed)
  titles = c("Null Simulation train/test","Power Simulation train/test","Null Simulation Inbag/OOB", "Power Simulation Inbag/OOB")
  path_names = c("null_sim_traintest","power_sim_traintest","null_sim_oob","power_sim_oob")
  for (i in 1:4){
    relevance = ifelse(i%%2,0,0.15)
    if (i<3){
      res_list = simulation_traintest(n,M,relevance)
      labels = c("all data","test in","test out","test_out smoothed") 
    }
    else{
      res_list = simulation_in_oob(n,M,relevance)
      labels = c("raw","inbag","oob","oob smoothed")
    }
    #fname = paste0(Sys.Date(),"-",path_names[i],".rda")
    save(res_list,file=file.path("../data/", paste0(Sys.Date(),"-", path_names[i],".rda",sep="")))
    png(file.path("../figures",paste0(Sys.Date(),"-",path_names[i],".png")))
    plot_figures(res_list,labels,titles[i])
    dev.off()
  }
}

#all_sims(120,100,seed=123)

