library('latex2exp')

plotinOutCorrs = function(test_in, test_out, cols2plot=c(1,2,4,5), 
                          xlab = expression(SHAP["in"]), 
                          ylab = expression(SHAP["out"]), nc=2, 
                          plotIt=c("base", "ggplot")[1]){
  
  p = length(cols2plot)
  nr = ceiling(p / nc)
  R2=r=rep(NA,p)
  
  if (plotIt=="base") {
    
    par(mfrow=c(nr,nc),cex=0.75, mar=c(3,4,3,1))
    for (j in cols2plot){
      try({
        fit = lm(test_out[,j] ~ test_in[,j] -1)
        R2[j] = round(summary(fit)$r.sq,2)
        r[j] = cor(fit$fitted.values, test_out[,j])
        #fit0 = lm(test_out[,j] ~ test_in[,j] -1)
        
        #if (plotIt) {
          #main title
          #mT = paste0(colnames(test_out)[j], " (R2 =", R2[j],")")
          mT = bquote(x[.(j)] ~ " (" ~ R^2 ~ " =" ~.(R2[j]) ~ ")"  )
          
          plot(test_in[,j],test_out[,j],col=rgb(0,0,1,0.5), pch=20,cex=0.75, 
               xlab = "", ylab="", main = mT);grid()
          title(ylab=ylab, line=2)#, cex.lab=1.2)
          title(xlab=xlab, line=2)#, cex.lab=1.2)#, family="Calibri Light")
          abline(fit,col=2, lwd=1.5)
          #abline(fit0,col="darkgreen")
        #}
      })
    }
  } else if (plotIt=="ggplot") {
    for (j in cols2plot){
      fit = lm(test_out[,j] ~ test_in[,j] -1)
      R2[j] = round(summary(fit)$r.sq,4)
      r[j] = cor(fit$fitted.values, test_out[,j])
    }
    print(R2)
    library(ggplot2)
      test_in_long <- tidyr::pivot_longer(test_in,
                       x1:x5,
                       names_to = "predictor",
                       values_to = "SHAP_in")
      test_out_long <- tidyr::pivot_longer(test_out,
                        x1:x5,
                        names_to = "predictor",
                        values_to = "SHAP_out")
      test_long <- cbind.data.frame(test_in_long, test_out_long)[,-3] #merge(beta1_long, SE1_long, by.x="predictor")
      test_long <- subset(test_long, predictor %in% paste0("x",cols2plot))
      
      p1 = ggplot(test_long, aes(x = SHAP_in, y = SHAP_out)) + 
        geom_point(colour="darkgreen", shape=20, alpha = 0.5) + geom_smooth(method="lm")#scale_y_sqrt()
      p1 = p1 + facet_wrap(~predictor,ncol=nc)#, scales = "free_y") #+ ggtitle(main)
      p1 + theme_bw()
  }
}

plotMDISHAPbias = function(NullSim, PowerSim){
  MDI_power_long <- tidyr::pivot_longer(as.data.frame(t(PowerSim$MDI_sim)),
                    x1:x5,
                    names_to = "predictor",
                    values_to = "MDI")
  MDI_power_long$sim = "power"
  colnames(shap_avs_sim1_power)=paste0("x",1:5)
  SHAP_power_long <- tidyr::pivot_longer(as.data.frame(shap_avs_sim1_power),
                     x1:x5,
                     names_to = "predictor",
                     values_to = "SHAP")
  SHAP_power_long$sim = "power"
  MDI_null_long <- tidyr::pivot_longer(as.data.frame(t(NullSim$MDI_sim)),
                    x1:x5,
                    names_to = "predictor",
                    values_to = "MDI")
  MDI_null_long$sim = "null"
  colnames(shap_avs_sim1_null)=paste0("x",1:5)
  SHAP_null_long <- tidyr::pivot_longer(as.data.frame(shap_avs_sim1_null),
                     x1:x5,
                     names_to = "predictor",
                     values_to = "SHAP")
  SHAP_null_long$sim = "null"

  power_long <- cbind.data.frame(MDI_power_long, SHAP_power_long)[,-c(4,6)]
  null_long <- cbind.data.frame(MDI_null_long, SHAP_null_long)[,-c(4,6)]
  SHAP_MDI_long <- rbind.data.frame(power_long, null_long) 
  
  SHAP_MDI_long <- tidyr::pivot_longer(SHAP_MDI_long,
                   c(MDI,SHAP),
                   names_to = "FI",
                   values_to = "score")
    
  p1 = ggplot(SHAP_MDI_long, aes(x = predictor, y = score)) + 
    geom_boxplot() #+ scale_y_sqrt()
  p1 = p1 + facet_grid(sim~ FI)#, scales = "free_y") #+ ggtitle(main)
  p1 + theme_bw()
}

in_out_long = function(in_out, r){
  res = in_out[[as.character(r)]]
  p = ncol(res$SHAPimp)
  SHAPimp=as.data.frame(res$SHAPimp)
  colnames(SHAPimp)=paste0("x",1:p)
  #browser()
  SHAPimp_long <- tidyr::pivot_longer(SHAPimp,
                                      paste0("x",1:p),
                                        names_to = "predictor",
                                        values_to = "SHAP")
  SHAPimp_long$r = r
  SHAPimp_long$type = "raw"
  
  SHAPimp_shrunk=as.data.frame(res$SHAPimp_shrunk)
  colnames(SHAPimp_shrunk)=paste0("x",1:p)
  SHAPimp_shrunk_long <- tidyr::pivot_longer(SHAPimp_shrunk,
                                      paste0("x",1:p),
                                      names_to = "predictor",
                                      values_to = "SHAP")
  SHAPimp_shrunk_long$r = r
  SHAPimp_shrunk_long$type = "shrunk"
  
  SHAPimp <- rbind.data.frame(SHAPimp_long, SHAPimp_shrunk_long)
  
  if ("SHAPimp_shrunk_oob" %in% names(res) ){
    SHAPimp_shrunk_oob=as.data.frame(res$SHAPimp_shrunk_oob)
    colnames(SHAPimp_shrunk_oob)=paste0("x",1:p)
    SHAPimp_shrunk_oob_long <- tidyr::pivot_longer(SHAPimp_shrunk_oob,
                                              paste0("x",1:p),
                                               names_to = "predictor",
                                               values_to = "SHAP")
    SHAPimp_shrunk_oob_long$r = r
    SHAPimp_shrunk_oob_long$type = "shrunk_oob"
    SHAPimp_shrunk_long$type = "shrunk_test" #in that case also change the previous name
    
    SHAPimp <- rbind.data.frame(SHAPimp_long, SHAPimp_shrunk_long,SHAPimp_shrunk_oob_long)
  }
  
  
  
  return(SHAPimp)
}

plotShrunkSHAP = function(in_out, 
                rVals = c(0,0.15,0.25),
                main = "",
                log_r = FALSE ##<< use logarithmic SNR as labels?
){
  if (missing(rVals)) rVals = as.numeric(names(in_out))
  nFtrs = ncol(in_out[[1]]$SHAPimp)
  SHAPimp = list()
  for (r in rVals){
    tmp = in_out_long(in_out, r)
    if (log_r) tmp$r  = round(log10(tmp$r),1)
    SHAPimp[[as.character(r)]] = tmp
  }
    
  SHAPimp <- do.call("rbind.data.frame", SHAPimp) 
  
  
    p1 = ggplot(SHAPimp, aes(x = predictor, y = SHAP)) + 
      geom_boxplot(outlier.size = 0.1,lwd=0.5, fatten=0.5) #+ scale_y_sqrt()
    p1 = p1 + facet_grid(r ~ type, scales = "free_y") #+ ggtitle(main)
    p1 = p1 + theme_bw() + ggtitle(main)
    if (nFtrs == 5)
      p1 = p1+ scale_x_discrete(labels = expression(x[1], x[2], x[3], x[4], x[5]))
    if (nFtrs == 7)
      p1 = p1+ scale_x_discrete(labels = expression(x[1], x[2], x[3], x[4], x[5], x[6], x[7]))
   return(p1)
      #stat_summary(geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = .75,col=3) 
}