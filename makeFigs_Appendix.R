source("src/helperFuns.R")

#recreating the simulation data takes a long time
rerun=FALSE

## Figure 7

n=120;M=100
rdaFile = paste0("data/in_out_xgb_n",n,".rda")
pars = list(n=n,M=M, alg = "xgb", genData = "strobl", desc="bla")

if (rerun){
  in_out=list()
  r=0.15;cat("r=", r, "\n")
  in_out[[as.character(r)]] = SHAPcorrs(r,n=n, M=M, alg = "xgb")
  
  r=0.25;cat("r=", r, "\n")
  in_out[[as.character(r)]] = SHAPcorrs(r,n=n, M=M, alg = "xgb")
  
  r=0;cat("r=", r, "\n")
  in_out[[as.character(r)]] = SHAPcorrs(r,n=n, M=M, alg = "xgb")
  
  save(in_out, pars, file = rdaFile)
} else {
  load(rdaFile)
}

plotShrunkSHAP(in_out)
fname=paste0("figures/ShrunkSHAP_xgb.pdf")
ggsave(fname, height=5, width=5)

## Figure 8, all features are continuous
k=0
#in_out_xgb_MARS_k0_n120_p7.pdf
in_out_cont1 = RunSim(n=120,M=100,p=7,k=k,alg = "xgb",genData = "MARS", main ="", 
                          rerun = FALSE, addSims = FALSE,height=4.5, width=5)

## Figure 9, first 3 are discrete
#rdaFile="data/in_out_xgb_MARS_k22200_n120_p5.rda"
p=5
k=rep(0,p);k[1:3]=2
in_out_discrete1 = RunSim(n=120,M=100,p=5,k=k,alg = "xgb",genData = "MARS", main ="", 
                 rerun = FALSE, addSims = FALSE,height=4.5, width=5, verbose = 1)
                          
## Figure 10
#rdaFile="data/May13/in_out_xgb_MARS_k2220000_n120_p7.rda"
p=7
k=rep(0,p);k[1:3]=2
in_out_discrete2 = RunSim(n=120,M=100,p=7,k=k,alg = "xgb",genData = "MARS", 
               main ="",rerun = FALSE, addSims = FALSE,height=4.5, width=5)
