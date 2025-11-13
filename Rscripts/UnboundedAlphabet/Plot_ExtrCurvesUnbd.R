
setwd("C:/Users/colom/BinomialCIs/Rscripts/UnboundedAlphabet")

# Librerie ----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
suppressWarnings(suppressPackageStartupMessages(library(progress)))
suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../../R/Rfunctions.R")



# 1 param Beta Proc Plot extrapolation curves -----------------------------------------------
seed = 1234
set.seed(seed)
seeds = sample(1:999999, size = 5000, replace = FALSE)


gammas = c(5,10,100,1000)
param_grid = expand.grid(gammas)

M = 5000
n = 2500
ngrid = ceiling(seq(1,n, length.out = 100))
BB = 50


save_plot_all = FALSE

i = 1
for(i in 1:nrow(param_grid)){
  cat("\n Start: ",i,"/",nrow(param_grid),"\n")
  if(save_plot_all)
    pdf(paste0("img/","BetaPr_ExtCurve",".pdf"))
  par(mfrow = c(1,1),bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,1,0))
  
  gamma = param_grid[i,1]
  trimmed_gamma = get_first3digits(gamma,5)
  prob_true_mat = r_BetaPr(BB, M, gamma, seed[i]) # (BB x M) matrix
    
  NumObs = matrix(0,nrow = length(ngrid), ncol = BB)
  for(nn in 1:length(ngrid)){
    data = apply(prob_true_mat, 2, function(pj) {rbinom(n = BB, size = ngrid[nn], prob = pj)} ) # (BB x M) matrix
    NumObs[nn, ] = rowSums(data > 0)
  }
    
  ExtrCrv = apply(NumObs, 1, quantile, probs = c(0.025,0.5,0.975))
  plot(x = 0, y = 0, type = "n",
       main = paste0("BetaPr"," - ", trimmed_gamma,"0,1"), 
       xlab = "#obs.", ylab = "#distincts",
       ylim = c(0,max(ExtrCrv)+1),
       xlim = c(0,n+1),
       pch = 1) # init plot
  polygon( c(ngrid, rev(ngrid)),
           c(ExtrCrv[1,], rev(ExtrCrv[3,])),
           col = "grey75",
           border = NA) # plot in-sample bands
    points(x = ngrid, y = ExtrCrv[2,], type = "l", lwd = 3) # plot mean obs
  if(save_plot_all)
    dev.off()
}



















# 3 param Beta Proc Plot extrapolation curves -----------------------------------------------
seed = 1234
set.seed(seed)
seeds = sample(1:999999, size = 5000, replace = FALSE)


gammas = c(10)
sigmas = c(0)
cs     = c(1000)
param_grid = expand.grid(gammas,sigmas,cs)

M = 5000
n = 2500
ngrid = ceiling(seq(1,n, length.out = 100))
BB = 50


save_plot_all = TRUE

i = 1
for(i in 1:nrow(param_grid)){
  cat("\n Start: ",i,"/",nrow(param_grid),"\n")
  if(save_plot_all)
    pdf(paste0("img/","3paramBetaPr_ExtCurve",".pdf"))
  par(mfrow = c(1,1),bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,1,0))
  
  gamma = param_grid[i,1]
  sigma = param_grid[i,2]
  c = param_grid[i,3]
  
  trimmed_gamma = get_first3digits(gamma,5)
  trimmed_sigma = get_first3digits(sigma,5)
  trimmed_c = get_first3digits(c,5)
  
  prob_true_mat =  r_3ParamBetaPr( BB, M, gamma, sigma, c, seed[i]) # (BB x M) matrix
  
  NumObs = matrix(0,nrow = length(ngrid), ncol = BB)
  for(nn in 1:length(ngrid)){
    data = apply(prob_true_mat, 2, function(pj) {rbinom(n = BB, size = ngrid[nn], prob = pj)} ) # (BB x M) matrix
    NumObs[nn, ] = rowSums(data > 0)
  }
  
  ExtrCrv = apply(NumObs, 1, quantile, probs = c(0.025,0.5,0.975))
  plot(x = 0, y = 0, type = "n",
       main = paste0("3paramBetaPr"," - ", trimmed_gamma,",",trimmed_sigma,",",trimmed_c), 
       xlab = "#obs.", ylab = "#distincts",
       ylim = c(0,max(ExtrCrv)+1),
       xlim = c(0,n+1),
       pch = 1) # init plot
  polygon( c(ngrid, rev(ngrid)),
           c(ExtrCrv[1,], rev(ExtrCrv[3,])),
           col = "grey75",
           border = NA) # plot in-sample bands
  points(x = ngrid, y = ExtrCrv[2,], type = "l", lwd = 3) # plot mean obs
  if(save_plot_all)
    dev.off()
}





# Brutta ------------------------------------------------------------------
Nrep = 10
Natoms = 1000
gamma = 10000
sigma = 0.999999
c = 10000
seed = 1234

pmat =  r_3ParamBetaPr( Nrep, Natoms, gamma, sigma, c, seed)


par(mfrow = c(1,1),bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,1,0))
plot(x = 0, y = 0, type = "n", 
     xlab = "Prob.", ylab = "Atoms",
     ylim = c(0,1), xlim = c(0,Natoms)) # init plot
for(i in 1:nrow(pmat)){
  points(x = 1:Natoms, y = sort(pmat[i,], decreasing = TRUE), 
         type = "b", lty = 2, lwd = 2, pch = 16)
}



Nrep = 10
Natoms = 100
gamma = 1000
sigma = 0
c = 100
seed = 1234

pmat =  r_BetaPr( Nrep, Natoms, gamma, seed)


par(mfrow = c(1,1),bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,1,0))
plot(x = 0, y = 0, type = "n", 
     xlab = "Prob.", ylab = "Atoms",
     ylim = c(0,1), xlim = c(0,Natoms)) # init plot
for(i in 1:nrow(pmat)){
  points(x = 1:Natoms, y = sort(pmat[i,], decreasing = TRUE), 
         type = "b", lty = 2, lwd = 2, pch = 16)
}
