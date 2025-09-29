
setwd("C:/Users/colom/BinomialCIs/Rscripts/BoundedAlphabet")

# Librerie ----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
suppressWarnings(suppressPackageStartupMessages(library(progress)))
suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../../R/Rfunctions.R")

# Custom functions ----------------------------------------------------------------

# Common parameters -------------------------------------------------------
alfa = 0.05
Brep = 2; Bor = 1000
seed = 42
set.seed(seed)
RunParallel = FALSE
if(!RunParallel)
  idx = 1 


# n fix, M varies ---------------------------------------------------------


## Options -----------------------------------------------------------------

n = 1000
Mgrid_max = 2000
Mgrid_step = 100

Mgrid = seq(10,Mgrid_max,by = Mgrid_step)
Nexp = length(Mgrid)

s = 1.02

exp_name = paste0("SimStudy_nfix_Zipfs_",idx)
save_exp = FALSE
file_name = paste0("save/",exp_name,".Rdat")
img_name = paste0("save/",exp_name,".pdf")

## Run  --------------------------------------------------------------------


run_n_fix = TRUE
if(run_n_fix){
  
  oracle_mat <- matrix(-Inf,Nexp,Bor) 
  ub_Bench_mat <- matrix(-Inf,Nexp,Brep) 
  oracle <- rep(NA,Nexp)
  
  pb <- progress_bar$new(total = Nexp)
  for(ii in 1:Nexp){
    cat("\n #Exp. = ",ii," \n")
    M = Mgrid[ii]
    BB = max(Bor,Brep)
    prob_true_mat = matrix(0,nrow = BB, ncol = M) # (BB x M) matrix
    prob_true_mat = t(apply(prob_true_mat, 1, function(x) { sim_TruncatedZipfs_features(M = M, s = s) } ))
    data = apply(prob_true_mat, 2, function(pj) {rbinom(n = BB, size = n, prob = pj)} ) # (BB x M) matrix
    # data[b,j]: number of obs. of j-th features in b-th repetition
    
    # Oracle
    Mmax = rep(NA,Bor) # needed to compute the oracle
    for(b in 1:Bor){
      n_i = data[b,] # number of obs. of all features in b-th repetition
      idx_obs = which(n_i > 0) # index of observed features
      Kn = length(idx_obs) # number of obs. features 
      data_obs = n_i[idx_obs] # get observed features
      if(Kn == M){
        # if here, the whole alphabet has been observed
        Mmax[b] = 0
      }else{
        # if here, some symbols have not been observed
        idx_unobs = which(n_i == 0)
        Mmax[b] = max(prob_true_mat[b,idx_unobs])
      }
      
      ## Get pmax = max{pj : Nj = 0} in b-th repetition
      oracle_mat[ii,b] = Mmax[b]  
      
    }
    
    ## Oracle 
    oracle[ii] = quantile(Mmax, probs = 1-alfa)
    
    # Estimators
    for(b in 1:Brep){
      n_i = data[b,] # number of obs. of all features in b-th repetition
      idx_obs = which(n_i > 0) # index of observed features
      Kn = length(idx_obs) # number of obs. features 
      data_obs = n_i[idx_obs] # get observed features
      
      ## Benchmark
      ub_Bench_mat[ii,b] = log(M/alfa)/n
    }
    
    pb$tick()
  }
  
  
  save_res = list( "ub_Bench_mat" = ub_Bench_mat,
                   "oracle_mat"  = oracle_mat,
                   "oracle" = oracle )
  
  if(save_exp)
    save(save_res, file = file_name)
  
}



## Final summary and plot ---------------------------------------------------

ub_Bench  = apply(ub_Bench_mat, 1, quantile, probs = c(0.025,0.5,0.975))

ymax = max(ub_Bench[3,],oracle) * 1.05 #18*1e-3
ymin = min(ub_Bench[1,],oracle) #5*1e-3
ylabs = round(seq(ymin*1e3,ymax*1e3,by = 1),3)


save_img = FALSE

if(save_img)
  pdf(img_name)
par( mfrow = c(1,1), mar = c(3,3,1,0.5), mgp=c(1.5,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "M", ylab = "1000 * bound",
     xlim = c(0,max(Mgrid) ) , ylim = c(ymin,ymax), 
     main = paste0(" "),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs*1e-3, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
points( x = Mgrid, y = ub_Bench[2,], 
        type = "l", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkred" ) 
points( x = Mgrid, y = oracle, 
        type = "l", 
        lwd = 3, pch = 16, lty = 1,
        col = "black" )
legend("bottomright",c("Benchmark","Proposed"), 
       lwd = 3, col = c("darkred","black"))
if(save_img)
  dev.off()


# M fix, n varies ---------------------------------------------------------


## Options -----------------------------------------------------------------

M = 1000
Ngrid_max = 2000
Ngrid_step = 100

Ngrid = seq(10,Ngrid_max,by = Ngrid_step)
Nexp = length(Ngrid)

s = 1.02

exp_name = paste0("SimStudy_Mfix_Zipfs_",idx)
save_exp = FALSE
file_name = paste0("save/",exp_name,".Rdat")
img_name = paste0("save/",exp_name,".pdf")

## Run  --------------------------------------------------------------------


run_M_fix = TRUE
if(run_M_fix){
  
  oracle_mat <- matrix(-Inf,Nexp,Bor) 
  ub_Bench_mat <- matrix(-Inf,Nexp,Brep) 
  oracle <- rep(NA,Nexp)
  
  # Generate true distribution
  BB = max(Bor,Brep)
  prob_true_mat = matrix(0,nrow = BB, ncol = M) # (BB x M) matrix
  prob_true_mat = t(apply(prob_true_mat, 1, function(x) { sim_TruncatedZipfs_features(M = M, s = s) } ))
  
  pb <- progress_bar$new(total = Nexp)
  for(ii in 1:Nexp){
    cat("\n #Exp. = ",ii," \n")
    n = Ngrid[ii]
    data = apply(prob_true_mat, 2, function(pj) {rbinom(n = BB, size = n, prob = pj)} ) # (BB x M) matrix
    # data[b,j]: number of obs. of j-th features in b-th repetition
    
    # Oracle
    Mmax = rep(NA,Bor) # needed to compute the oracle
    for(b in 1:Bor){
      n_i = data[b,] # number of obs. of all features in b-th repetition
      idx_obs = which(n_i > 0) # index of observed features
      Kn = length(idx_obs) # number of obs. features 
      data_obs = n_i[idx_obs] # get observed features
      if(Kn == M){
        # if here, the whole alphabet has been observed
        Mmax[b] = 0
      }else{
        # if here, some symbols have not been observed
        idx_unobs = which(n_i == 0)
        Mmax[b] = max(prob_true_mat[b,idx_unobs])
      }
      
      ## Get pmax = max{pj : Nj = 0} in b-th repetition
      oracle_mat[ii,b] = Mmax[b]  
      
    }
    
    ## Oracle 
    oracle[ii] = quantile(Mmax, probs = 1-alfa)
    
    # Estimators
    for(b in 1:Brep){
      n_i = data[b,] # number of obs. of all features in b-th repetition
      idx_obs = which(n_i > 0) # index of observed features
      Kn = length(idx_obs) # number of obs. features 
      data_obs = n_i[idx_obs] # get observed features
      
      ## Benchmark
      ub_Bench_mat[ii,b] = log(M/alfa)/n
    }
    
    pb$tick()
  }
  
  
  save_res = list( "ub_Bench_mat" = ub_Bench_mat,
                   "oracle_mat"  = oracle_mat,
                   "oracle" = oracle )
  
  if(save_exp)
    save(save_res, file = file_name)
  
}



## Final summary and plot ---------------------------------------------------

ub_Bench  = apply(ub_Bench_mat, 1, quantile, probs = c(0.025,0.5,0.975))

ymax = max(ub_Bench[3,],oracle) * 1.05 #18*1e-3
ymin = min(ub_Bench[1,],oracle) #5*1e-3
# ylabs = round(seq(ymin*1e3,ymax*1e3,by = 1),0)
ylabs = round(seq(ymin,ymax,length.out = 5),3)


save_img = FALSE

if(save_img)
  pdf(img_name)
par( mfrow = c(1,1), mar = c(4,4,1,0.5), mgp=c(3,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "M", ylab = "bound",
     xlim = c(0,max(Ngrid) ) , ylim = c(ymin,ymax), 
     main = paste0(" "),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
points( x = Ngrid, y = ub_Bench[2,], 
        type = "l", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkred" ) 
points( x = Ngrid, y = oracle, 
        type = "l", 
        lwd = 3, pch = 16, lty = 1,
        col = "black" )
legend("topright",c("Benchmark","Proposed"), 
       lwd = 3, col = c("darkred","black"))
if(save_img)
  dev.off()