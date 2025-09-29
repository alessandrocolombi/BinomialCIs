
setwd("C:/Users/colom/BinomialCIs/Rscripts/UnboundedAlphabet")

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
Brep = 100; Bor = 10000
seed = 42
set.seed(seed)
RunParallel = FALSE
if(!RunParallel)
  idx = 1 

# n varies ---------------------------------------------------------

eps_truncation = 1e-10

## Options -----------------------------------------------------------------

Ngrid_max = 200
Ngrid_step = 20

Ngrid = seq(10,Ngrid_max,by = Ngrid_step)
Nexp = length(Ngrid)

s = 1.025

exp_name = paste0("SimStudy_Mfix_Zipfs_",idx)
save_exp = FALSE
file_name = paste0("save/",exp_name,".Rdat")
img_name = paste0("save/",exp_name,".pdf")

## Run  --------------------------------------------------------------------

beta = 1e-5

run_M_fix = TRUE
if(run_M_fix){
  
  oracle_mat <- n_oracle_mat <- matrix(-Inf,Nexp,Bor) 
  ub_Inter_mat <- ub_rnorm_mat <- matrix(-Inf,Nexp,Brep) 
  nub_Inter_mat <- nub_rnorm_mat <- matrix(-Inf,Nexp,Brep) 
  oracle <- n_oracle <- rep(NA,Nexp)
  
  # Generate true distribution
  BB = max(Bor,Brep)
  
  prob_true_list = vector("list", length = Nexp) # prob_true_list[n] truncate Zipfs according to n
  prob_true_list = lapply(Ngrid, function(n) { 
    x = sim_zipfs_features(s = s, n = n, eps = eps_truncation)
    Mn = length(x)
    prob_true_mat = matrix(0,nrow = BB, ncol = Mn) # (BB x Mn) matrix
    t(apply(prob_true_mat, 1, function(y) { x } ))
  } )
  
  pb <- progress_bar$new(total = Nexp)
  for(ii in 1:Nexp){
    cat("\n #Exp. = ",ii," \n")
    n = Ngrid[ii]
    Mn = ncol(prob_true_list[[ii]])
    
    data = apply(prob_true_list[[ii]], 2, function(pj) {rbinom(n = BB, size = n, prob = pj)} ) # (BB x M) matrix
    # data[b,j]: number of obs. of j-th features in b-th repetition
    
    # Oracle
    Mmax = rep(NA,Bor) # needed to compute the oracle
    for(b in 1:Bor){
      n_i = data[b,] # number of obs. of all features in b-th repetition
      idx_obs = which(n_i > 0) # index of observed features
      Kn = length(idx_obs) # number of obs. features 
      data_obs = n_i[idx_obs] # get observed features
      if(Kn == Mn){
        # if here, the whole alphabet has been observed
        Mmax[b] = 0
      }else{
        # if here, some symbols have not been observed
        idx_unobs = which(n_i == 0)
        Mmax[b] = max(prob_true_list[[ii]][b,idx_unobs])
      }
      
      ## Get pmax = max{pj : Nj = 0} in b-th repetition
      oracle_mat[ii,b] = Mmax[b]  
      n_oracle_mat[ii,b] = n * oracle_mat[ii,b]
    }
    
    ## Oracle 
    oracle[ii] = quantile(Mmax, probs = 1-alfa)
    n_oracle[ii] = n * oracle[ii]
    
    # Estimators
    for(b in 1:Brep){
      n_i = data[b,] # number of obs. of all features in b-th repetition
      idx_obs = which(n_i > 0) # index of observed features
      Kn = length(idx_obs) # number of obs. features 
      data_obs = n_i[idx_obs] # get observed features
      # Compute Shat
      Shat = sum( data_obs )/n
      
      ## Intersection
      ub_Inter_mat[ii,b] = compute_UB_intersection(n, alfa, beta, Shat)
      nub_Inter_mat[ii,b] = n * ub_Inter_mat[ii,b]
      ## r-norm
      Sstar = ( sqrt( -log(beta)/(2*n) ) + sqrt( Shat + (-log(beta)/(2*n)) ) )^2
      rn = log( Sstar / (-log(1-alfa+beta)) ) + log(n) - log(log(n))
      ub_rnorm_mat[ii,b] = compute_UB_rnorm(n, alfa, beta, rn, Shat)
      nub_rnorm_mat[ii,b] = n * ub_rnorm_mat[ii,b]
    }
    
    pb$tick()
  }
  
  
  save_res = list( "ub_Inter_mat" = ub_Inter_mat,
                   "nub_Inter_mat" = nub_Inter_mat,
                   "ub_rnorm_mat" = ub_rnorm_mat,
                   "nub_rnorm_mat" = nub_rnorm_mat,
                   "oracle_mat"  = oracle_mat,
                   "n_oracle_mat"  = n_oracle_mat,
                   "oracle" = oracle, 
                   "n_oracle" = n_oracle
                   )
  
  if(save_exp)
    save(save_res, file = file_name)
  
}



## Final summary and plot ---------------------------------------------------

ub_Inter  = apply(ub_Inter_mat, 1, quantile, probs = c(0.025,0.5,0.975))
ub_rnorm  = apply(ub_rnorm_mat, 1, quantile, probs = c(0.025,0.5,0.975))

nub_Inter  = apply(nub_Inter_mat, 1, quantile, probs = c(0.025,0.5,0.975))
nub_rnorm  = apply(nub_rnorm_mat, 1, quantile, probs = c(0.025,0.5,0.975))




ymax = max( c( nub_Inter[3,],nub_rnorm[3,],n_oracle)) * 1.05
ymin = min( c( nub_Inter[1,],nub_rnorm[1,],n_oracle)) 
ylabs = round(seq(ymin,ymax,length.out = 5),1)


save_img = FALSE

if(save_img)
  pdf(img_name)
par( mfrow = c(1,1), mar = c(3.5,3.5,1,0.5), mgp=c(2.5,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "n", ylab = "n x Length",
     xlim = c(0,max(Ngrid) ) , ylim = c(ymin,ymax), 
     main = paste0(" "),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
points( x = Ngrid, y = nub_Inter[2,], 
        type = "l", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkred" )
points( x = Ngrid, y = nub_rnorm[2,], 
        type = "l", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkorange" ) 
points( x = Ngrid, y = n_oracle, 
        type = "l", 
        lwd = 3, pch = 16, lty = 1,
        col = "black" )
legend("topleft",c("Intersection","r-norm","Oracle"), 
       lwd = 3, col = c("darkred","darkorange","black"))
if(save_img)
  dev.off()



# Nota per me: se per n grande vedi l'oracle droppare a zero, 
# allora vuol dire che il livello di truncation della zipfs è troppo grande 
# e vedo tutti i simboli