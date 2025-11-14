
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
beta = 1e-5
Brep = 10; Bor = 10
seed = 42
set.seed(seed)
RunParallel = FALSE
if(!RunParallel)
  idx = 1 


# n fix, M varies ---------------------------------------------------------


## Options -----------------------------------------------------------------

n = 2000
M = 5000
Mfact_grid = c(1,2,10,100,1000,10000)
Nexp = length(Mfact_grid)


exp_name = paste0("SS_overM_nfix_Unif_",idx)
save_exp = FALSE
file_name = paste0("save/",exp_name,".Rdat")
img_name = paste0("img/overM/",exp_name,".pdf")

## Run  --------------------------------------------------------------------


run_n_fix = TRUE
if(run_n_fix){
  
  oracle_mat <- matrix(-Inf,Nexp,Bor) 
  ub_An_mat <- ub_Ubb_rnorm_mat <- matrix(-Inf,Nexp,Brep) 
  
  pb <- progress_bar$new(total = Nexp); ii = 1
  for(ii in 1:Nexp){
    cat("\n #Exp. = ",ii," \n")
    BB = max(Bor,Brep)
    prob_true_mat = matrix(0,nrow = BB, ncol = M) # (BB x M) matrix
    prob_true_mat = t(apply(prob_true_mat, 1, function(x) { sim_TruncatedUnif_features(M = M, Meff = M, pmax = 0.1) } ))
    data = apply(prob_true_mat, 2, function(pj) {rbinom(n = BB, size = n, prob = pj)} ) # (BB x M) matrix
    # data[b,j]: number of obs. of j-th features in b-th repetition
    
    # Estimators
    b = 1
    for(b in 1:Brep){
      n_i = data[b,] # number of obs. of all features in b-th repetition
      idx_obs = which(n_i > 0) # index of observed features
      Kn = length(idx_obs) # number of obs. features 
      data_obs = n_i[idx_obs] # get observed features
      
      if(Kn > 0){
        Mguess = Kn * Mfact_grid[ii]
      }else{
        Mguess = 1
      }      
      ## UB Analytical
      bn = log(n)
      n_i_guess = c( n_i[idx_obs], rep(0,Mguess-Kn)  )
      ub_An_mat[ii,b] = compute_UB_analytical(n,n_i_guess,Mguess,bn,alfa,FALSE)
      
      # Unbounded alphabet:
      Shat = sum( n_i )/n
      
      ## r-norm
      Sstar = ( sqrt( -log(beta)/(2*n) ) + sqrt( Shat + (-log(beta)/(2*n)) ) )^2
      rn = log( Sstar / (-log(1-alfa+beta)) ) + log(n) - log(log(n))
      ub_Ubb_rnorm_mat[ii,b] = compute_UB_rnorm(n, alfa, beta, rn, Shat)
    }
    
    pb$tick()
  }
  
  
}



## Final summary and plot ---------------------------------------------------

ub_An         = apply(ub_An_mat, 1, quantile, probs = c(0.025,0.5,0.975))
ub_Ubb_rnorm  = apply(ub_Ubb_rnorm_mat, 1, quantile, probs = c(0.025,0.5,0.975))

ymax = max(ub_An[3,],ub_Ubb_rnorm[3,]) * 1.05 #18*1e-3
ymin = min(ub_An[1,],ub_Ubb_rnorm[1,]) #5*1e-3
if(ymin < 0)
  ymin = 0
ylabs = round(seq(ymin*1e3,ymax*1e3,by = 1),1)


save_img = FALSE

if(save_img)
  pdf(img_name)
par( mfrow = c(1,1), mar = c(3,3,1,0.5), mgp=c(1.5,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "log(kappa)", ylab = "1000 * bound",
     xlim = range(log(Mfact_grid)) , ylim = c(ymin,ymax), 
     main = paste0(" "),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs*1e-3, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
points( x = log(Mfact_grid), y = ub_An[2,], 
        type = "b", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkgreen" ) 
points( x = log(Mfact_grid), y = ub_Ubb_rnorm[2,], 
        type = "l", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkorange" ) 
legend("bottomright",c("Bounded","Unbounded"), 
       lwd = 3, col = c("darkgreen","darkorange"))
if(save_img)
  dev.off()


