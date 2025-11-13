
setwd("C:/Users/colom/BinomialCIs/Rscripts/UnboundedAlphabet")

# Librerie ----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
suppressWarnings(suppressPackageStartupMessages(library(progress)))
suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../../R/Rfunctions.R")

# Sim specific param ----------------------------------------------------------------
RunParallel = FALSE
if(!RunParallel){
  idx = 1 
  param = 1.5
}

names = c("Zipfs","3PBetaPr","Constant","Uniform")
name = names[2] # choose the name here


# Common parameters -------------------------------------------------------
alfa = 0.05
beta = 1e-5
Brep = 100; Bor = 100
seed = 42
set.seed(seed)


# M fix, n varies ---------------------------------------------------------
save_cov = FALSE

## Options -----------------------------------------------------------------

M = 20000
Ngrid_max =  5000
Ngrid_min =  1000
Ngrid_step = 250

Ngrid = seq(Ngrid_min,Ngrid_max,by = Ngrid_step)
Nexp = length(Ngrid)

trimmed_param = get_first3digits(param,4)
exp_name = paste0("SSUbb_Mfix_",name,"_",trimmed_param)
save_exp = FALSE
file_name = paste0("save/",exp_name,".Rdat")
img_name = paste0("img/",name,"/",exp_name,".pdf")
cov_name = paste0("save/Coverage/",name,"/",exp_name,".Rdat")

## Run  --------------------------------------------------------------------

run_M_fix = TRUE
if(run_M_fix){
  
  oracle_mat <- matrix(-Inf,Nexp,Bor) 
  ub_Bench_mat <- ub_An_mat <- ub_Ubb_mat <- matrix(-Inf,Nexp,Brep)
  oracle <- rep(NA,Nexp)
  
  # Coverage
  cov_Bench_mat <- cov_An_mat <- cov_Ubb_mat <- matrix(0,Nexp,Brep) 
  
  # Generate true distribution
  BB = max(Bor,Brep)
  gamma = 10; c = 1000; sigma = 0
  prob_true_mat = r_3ParamBetaPr( BB, M, gamma, sigma, c, seed)
  
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
      pmax = -Inf    
      if(any(n_i == 0)){
        pmax = max( prob_true_mat[b,n_i == 0] ) # get pmax true
      }
      
      ## Benchmark
      Mguess = 10*Kn
      ub_Bench_mat[ii,b] = log(Mguess/alfa)/n
      if(pmax <= ub_Bench_mat[ii,b]){
        cov_Bench_mat[ii,b] = 1
      }
      
      ## LB Analytical
      # lb_An_mat[ii,b] = compute_LB_analytical(n,n_i,alfa)
      
      ## UB Analytical
      bn = log(n)
      n_i_guess = c( n_i[idx_obs], rep(0,Mguess-Kn)  )
      ub_An_mat[ii,b] = compute_UB_analytical(n,n_i_guess,Mguess,bn,alfa,FALSE)
      
      if(pmax <= ub_An_mat[ii,b]){
        cov_An_mat[ii,b] = 1
      }
      
      # Unbounded alphabet:
      Shat = sum( n_i )/n
      
      ## Intersection
      # ub_Ubb_inter_mat[ii,b] = compute_UB_intersection(n, alfa, beta, Shat)
      
      ## r-norm
      Sstar = ( sqrt( -log(beta)/(2*n) ) + sqrt( Shat + (-log(beta)/(2*n)) ) )^2
      rn = log( Sstar / (-log(1-alfa+beta)) ) + log(n) - log(log(n))
      ub_Ubb_mat[ii,b] = compute_UB_rnorm(n, alfa, beta, rn, Shat)
      if(pmax <= ub_Ubb_mat[ii,b]){
        cov_Ubb_mat[ii,b] = 1
      }
    }
    
    pb$tick()
  }
  
  
  save_res = list( "ub_Bench_mat" = ub_Bench_mat,
                   "oracle_mat"   = oracle_mat,
                   "ub_An_mat"    = ub_An_mat,
                   "ub_Ubb_mat"   = ub_Ubb_mat,
                   "oracle" = oracle )
  
  if(save_exp)
    save(save_res, file = file_name)
  
}

## Sample coverage ---------------------------------------------------------
cov_mat = Ngrid
cov_mat = cbind(cov_mat, colMeans( t(cov_Bench_mat) ))
cov_mat = cbind(cov_mat, colMeans( t(cov_Ubb_mat)))
colnames(cov_mat) = c("M","Benchmark","Unbounded")
cov_mat

if(save_cov)
  save(cov_mat, file = cov_name)


## Final summary and plot ---------------------------------------------------
ub_Bench  = apply(ub_Bench_mat, 1, quantile, probs = c(0.025,0.5,0.975))
ub_An     = apply(ub_An_mat, 1, quantile, probs = c(0.025,0.5,0.975))
# ub_Ubb_inter = apply(ub_Ubb_inter_mat, 1, quantile, probs = c(0.025,0.5,0.975))
ub_Ubb_rnorm  = apply(ub_Ubb_mat, 1, quantile, probs = c(0.025,0.5,0.975))

ymax = max(ub_Bench[3,],ub_An[3,],ub_Ubb_rnorm[3,],oracle) * 1.05 #18*1e-3
ymin = min(ub_Bench[1,],ub_An[1,],ub_Ubb_rnorm[1,],oracle) #5*1e-3
if(ymin < 0)
  ymin = 0
ylabs = round(seq(ymin*1e3,ymax*1e3,by = 1),1)

save_img = TRUE

if(save_img)
  pdf(img_name)
par( mfrow = c(1,1), mar = c(4,4,1,0.5), mgp=c(2.5,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "n", ylab = "1000*bound",
     xlim = c(min(Ngrid),max(Ngrid) ) , ylim = c(ymin,ymax), 
     main = paste0(" "),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs*1e-3, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
points( x = Ngrid, y = ub_Bench[2,], 
        type = "l", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkred" ) 
points( x = Ngrid, y = ub_An[2,], 
        type = "l", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkgreen" ) 
points( x = Ngrid, y = ub_Ubb_rnorm[2,],
        type = "l",
        lwd = 3, pch = 16, lty = 1,
        col = "darkorange" )
# points( x = Ngrid, y = lb_An[2,], 
#         type = "l", 
#         lwd = 3, pch = 16, lty = 1,
#         col = "grey45" ) 
# points( x = Ngrid, y = oracle, 
#         type = "l", 
#         lwd = 3, pch = 16, lty = 1,
#         col = "black" )
#        lwd = 3, col = c("darkred","darkgreen","darkblue","grey45","black"))
legend("topright",c("Benchmark","Bounded","Unbounded"), 
       lwd = 3, col = c("darkred","darkgreen","darkorange"))
if(save_img)
  dev.off()
