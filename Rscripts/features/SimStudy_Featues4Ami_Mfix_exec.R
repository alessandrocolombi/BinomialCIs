setwd("C:/Users/colom/bnp_upperbounds/Rscripts/features")
# setwd("/home/lucia.paci/Lucia/Ale/bnp_upperbounds/Rscripts/features")

# Librerie ----------------------------------------------------------------

suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
suppressWarnings(suppressPackageStartupMessages(library(progress)))
suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../../R/Rfunctions.R")

# Funzioni ----------------------------------------------------------------
sim_unif_features = function(M = M){
  runif(n = M)
}
# Common parameters -------------------------------------------------------
alfa = 0.05
Bor  = 100 # number of repetitions
Bnoi <- Bfr <- Bor
M = 10000
Ngrid_max = 10000
Ngrid = seq(1,Ngrid_max,by = 100)
Nexp = length(Ngrid)

seed = 42
set.seed(seed)


# Ami sim study -------------------------------------------------------------------
run_Ami_simMfix = TRUE
if(run_Ami_simMfix){
  
  ub_Freq_mat <- ub_Bench_mat <- matrix(-Inf,Nexp,Bnoi) # Frequentist and Benchmark
  oracle_mat  <- matrix(-Inf,Nexp,Bor) # Oracle
  oracle <- rep(NA,Nexp)
  
  prob_true_mat = matrix(runif(n=Bor*M), nrow = Bor, ncol = M)
  
  pb <- progress_bar$new(total = Nexp)
  for(ii in 1:Nexp){
    cat("\n #Exp. = ",ii," \n")
    n = Ngrid[ii]
    data = apply(prob_true_mat, 2, function(pj) {rbinom(n = Bor, size = n, prob = pj)} )
    # Oracle
    Mmax = rep(NA,Bor)
    for(b in 1:Bor){
      n_i = data[b,]
      idx_obs = which(n_i > 0)
      Kn = length(idx_obs)
      data_obs = n_i[idx_obs]
      if(Kn == M){
        Mmax[b] = 0
      }else{
        idx_unobs = which(n_i == 0)
        Mmax[b] = max(prob_true_mat[b,idx_unobs])
      }
      oracle_mat[ii,b] = Mmax[b]
    }
    oracle[ii] = quantile(Mmax, probs = 1-alfa)
    
    # Frequentist estimators
    for(b in 1:Bfr){
      # Read data
      n_i = data[b,]
      idx_obs = which(n_i > 0)
      Kn = length(idx_obs)
      data_obs = n_i[idx_obs]
      Shat = sum( data_obs/n )
      
      # Minimization wrt beta
      start_params <- c(beta = alfa/2)
      fit <- optim(par = start_params, fn = lf_beta, 
                   method = "L-BFGS-B",
                   n = n, alfa = alfa, Shat = Shat,
                   lower = 1e-16, upper = alfa - 1e-10 ) 
      beta_opt = fit$par[1]
      # Upper bound (freq)
      ub_Freq_mat[ii,b] = compute_UBFreq_BeBe(n, alfa, beta_opt, Shat)
      
      # Oracle estimator:
      ub_Bench_mat[ii,b] = log(M/alfa)/n
    }
    pb$tick()
  }
  
  
  save_res = list("ub_Freq_mat"  = ub_Freq_mat,
                  "ub_Bench_mat" = ub_Bench_mat,
                  "oracle_mat"  = oracle_mat,
                  "oracle" = oracle)
  
  # file_name = paste0("../save/SimStudyFeatures_AmiSim2N10k_",idx,".Rdat")
  # save(save_res, file = file_name)
  
}



# Compute final summary ---------------------------------------------------
ub_Freq   = apply(ub_Freq_mat, 1, quantile, probs = c(0.025,0.5,0.975))
ub_Bench  = apply(ub_Bench_mat, 1, quantile, probs = c(0.025,0.5,0.975))
ub_oracle = apply(oracle_mat, 1, quantile, probs = 1-alfa)

# Plot curves -------------------------------------------------------------------
ymax = 0.025
ymin = 0
ylabs = round(seq(0,0.025,length.out = 6),3)


save_img = FALSE
img_name = paste0("../img/SSFeatures_SimAmi2M10k.pdf")

if(save_img)
  pdf(img_name)
par( mfrow = c(1,1), mar = c(3,3,1,0.5), mgp=c(1.5,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "n", ylab = "bound",
     xlim = c(0,max(Ngrid) ) , ylim = c(ymin,ymax), 
     main = paste0(" "),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
points( x = Ngrid, y = ub_oracle,
        type = "b", lwd = 1, pch = 16, lty = 2,
        col = "black" )
points( x = Ngrid, y = ub_Freq[2,], 
        type = "b", lwd = 1, pch = 16, lty = 2,
        col = "darkorange" ) 
points( x = Ngrid, y = ub_Bench[2,], 
        type = "b", lwd = 1, pch = 16, lty = 2,
        col = "darkred" ) 
legend("topright",c("Oracle","Benchmark","Proposed"), 
       lwd = 3, col = c("black","darkred","darkorange"))
if(save_img)
  dev.off()
