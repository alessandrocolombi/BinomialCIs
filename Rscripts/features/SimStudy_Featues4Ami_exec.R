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
n = 1000

Mgrid_max = 10000
Mgrid = seq(1,Mgrid_max,by = 100)
Nexp = length(Mgrid)

seed = 42
set.seed(seed)
idx = 1 # Ocio!! togli se va in parallelo

# Ami sim study -------------------------------------------------------------------
run_Ami_sim1 = TRUE
if(run_Ami_sim1){
  
  ub_Freq_mat <- ub_Bench_mat <- matrix(-Inf,Nexp,Bnoi) # Frequentist and Benchmark
  oracle_mat  <- matrix(-Inf,Nexp,Bor) # Oracle
  oracle <- rep(NA,Nexp)
  
  pb <- progress_bar$new(total = Nexp)
  for(ii in 1:Nexp){
    cat("\n #Exp. = ",ii," \n")
    M = Mgrid[ii]
    prob_true_mat = matrix(runif(n=Bor*M), nrow = Bor, ncol = M)
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
      cat("\n Shat = ",Shat,"\n")
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
  
  file_name = paste0("../save/SimStudyFeatures_AmiSim1M10k_",idx,".Rdat")
  # save(save_res, file = file_name)
  
}



# Compute final summary ---------------------------------------------------

ub_Freq   = apply(ub_Freq_mat, 1, quantile, probs = c(0.025,0.5,0.975))
ub_Bench  = apply(ub_Bench_mat, 1, quantile, probs = c(0.025,0.5,0.975))
ub_oracle = apply(oracle_mat, 1, quantile, probs = 1-alfa)

# Plot curves -------------------------------------------------------------------
ymax = 18*1e-3
ymin = 5*1e-3
ylabs = round(seq(5,16,by = 1),3)


save_img = FALSE
img_name = paste0("../img/SSFeatures_SimAmi1M500k.pdf")

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
# points( x = Mgrid, y = ub_oracle, 
#         type = "b", lwd = 1, pch = 16, lty = 2,
#         col = "black" ) 
points( x = Mgrid, y = ub_Freq[2,], 
        type = "b", lwd = 1, pch = 16, lty = 2,
        col = "darkorange" ) 
points( x = Mgrid, y = ub_Bench[2,], 
        type = "b", lwd = 1, pch = 16, lty = 2,
        col = "darkred" ) 
legend("bottomright",c("Benchmark","Proposed"), 
       lwd = 3, col = c("darkred","darkorange"))
if(save_img)
  dev.off()



# Ami sim study Zipfs -------------------------------------------------------------------

Mgrid_max = 10000
Mgrid = seq(1,Mgrid_max,by = 100)
Nexp = length(Mgrid)

run_Ami_sim2 = TRUE
if(run_Ami_sim2){
  
  ub_Freq_mat <- ub_Bench_mat <- matrix(-Inf,Nexp,Bnoi) # Frequentist and Benchmark
  oracle_mat  <- matrix(-Inf,Nexp,Bor) # Oracle
  oracle <- rep(NA,Nexp)
  
  pb <- progress_bar$new(total = Nexp)
  for(ii in 1:Nexp){
    s = 1.01
    cat("\n #Exp. = ",ii," \n")
    M = Mgrid[ii]
    # Gen. true distribution
    ptrue = sim_zipfs_features(M = M, s = s)
    # Gen. all datasets
    data = lapply(ptrue, function(pj) rbinom(n = Bor, size = n, prob = pj) )
    data = do.call(cbind, data) # Bor x M matrix. data[b,] are the probabilities for the b-th replicate
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
        Mmax[b] = max(ptrue[idx_unobs])
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
  
  file_name = paste0("../save/SimStudyFeatures_AmiSim2M10k_",idx,".Rdat")
  # save(save_res, file = file_name)
  
}



# Compute final summary ---------------------------------------------------
ub_Freq   = apply(ub_Freq_mat, 1, quantile, probs = c(0.025,0.5,0.975))
ub_Bench  = apply(ub_Bench_mat, 1, quantile, probs = c(0.025,0.5,0.975))
ub_oracle = apply(oracle_mat, 1, quantile, probs = 1-alfa)

# Plot curves -------------------------------------------------------------------
ymax = 13*1e-3
ymin = 5*1e-3
ylabs = round(seq(5,13,by = 1),3)


save_img = FALSE
img_name = paste0("../img/SSFeatures_SimAmi1_zipfs_M10k.pdf")

if(save_img)
  pdf(img_name)
par( mfrow = c(1,1), mar = c(3,3,1,0.5), mgp=c(1.5,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "M", ylab = "1000 * bound",
     xlim = c(0,max(Mgrid) ) , ylim = c(ymin,ymax), 
     main = paste0(" Zipfs "),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs*1e-3, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
# points( x = Mgrid, y = ub_oracle, 
#         type = "b", lwd = 1, pch = 16, lty = 2,
#         col = "black" ) 
points( x = Mgrid, y = ub_Freq[2,], 
        type = "b", lwd = 1, pch = 16, lty = 2,
        col = "darkorange" ) 
points( x = Mgrid, y = ub_Bench[2,], 
        type = "b", lwd = 1, pch = 16, lty = 2,
        col = "darkred" ) 
legend("bottomright",c("Benchmark","Proposed"), 
       lwd = 3, col = c("darkred","darkorange"))
if(save_img)
  dev.off()

setwd("C:/Users/colom/bnp_upperbounds/Rscripts/features")
# setwd("/home/lucia.paci/Lucia/Ale/bnp_upperbounds/Rscripts/features")



# Few data points -------------------------------------------------------
n = 10

# Ami sim study -------------------------------------------------------------------
run_Ami_sim1 = TRUE
if(run_Ami_sim1){
  
  ub_Freq_mat <- ub_Bench_mat <- matrix(-Inf,Nexp,Bnoi) # Frequentist and Benchmark
  oracle_mat  <- matrix(-Inf,Nexp,Bor) # Oracle
  oracle <- rep(NA,Nexp)
  
  pb <- progress_bar$new(total = Nexp)
  for(ii in 1:Nexp){
    cat("\n #Exp. = ",ii," \n")
    M = Mgrid[ii]
    prob_true_mat = matrix(runif(n=Bor*M), nrow = Bor, ncol = M)
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
  
  file_name = paste0("../save/SimStudyFeatures_AmiSim1M10k_",idx,".Rdat")
  # save(save_res, file = file_name)
  
}

# Compute final summary ---------------------------------------------------
ub_Freq   = apply(ub_Freq_mat, 1, quantile, probs = c(0.025,0.5,0.975))
ub_Bench  = apply(ub_Bench_mat, 1, quantile, probs = c(0.025,0.5,0.975))
ub_oracle = apply(oracle_mat, 1, quantile, probs = 1-alfa)

# Plot curves -------------------------------------------------------------------
ymax = max(c(ub_Bench[2,],ub_Freq[2,]))
ymin = 0
ylabs = round(seq(ymin,ymax,length.out = 5),3)


save_img = TRUE
img_name = paste0("../img/SSFeatures_SimAmi1n10M10k.pdf")

if(save_img)
  pdf(img_name)
par( mfrow = c(1,1), mar = c(3,3,1,0.5), mgp=c(1.5,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "M", ylab = "bound",
     xlim = c(0,max(Mgrid) ) , ylim = c(ymin,ymax), 
     main = paste0(" "),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
# points( x = Mgrid, y = ub_oracle, 
#         type = "b", lwd = 1, pch = 16, lty = 2,
#         col = "black" ) 
points( x = Mgrid, y = ub_Freq[2,], 
        type = "b", lwd = 1, pch = 16, lty = 2,
        col = "darkorange" ) 
points( x = Mgrid, y = ub_Bench[2,], 
        type = "b", lwd = 1, pch = 16, lty = 2,
        col = "darkred" ) 
legend("bottomright",c("Benchmark","Proposed"), 
       lwd = 3, col = c("darkred","darkorange"))
if(save_img)
  dev.off()


