# setwd("C:/Users/colom/bnp_upperbounds/Rscripts/species")
setwd("/home/lucia.paci/Lucia/Ale/bnp_upperbounds/Rscripts/species")

# Librerie ----------------------------------------------------------------

suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
suppressWarnings(suppressPackageStartupMessages(library(progress)))
suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../../R/Rfunctions.R")

# Funzioni ----------------------------------------------------------------

# Common parameters -------------------------------------------------------
alfa = 0.05
Bor  = 3 # nel suo codice, painsky ripete 10k volte, noi 5k?
Bnoi = Bor 
n = 1000
Rmax = 100; RmaxFD = 50
Mgrid = seq(10,1010,by = 50)
Nexp = length(Mgrid)
M_max = 500

seed = 42*idx
set.seed(seed)

# Zipfs -------------------------------------------------------------------
run_zipfs = TRUE

if(run_zipfs){
  s = 1.01
  
  pain_mat = matrix(ub_pain(n = n, Rmax = Rmax, alfa = alfa), Nexp, Bnoi)
  lub_Mrk_mat <- lub_FD_mat <- oracle_mat <- matrix(NA,Nexp,Bnoi)
  
  pb <- progress_bar$new(total = Nexp)
  for(ii in 1:Nexp){
    cat("\n #Exp. = ",ii," \n")
    M = Mgrid[ii]
    # Gen. true distribution
    ptrue = sim_zipfs(M = M, s = s)
    # plot(ptrue)
    
    # Oracle
    Mmax = rep(NA,Bor)
    for(b in 1:Bor){
      # Gen. data
      data = sample(1:M, size = n, replace = TRUE, prob = ptrue)
      n_i = tabulate(data, nbins = M)
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
      
      # oracle[ii] = quantile(Mmax, probs = 1-alfa)
      
      # Markov - Cantelli - FD
      # Param. estimation (PYP)
      start_params <- c(alpha = 0.1, theta = 1)
      fit <- optim(par = start_params, fn = llik_pyp, 
                   method = "L-BFGS-B",
                   lower = c(0, -1), upper = c(1-1e-10, Inf)) 
      alpha_mle = fit$par[1]
      theta_mle = fit$par[2]
      
      # Upper bound (PYP)
      lub_Mrk_mat[ii,b] = compute_log_UBMarkov( Rmax, alpha_mle, theta_mle, 
                                                Kn, n, alfa )
      
      # Param. estimation (FD)
      start_params <- c(gamma = 0.1, Lambda = Kn)
      fit <- optim(par = start_params, fn = llik_FD, 
                   method = "L-BFGS-B",
                   lower = c(1e-5, 1e-5), upper = c(Inf, Inf)) 
      gamma_mle = fit$par[1]
      Lambda_mle = fit$par[2]
      
      # Upper bound (FD)
      lub_FD_mat[ii,b] = compute_log_UBMarkov_FD( RmaxFD, gamma_mle, Lambda_mle, Kn, n, alfa, M_max )
    }
    pb$tick()
  }
  
  save_res = list("lub_Mrk_mat" = lub_Mrk_mat, 
                  "lub_FD_mat"  = lub_FD_mat,
                  "pain_mat"    = pain_mat,
                  "oracle_mat" = oracle_mat)
  
  
  file_name = paste0("../save/SimStudySpecies_zipfs_",idx,".Rdat")
  save(save_res, file = file_name) 
  
}

