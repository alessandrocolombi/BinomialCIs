setwd("C:/Users/colom/bnp_upperbounds/Rscripts/features")

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
Bor  = 50 # nel suo codice, painsky ripete 10k volte, noi 5k?
Bnoi <- Bfr <- 50
n = 1000
Rmax = 100; 
Mgrid = seq(10,1010,by = 50)
Nexp = length(Mgrid)
M_max = 500

seed = 42*idx
set.seed(seed)


# Zipfs -------------------------------------------------------------------
s = 1.01
var_gamma = 10
var_nb    = 10

lub_PP_mat <- lub_PP2_mat <- lub_PP3_mat <- matrix(-Inf,Nexp,Bnoi) # Poisson process
lub_MixPois_mat <- lub_MixBin_mat <- matrix(-Inf,Nexp,Bnoi) # Mixed Poisson/Binomial
ub_Freq_mat <- matrix(-Inf,Nexp,Bnoi)
oracle_mat  <- matrix(-Inf,Nexp,Bor)
oracle <- rep(NA,Nexp)

pb <- progress_bar$new(total = Nexp)
for(ii in 1:Nexp){
  cat("\n #Exp. = ",ii," \n")
  M = Mgrid[ii]
  # Gen. true distribution
  ptrue = sim_zipfs_features(M = M, s = s)
  # plot(ptrue)
  # Gen. all datasets
  data = lapply(ptrue, function(pj) rbinom(n = Bor, size = n, prob = pj) )
  data = do.call(cbind, data)
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
  # Proposed estimators
  
  for(b in 1:Bnoi){
    # Gen. data
    n_i = data[b,]
    idx_obs = which(n_i > 0)
    Kn = length(idx_obs)
    data_obs = n_i[idx_obs]
    
    # # Param. estimation (PP, gamma = 1)
    # gamma = 1
    # start_params <- c(alpha = 0.1, u = 1)
    # fit <- optim(par = start_params, fn = llik_PP, 
    #              method = "L-BFGS-B",
    #              lower = c(0, 1e-16), upper = c(1-1e-10, Inf)) 
    # alpha_mle = fit$par[1]
    # c_mle = fit$par[2] - alpha_mle
    # 
    # # Upper bound (PP, gamma = 1)
    # lub_PP_mat[ii,b] = compute_log_UBMarkov_BeBePois( Rmax, alpha_mle, c_mle, gamma, n, alfa )
    # 
    # # Param. estimation (PP, alpha = 0)
    # alpha = 0
    # start_params <- c(gamma = 1, c = 1)
    # fit <- optim(par = start_params, fn = llik_PP2, 
    #              method = "L-BFGS-B",
    #              lower = c(1e-16, 1e-16), upper = c(Inf, Inf)) 
    # gamma_mle = fit$par[1]
    # c_mle     = fit$par[2] 
    # 
    # # Upper bound (PP, alpha = 0)
    # lub_PP2_mat[ii,b] = compute_log_UBMarkov_BeBePois( Rmax, alpha, c_mle, gamma_mle, n, alfa )
    # Param. estimation (3 params PP)
    start_params <- c(alpha = 0.1, gamma= 1, u = 1)
    fit <- optim(par = start_params, fn = llik_PP3Parm, 
                 method = "L-BFGS-B",
                 lower = c(1e-16, 1e-16, 1e-16), 
                 upper = c(1-1e-10, Inf, Inf)) 
    alpha_mle = fit$par[1]
    gamma_mle = fit$par[2]
    c_mle     = fit$par[3] - alpha_mle
    
    # Upper bound (3 params PP)
    lub_PP3_mat[ii,b] = compute_log_UBMarkov_BeBePois( Rmax, alpha_mle, c_mle, gamma_mle, n, alfa )
    # Param. estimation (Mixed Poisson)
    start_params <- c(alpha = 0.1, u = 1, mu_gamma = 1)
    fit <- optim(par = start_params, fn = llik_MixPois,
                 method = "L-BFGS-B",
                 lower = c(1e-16, 1e-16, 1e-16), upper = c(1-1e-10, Inf, Inf))
    alpha_mle = fit$par[1]
    c_mle = fit$par[2] - alpha_mle
    mugamma_mle = fit$par[3]
    
    gamma_hyperparams = gamma_shape_rate(mugamma_mle,var_gamma)
    u = gamma_hyperparams$shape
    v = gamma_hyperparams$rate
    # Upper bound (Mixed Poisson)
    lub_MixPois_mat[ii,b] =  compute_log_UBMarkov_BeBeMixPois( Rmax, alpha_mle, c_mle,
                                                               n, Kn, u, v, alfa)
    # Param. estimation (Mixed Binomial)
    start_params <- c(a = 1, b = 1, mu_nb = 1)
    fit <- optim(par = start_params, fn = llik_MixBin,
                 method = "L-BFGS-B",
                 lower = c(1e-10, 1e-10, 1e-10), upper = c(Inf, Inf, var_nb-1e-10))
    a_mle = fit$par[1]
    b_mle = fit$par[2] 
    munb_mle = fit$par[3]
    
    nb_hyperparams = NegBin_params(munb_mle,var_nb)
    r_nb = nb_hyperparams$r
    p_nb = nb_hyperparams$p
    # Upper bound (Mixed Binomial)
    lub_MixBin_mat[ii,b] =  compute_log_UBMarkov_BeBeMixNBin( Rmax, a_mle, b_mle,
                                                              n, Kn, r_nb, p_nb, alfa)
  }
  
  # Frequentist estimators
  for(b in 1:Bfr){
    # Gen. data
    n_i = data[b,]
    idx_obs = which(n_i > 0)
    Kn = length(idx_obs)
    data_obs = n_i[idx_obs]
    Shat = sum( data_obs/n )
    
    # Minimization wrt beta
    start_params <- c(beta = alfa/2)
    fit <- optim(par = start_params, fn = lf_beta, 
                 method = "L-BFGS-B",
                 lower = 1e-16, upper = alfa - 1e-10 ) 
    beta_opt = fit$par[1]
    # Upper bound (freq)
    ub_Freq_mat[ii,b] = compute_UBFreq_BeBe(n, alfa, beta_opt, Shat)
  }
  pb$tick()
}


save_res = list("lub_PP3_mat" = lub_PP3_mat,
                "lub_MixPois_mat" = lub_MixPois_mat,
                "lub_MixBin_mat"  = lub_MixBin_mat,
                "ub_Freq_mat"  = ub_Freq_mat,
                "oracle_mat"  = oracle_mat,
                "oracle" = oracle)

file_name = paste0("../save/SimStudyFeatures_zipfs_",idx,".Rdat")
save(save_res, file = file_name)


