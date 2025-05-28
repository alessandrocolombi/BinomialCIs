setwd("C:/Users/colom/bnp_upperbounds/Rscripts/species")

# Librerie ----------------------------------------------------------------

suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
suppressWarnings(suppressPackageStartupMessages(library(progress)))
suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../../R/Rfunctions.R")

# Funzioni ----------------------------------------------------------------

worst_uniform <- function(m, n, r_unbounded, alpha){
  # Compute initial m_alpha
  m_a <- ceiling((n + r_unbounded - 1) / (r_unbounded - 1))
  
  if (m_a < m) {
    p <- rep(0, m)
    done <- FALSE
    CI_max <- 0
    
    while (!done) {
      p <- rep(0, m)
      p[1:m_a] <- 1 / m_a
      p <- p / sum(p)
      
      num_of_exp <- 100000
      real_vec_current <- numeric(num_of_exp)
      
      for (exp_ind in 1:num_of_exp) {
        c <- rmultinom(1, n, p)
        if (any(c == 0)) {
          real_vec_current[exp_ind] <- max(p[c == 0])
        }
      }
      
      CI_current <- quantile(real_vec_current, probs = 1 - alpha)
      
      if (CI_current > CI_max) {
        CI_max <- CI_current
        m_a <- m_a - 1
      } else {
        done <- TRUE
        m_a <- m_a + 1
        p <- rep(0, m)
        p[1:m_a] <- 1 / m_a
        p <- p / sum(p)
      }
    }
    
  } else {
    p <- rep(1 / m, m)
    CI_max <- NULL
  }
  
  return(list(p = p, m_a = m_a, CI_max = CI_max))
}
sim_worstunif = function(M,n,Rmax,alfa){
  runb = compute_r_unbounded(n=n, Rmax = Rmax, alfa = alfa)
  worst_uniform(m = M, n = n, r_unbounded = runb, alpha = alfa)$p
}

# Common parameters -------------------------------------------------------
alfa = 0.05
Bor  = 10000 # nel suo codice, painsky ripete 10k volte, noi 5k?
Bnoi = 10 
n = 1000
Rmax = 100; RmaxFD = 50
Mgrid = seq(10,1000,by = 50)
Nexp = length(Mgrid)
M_max = 500
seed = 42
set.seed(seed)
# Zipfs -------------------------------------------------------------------
s = 1.01

pain = rep( ub_pain(n = n, Rmax = Rmax, alfa = alfa), Nexp )
lub_Mrk_mat <- lub_Cnt_mat <- lub_FD_mat <- matrix(NA,Nexp,Bnoi)
oracle <- rep(NA,Nexp)

pb <- progress_bar$new(total = Nexp)
for(ii in 1:Nexp){
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
  }
  oracle[ii] = quantile(Mmax, probs = 1-alfa)
  
  # Markov - Cantelli - FD
  for(b in 1:Bnoi){
    # Gen. data
    data = sample(1:M, size = n, replace = TRUE, prob = ptrue)
    n_i = tabulate(data, nbins = M)
    idx_obs = which(n_i > 0)
    Kn = length(idx_obs)
    data_obs = n_i[idx_obs]
    
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
    lub_Cnt_mat[ii,b] = compute_log_UBCantelli( Rmax, alpha_mle, theta_mle, 
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
                "lub_Cnt_mat" = lub_Cnt_mat,
                "lub_FD_mat"  = lub_FD_mat,
                "oracle" = oracle)

save(save_res, file = "../save/SimStudyPain_zipfs.Rdat") 

ub_Mrk = exp(apply(lub_Mrk_mat, 1, quantile, probs = c(0.025,0.5,0.975)))
ub_Cnt = exp(apply(lub_Cnt_mat, 1, quantile, probs = c(0.025,0.5,0.975)))
ub_FD = exp(apply(lub_FD_mat, 1, quantile, probs = c(0.025,0.5,0.975)))

ymax = (11/10) * (max(oracle,pain,ub_Mrk[2,],ub_Cnt[2,],ub_FD[2,]))
ymin = 0
ylabs = round(seq(0,max(oracle,pain),length.out = 5),3)

par( mfrow = c(1,1), mar = c(3,3,1,0.5), mgp=c(1.5,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "M", ylab = " ",
     xlim = c(0,max(Mgrid) ) , ylim = c(ymin,ymax), 
     main = paste0("Zipfs"),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
points( x = Mgrid, y = oracle, 
        type = "l", lwd = 2, 
        col = "black" ) 
points( x = Mgrid, y = pain, 
        type = "l", lwd = 2, 
        col = "blue" ) 
points( x = Mgrid, y = ub_Mrk[2,], 
        type = "l", lwd = 2, 
        col = "darkred" ) 
points( x = Mgrid, y = ub_Cnt[2,], 
        type = "l", lwd = 2, 
        col = "darkgreen" ) 
points( x = Mgrid, y = ub_FD[2,], 
        type = "l", lwd = 2, 
        col = "darkorange" ) 
legend("bottomright",c("Oracle","Painsky","Markov","Cantelli","FD"), 
       lwd = 3, col = c("black","blue","darkred","darkgreen","darkorange"))

# Geometric -------------------------------------------------------------------
a = 0.4

pain = rep( ub_pain(n = n, Rmax = Rmax, alfa = alfa), Nexp )
lub_Mrk_mat <- lub_Cnt_mat <- lub_FD_mat <- matrix(NA,Nexp,Bnoi)
oracle <- rep(NA,Nexp)

pb <- progress_bar$new(total = Nexp)
for(ii in 1:Nexp){
  M = Mgrid[ii]
  # Gen. true distribution
  ptrue = sim_geom(M = M, a = a)
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
  }
  oracle[ii] = quantile(Mmax, probs = 1-alfa)
  
  # Markov - Cantelli - FD
  for(b in 1:Bnoi){
    # Gen. data
    data = sample(1:M, size = n, replace = TRUE, prob = ptrue)
    n_i = tabulate(data, nbins = M)
    idx_obs = which(n_i > 0)
    Kn = length(idx_obs)
    data_obs = n_i[idx_obs]
    
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
    lub_Cnt_mat[ii,b] = compute_log_UBCantelli( Rmax, alpha_mle, theta_mle, 
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
                "lub_Cnt_mat" = lub_Cnt_mat,
                "lub_FD_mat"  = lub_FD_mat,
                "oracle" = oracle)
save(save_res, file = "save/SimStudyPain_geom.Rdat")


ub_Mrk = exp(apply(lub_Mrk_mat, 1, quantile, probs = c(0.025,0.5,0.975)))
ub_Cnt = exp(apply(lub_Cnt_mat, 1, quantile, probs = c(0.025,0.5,0.975)))
ub_FD = exp(apply(lub_FD_mat, 1, quantile, probs = c(0.025,0.5,0.975)))

ymax = (11/10) * (max(oracle,pain,ub_Mrk[2,],ub_Cnt[2,],ub_FD[2,]))
ymin = 0
ylabs = round(seq(0,max(oracle,pain),length.out = 5),3)

par( mfrow = c(1,1), mar = c(3,3,1,0.5), mgp=c(1.5,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "M", ylab = " ",
     xlim = c(0,max(Mgrid) ) , ylim = c(ymin,ymax), 
     main = paste0("Geometric"),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
points( x = Mgrid, y = oracle, 
        type = "l", lwd = 2, 
        col = "black" ) 
points( x = Mgrid, y = pain, 
        type = "l", lwd = 2, 
        col = "blue" ) 
points( x = Mgrid, y = ub_Mrk[2,], 
        type = "l", lwd = 2, 
        col = "darkred" ) 
points( x = Mgrid, y = ub_Cnt[2,], 
        type = "l", lwd = 2, 
        col = "darkgreen" ) 
points( x = Mgrid, y = ub_FD[2,], 
        type = "l", lwd = 2, 
        col = "darkorange" ) 
legend("bottomright",c("Oracle","Painsky","Markov","Cantelli","FD"), 
       lwd = 3, col = c("black","blue","darkred","darkgreen","darkorange"))

# Negative binomial -------------------------------------------------------------------
l = 1; r = 0.003
pain = rep( ub_pain(n = n, Rmax = Rmax, alfa = alfa), Nexp )
lub_Mrk_mat <- lub_Cnt_mat <- lub_FD_mat <- matrix(NA,Nexp,Bnoi)
oracle <- rep(NA,Nexp)

pb <- progress_bar$new(total = Nexp)
for(ii in 1:Nexp){
  M = Mgrid[ii]
  # Gen. true distribution
  ptrue = sim_negbin(M = M, l = l, r = r)
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
  }
  oracle[ii] = quantile(Mmax, probs = 1-alfa)
  
  # Markov - Cantelli - FD
  for(b in 1:Bnoi){
    # Gen. data
    data = sample(1:M, size = n, replace = TRUE, prob = ptrue)
    n_i = tabulate(data, nbins = M)
    idx_obs = which(n_i > 0)
    Kn = length(idx_obs)
    data_obs = n_i[idx_obs]
    
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
    lub_Cnt_mat[ii,b] = compute_log_UBCantelli( Rmax, alpha_mle, theta_mle, 
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
                "lub_Cnt_mat" = lub_Cnt_mat,
                "lub_FD_mat"  = lub_FD_mat,
                "oracle" = oracle)
save(save_res, file = "save/SimStudyPain_negbin.Rdat")

ub_Mrk = exp(apply(lub_Mrk_mat, 1, quantile, probs = c(0.025,0.5,0.975)))
ub_Cnt = exp(apply(lub_Cnt_mat, 1, quantile, probs = c(0.025,0.5,0.975)))
ub_FD = exp(apply(lub_FD_mat, 1, quantile, probs = c(0.025,0.5,0.975)))

ymax = (11/10) * (max(oracle,pain,ub_Mrk[2,],ub_Cnt[2,],ub_FD[2,]))
ymin = 0
ylabs = round(seq(0,max(oracle,pain),length.out = 5),3)

par( mfrow = c(1,1), mar = c(3,3,1,0.5), mgp=c(1.5,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "M", ylab = " ",
     xlim = c(0,max(Mgrid) ) , ylim = c(ymin,ymax), 
     main = paste0("Neg. Binomial"),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
points( x = Mgrid, y = oracle, 
        type = "l", lwd = 2, 
        col = "black" ) 
points( x = Mgrid, y = pain, 
        type = "l", lwd = 2, 
        col = "blue" ) 
points( x = Mgrid, y = ub_Mrk[2,], 
        type = "l", lwd = 2, 
        col = "darkred" ) 
points( x = Mgrid, y = ub_Cnt[2,], 
        type = "l", lwd = 2, 
        col = "darkgreen" ) 
points( x = Mgrid, y = ub_FD[2,], 
        type = "l", lwd = 2, 
        col = "darkorange" ) 
legend("bottomright",c("Oracle","Painsky","Markov","Cantelli","FD"), 
       lwd = 3, col = c("black","blue","darkred","darkgreen","darkorange"))

# Uniform -------------------------------------------------------------------
pain = rep( ub_pain(n = n, Rmax = Rmax, alfa = alfa), Nexp )
lub_Mrk_mat <- lub_Cnt_mat <- lub_FD_mat <- matrix(NA,Nexp,Bnoi)
oracle <- rep(NA,Nexp)

pb <- progress_bar$new(total = Nexp)
for(ii in 1:Nexp){
  M = Mgrid[ii]
  # Gen. true distribution
  ptrue = sim_unif(M = M)
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
  }
  oracle[ii] = quantile(Mmax, probs = 1-alfa)
  
  # Markov - Cantelli - FD
  for(b in 1:Bnoi){
    # Gen. data
    data = sample(1:M, size = n, replace = TRUE, prob = ptrue)
    n_i = tabulate(data, nbins = M)
    idx_obs = which(n_i > 0)
    Kn = length(idx_obs)
    data_obs = n_i[idx_obs]
    
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
    lub_Cnt_mat[ii,b] = compute_log_UBCantelli( Rmax, alpha_mle, theta_mle, 
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
                "lub_Cnt_mat" = lub_Cnt_mat,
                "lub_FD_mat"  = lub_FD_mat,
                "oracle" = oracle)
save(save_res, file = "save/SimStudyPain_unif.Rdat")


ub_Mrk = exp(apply(lub_Mrk_mat, 1, quantile, probs = c(0.025,0.5,0.975)))
ub_Cnt = exp(apply(lub_Cnt_mat, 1, quantile, probs = c(0.025,0.5,0.975)))
ub_FD = exp(apply(lub_FD_mat, 1, quantile, probs = c(0.025,0.5,0.975)))

ymax = (11/10) * (max(oracle,pain,ub_Mrk[2,],ub_Cnt[2,],ub_FD[2,]))
ymin = 0
ylabs = round(seq(0,max(oracle,pain),length.out = 5),3)

par( mfrow = c(1,1), mar = c(3,3,1,0.5), mgp=c(1.5,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "M", ylab = " ",
     xlim = c(0,max(Mgrid) ) , ylim = c(ymin,ymax), 
     main = paste0("Uniform"),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
points( x = Mgrid, y = oracle, 
        type = "l", lwd = 2, 
        col = "black" ) 
points( x = Mgrid, y = pain, 
        type = "l", lwd = 2, 
        col = "blue" ) 
points( x = Mgrid, y = ub_Mrk[2,], 
        type = "l", lwd = 2, 
        col = "darkred" ) 
points( x = Mgrid, y = ub_Cnt[2,], 
        type = "l", lwd = 2, 
        col = "darkgreen" ) 
points( x = Mgrid, y = ub_FD[2,], 
        type = "l", lwd = 2, 
        col = "darkorange" ) 
legend("bottomright",c("Oracle","Painsky","Markov","Cantelli","FD"), 
       lwd = 3, col = c("black","blue","darkred","darkgreen","darkorange"))

# Worst Unif. -------------------------------------------------------------------
pain = rep( ub_pain(n = n, Rmax = Rmax, alfa = alfa), Nexp )
oracle <- rep(NA,Nexp)

pb <- progress_bar$new(total = Nexp)
for(ii in 1:Nexp){
  M = Mgrid[ii]
  # Gen. true distribution
  ptrue = sim_worstunif(M = M, n = n, Rmax = Rmax, alfa = alfa)
  # plot(ptrue)
  Mmax = rep(NA,Bor)
  for(b in 1:Bor){
    # Gen. data
    data = sample(1:M, size = n, replace = TRUE, prob = ptrue)
    n_i = tabulate(data, nbins = M)
    idx_obs = which(n_i > 0)
    Kn = length(idx_obs)
    if(Kn == M){
      Mmax[b] = 0
    }else{
      idx_unobs = which(n_i == 0)
      Mmax[b] = max(ptrue[idx_unobs])
    }
  }
  oracle[ii] = quantile(Mmax, probs = 1-alfa)
  pb$tick()
}

ymax = (11/10) * (max(oracle,pain))
ymin = 0
# ylabs = round(as.numeric(quantile(oracle, probs = c(0,0.25,0.5,0.75,1))), 4)
ylabs = round(seq(0,max(oracle,pain),length.out = 5),3)

par( mfrow = c(1,1), mar = c(3,3,1,0.5), mgp=c(1.5,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "M", ylab = " ",
     xlim = c(0,max(Mgrid) ) , ylim = c(ymin,ymax), 
     main = paste0("Worst uniform"),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
points( x = Mgrid, y = oracle, 
        type = "l", lwd = 2, 
        col = "black" ) 
points( x = Mgrid, y = pain, 
        type = "l", lwd = 2, 
        col = "blue" ) 

# Brutta ------------------------------------------------------------------
stop("FERMO IO, qua è la brutta")
a = 2; b = 2
M = 10
pex = sim_betabin(M,a,b)
plot(x = 0:(M-1), pex)

Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
M = 100
n = 1000 
s = 1.01
ptrue = sim_zipfs(M = M, s = s)
plot(ptrue)
data = sample(1:M, size = n, replace = TRUE, prob = ptrue)
n_i = tabulate(data, nbins = M)
idx_obs = which(n_i > 0)
Kn = length(idx_obs)
data_obs = n_i[idx_obs]

alpha = 0; theta = 10
log_eppfPYP( n, Kn, data_obs, alpha, theta )



# Define the log-likelihood function (you must replace this with your actual function)
llik <- function(x) {
  alpha <- x[1]
  theta <- x[2]
  -log_eppfPYP( n, Kn, data_obs, alpha, theta )
}
# Starting values 
start_params <- c(alpha = 0.1, theta = 1)
# Use optim to maximise the log-likelihood
fit <- optim(par = start_params, fn = llik, method = "L-BFGS-B",
             lower = c(0, -1), upper = c(1, Inf))  # Set bounds if needed
fit

alpha_mle = fit$par[1]
theta_mle = fit$par[2]

Nrep = 10; Natoms = 100
sim_SB =  r_SB( Nrep, Natoms, alpha_mle, theta_mle, 42)

plot(0,0,type = "n", xlim = c(0,100), ylim = c(0,0.2))
for(ii in 1:Nrep){
  points(x = 1:Natoms, y = sort(sim_SB[ii,], decreasing = TRUE), type = "b", pch = 16, lwd = 2)  
}
points(x = 1:Natoms, y = ptrue, type = "p", pch = 16, col = "red")


# Geometric -------------------------------------------------------------------
a = 0.4

pain = rep( ub_pain(n = n, Rmax = Rmax, alfa = alfa), Nexp )
oracle <- rep(NA,Nexp)

pb <- progress_bar$new(total = Nexp)
for(ii in 1:Nexp){
  M = Mgrid[ii]
  # Gen. true distribution
  ptrue = sim_geom(M = M, a = a)
  # plot(ptrue)
  Mmax = rep(NA,B)
  for(b in 1:B){
    # Gen. data
    data = sample(1:M, size = n, replace = TRUE, prob = ptrue)
    n_i = tabulate(data, nbins = M)
    idx_obs = which(n_i > 0)
    Kn = length(idx_obs)
    if(Kn == M){
      Mmax[b] = 0
    }else{
      idx_unobs = which(n_i == 0)
      Mmax[b] = max(ptrue[idx_unobs])
    }
  }
  oracle[ii] = quantile(Mmax, probs = 1-alfa)
  pb$tick()
}

ymax = (11/10) * (max(oracle,pain))
ymin = 0
# ylabs = round(as.numeric(quantile(oracle, probs = c(0,0.25,0.5,0.75,1))), 4)
ylabs = round(seq(0,max(oracle,pain),length.out = 5),3)

par( mfrow = c(1,1), mar = c(3,3,1,0.5), mgp=c(1.5,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "M", ylab = " ",
     xlim = c(0,max(Mgrid) ) , ylim = c(ymin,ymax), 
     main = paste0("Geometric"),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
points( x = Mgrid, y = oracle, 
        type = "l", lwd = 2, 
        col = "black" ) 
points( x = Mgrid, y = pain, 
        type = "l", lwd = 2, 
        col = "blue" ) 
# Negative binomial -------------------------------------------------------------------
l = 1; r = 0.003

pain = rep( ub_pain(n = n, Rmax = Rmax, alfa = alfa), Nexp )
oracle <- rep(NA,Nexp)

pb <- progress_bar$new(total = Nexp)
for(ii in 1:Nexp){
  M = Mgrid[ii]
  # Gen. true distribution
  ptrue = sim_negbin(M = M, l = l, r = r)
  # plot(ptrue)
  Mmax = rep(NA,B)
  for(b in 1:B){
    # Gen. data
    data = sample(1:M, size = n, replace = TRUE, prob = ptrue)
    n_i = tabulate(data, nbins = M)
    idx_obs = which(n_i > 0)
    Kn = length(idx_obs)
    if(Kn == M){
      Mmax[b] = 0
    }else{
      idx_unobs = which(n_i == 0)
      Mmax[b] = max(ptrue[idx_unobs])
    }
  }
  oracle[ii] = quantile(Mmax, probs = 1-alfa)
  pb$tick()
}

ymax = (11/10) * (max(oracle,pain))
ymin = 0
# ylabs = round(as.numeric(quantile(oracle, probs = c(0,0.25,0.5,0.75,1))), 4)
ylabs = round(seq(0,max(oracle,pain),length.out = 5),3)

par( mfrow = c(1,1), mar = c(3,3,1,0.5), mgp=c(1.5,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "M", ylab = " ",
     xlim = c(0,max(Mgrid) ) , ylim = c(ymin,ymax), 
     main = paste0("Neg. Binomial"),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
points( x = Mgrid, y = oracle, 
        type = "l", lwd = 2, 
        col = "black" ) 
points( x = Mgrid, y = pain, 
        type = "l", lwd = 2, 
        col = "blue" ) 
# Beta binomial -------------------------------------------------------------------
a = 2; b = 2

pain = rep( ub_pain(n = n, Rmax = Rmax, alfa = alfa), Nexp )
oracle <- rep(NA,Nexp)

# pdf("prova.pdf")
pb <- progress_bar$new(total = Nexp)
for(ii in 1:Nexp){
  M = Mgrid[ii]
  # Gen. true distribution
  ptrue = sim_betabin(M = M, a = a, b = b)
  plot(ptrue, main = paste0("M = ",M))
  Mmax = rep(NA,B)
  for(b in 1:B){
    # Gen. data
    data = sample(1:M, size = n, replace = TRUE, prob = ptrue)
    n_i = tabulate(data, nbins = M)
    idx_obs = which(n_i > 0)
    Kn = length(idx_obs)
    if(Kn == M){
      Mmax[b] = 0
    }else{
      idx_unobs = which(n_i == 0)
      Mmax[b] = max(ptrue[idx_unobs])
    }
  }
  oracle[ii] = quantile(Mmax, probs = 1-alfa)
  pb$tick()
}
# dev.off()

ymax = (11/10) * (max(oracle,pain))
ymin = 0
# ylabs = round(as.numeric(quantile(oracle, probs = c(0,0.25,0.5,0.75,1))), 4)
ylabs = round(seq(0,max(oracle,pain),length.out = 5),3)

par( mfrow = c(1,1), mar = c(3,3,1,0.5), mgp=c(1.5,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "M", ylab = " ",
     xlim = c(0,max(Mgrid) ) , ylim = c(ymin,ymax), 
     main = paste0("Beta binomial"),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
points( x = Mgrid, y = oracle, 
        type = "l", lwd = 2, 
        col = "black" ) 
points( x = Mgrid, y = pain, 
        type = "l", lwd = 2, 
        col = "blue" ) 
