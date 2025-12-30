setwd("C:/Users/colom/BinomialCIs/Rscripts/SeedsData")

# Librerie ----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
suppressWarnings(suppressPackageStartupMessages(library(progress)))
suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../../R/Rfunctions.R")

## ------------------------------------------------------------
## Packages and data
## ------------------------------------------------------------
library(future.apply)

## ------------------------------------------------------------
## 1. Load functions
## ------------------------------------------------------------
eps=0.01;alpha= 0.05;C_target= 0.95;g_Chao2009= 0.95;beta= 1e-5
stopping_rule_run <- function(data, seed,
                              eps            = 0.01,
                              alpha          = 0.05,
                              C_target       = 0.95,
                              g_Chao2009     = 0.95,
                              beta           = 1e-5)
{
  ## Functions
  suppressWarnings(suppressPackageStartupMessages(library(tibble)))
  suppressWarnings(suppressPackageStartupMessages(library(parallel)))
  suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
  suppressWarnings(suppressPackageStartupMessages(library(progress)))
  suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
  Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
  source("../../R/Rfunctions.R")
  
  
  set.seed(seed)
  n = nrow(data) # total num. obs.
  Kn = ncol(data) # total num. distinct
  ordered_idx = sample(1:n, size = n) # choose ordering of obs.
  
  # Stopping flags and outputs
  stopped_bounded   <- FALSE
  stopped_unbounded <- FALSE
  stopped_cov       <- FALSE
  stopped_Chao2009  <- FALSE
  
  Nstop_bounded   <- NA_integer_
  Nstop_unbounded <- NA_integer_
  Nstop_cov       <- NA_integer_
  Nstop_Chao2009  <- NA_integer_
  
  ## ------------------------------------------------------------
  ## Run loop up to n_max = n
  ## ------------------------------------------------------------
  n_max = n
  ni = 2
  for(ni in 2:(n_max-1)) {
    # Allow for a non-multiple n_max if needed
    remaining <- n_max - ni
    if (remaining <= 0L) break
    
    ## ---- Observed abundance vector (true + error species) ----
    idx_species_i = ordered_idx[1:ni] # select obs. up to time ni
    data_i = matrix(data[idx_species_i,], nrow = ni, ncol = Kn) # get data
    Nj_i = colSums(data_i) # compute frequencies
    Nj_i = Nj_i[Nj_i > 0]
    Kobs_i = length(which(Nj_i > 0))
    if( Kobs_i == 0L) next   # nothing observed yet
    
    ## ---- 5.1 Bounded-alphabet CI rule (on observed data) ----
    if (!stopped_bounded) {
      b_n <- log(ni)
      Mguess = 10 * Kobs_i
      Nj_guess = c(Nj_i, rep(0,Mguess - length(Nj_i) ))
      U_bounded <- compute_UB_analytical(ni, Nj_guess, Mguess, b_n, alpha, FALSE)
      U_bounded <- min(1,U_bounded)
      if (!is.na(U_bounded) && U_bounded < eps) {
        stopped_bounded <- TRUE
        Nstop_bounded   <- ni
      }
    }
    
    ## ---- 5.2 Unbounded-alphabet CI rule (on observed data) ----
    if (!stopped_unbounded) {
      Shat  <- sum(Nj_i) / ni
      Sstar <- ( sqrt( -log(beta) / (2 * ni) ) +
                   sqrt( Shat + (-log(beta) / (2 * ni)) ) )^2
      r_n   <- log( Sstar / (-log(1 - alpha + beta)) ) +
               log(ni) - log(log(ni))
      U_unbounded <- compute_UB_rnorm(ni, alpha, beta, r_n, Shat)
      
      if (!is.na(U_unbounded) && U_unbounded <= eps) {
        stopped_unbounded <- TRUE
        Nstop_unbounded   <- ni
      }
    }
    
    ## ---- 5.3 Chao–Jost coverage-based rule (same counts_obs) ----
    if (!stopped_cov) {
      C_hat <- coverage_ChaoJost(ni, Nj_i)
      C_hat = max(0,C_hat)
      if (!is.na(C_hat) && C_hat >= C_target) {
        stopped_cov <- TRUE
        Nstop_cov   <- ni
      }
    }
    
    ## ---- 5.4 Chao2009 Eq.15 stopping rule (same counts_obs) ----
    if (!stopped_Chao2009) {
      mg <- Nsample_Chao2009(ni, Nj_i, g_Chao2009)
      if (!is.na(mg) && mg <= 1) {
        stopped_Chao2009 <- TRUE
        Nstop_Chao2009   <- ni
      }
    }
    
    # Early exit if all four rules have stopped
    if (stopped_bounded && stopped_unbounded && stopped_cov && stopped_Chao2009) break
  }
  
  ## ------------------------------------------------------------
  ## Post-loop: handle rules that *never* stopped by n_max
  ## ------------------------------------------------------------
  if (!stopped_bounded) {
    Nstop_bounded          <- n_max
  }
  if (!stopped_unbounded) {
    Nstop_unbounded        <- n_max
  }
  if (!stopped_cov) {
    Nstop_cov              <- n_max
  }
  if (!stopped_Chao2009) {
    Nstop_Chao2009         <- n_max
  }
  
  return(list("Nstop_bounded" = Nstop_bounded,
              "Nstop_unbounded" = Nstop_unbounded,
              "Nstop_cov" = Nstop_cov,
              "Nstop_Chao2009" = Nstop_Chao2009))
}
  

# Function to run simulation
Nrep = 4; num_cores = 2; seed0 = 123
stopping_rule = function( data, Nrep, num_cores, seed0,
                          eps            = 0.01,
                          alpha          = 0.05,
                          C_target       = 0.95,
                          g_Chao2009     = 0.95,
                          beta           = 1e-5)
{
  res_mat = matrix(NA,nrow = Nrep, ncol = 4)
  colnames(res_mat) = c("Nstop_bounded","Nstop_unbounded","Nstop_cov","Nstop_Chao2009")
  
  ## Parallel run (no prints allowed)
  seeds = sample(1:999999, size = Nrep)
  cluster <- makeCluster(num_cores, type = "SOCK")
  doSNOW::registerDoSNOW(cluster)
  clusterExport(cluster, list("compute_UB_analytical"),
                envir = environment())
  inner_result = parLapply( cl = cluster, x = seeds,
                            fun = stopping_rule_run,
                            data = data,
                            eps = eps, alpha = alpha, 
                            C_target = C_target, g_Chao2009 = g_Chao2009,
                            beta = beta)
  stopCluster(cluster)
  
  for(hh in 1:4){
    res_mat[,hh] = sapply(inner_result, function(x) x[[hh]])
  }
  
  return(res_mat)
}


## ------------------------------------------------------------
## 2. Run
## ------------------------------------------------------------
data(DiversityData)
data = t(DiversityData$Inci_raw)
n = nrow(data)
Nj = c(apply(data, 2, sum))
Kobs = length( which(Nj > 0) )
TabNj = tabulate(Nj, nbins = n) 

eps            = 0.1
alpha          = 0.05
C_target       = 0.95
g_Chao2009     = 0.95
beta           = 1e-5

Nrep = 100
seed0 = 4224
num_cores = 5
res = stopping_rule( data, Nrep, num_cores, seed0, 
                     eps, alpha, C_target, g_Chao2009, beta)
res
n



colMeans(res)















