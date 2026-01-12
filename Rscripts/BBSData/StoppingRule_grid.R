wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100)
choose_wd = wd_vec[1] # <--- modify here
wd = paste0(choose_wd,"BinomialCIs/Rscripts/BBSData/")
setwd(wd)


# Librerie ----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
suppressWarnings(suppressPackageStartupMessages(library(progress)))
# suppressWarnings(suppressPackageStartupMessages(library(VGAM)))

Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../../R/Rfunctions.R")

## ------------------------------------------------------------
## 1. Load functions
## 7------------------------------------------------------------

## Bounded and Unbounded UB stopping rules
SR_grid_multiple_run <- function(eps, data, seed0, Nrep, alpha, beta)
{
  ## Functions
  suppressWarnings(suppressPackageStartupMessages(library(tibble)))
  suppressWarnings(suppressPackageStartupMessages(library(parallel)))
  suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
  suppressWarnings(suppressPackageStartupMessages(library(progress)))
  suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
  Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
  source("../../R/Rfunctions.R")
  
  set.seed(seed0)
  seeds = sample(1:999999, size = Nrep)
  res = matrix(nrow = Nrep, ncol = 2)
  colnames(res) = c("Bounded","Unbounded")
  
  for(ii in 1:Nrep){
    seed = seeds[ii]
    res[ii,] = SR_grid_single_run(eps, data, seed, alpha, beta)
  }
  return(res)
}

SR_grid_single_run <- function(eps, data, seed, alpha, beta)
{
  ## Functions
  # suppressWarnings(suppressPackageStartupMessages(library(tibble)))
  # suppressWarnings(suppressPackageStartupMessages(library(parallel)))
  # suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
  # suppressWarnings(suppressPackageStartupMessages(library(progress)))
  # suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
  # Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
  # source("../../R/Rfunctions.R")
  
  
  set.seed(seed) # set seed
  
  n = nrow(data) # total num. obs.
  Kn = ncol(data) # total num. distinct
  ordered_idx = sample(1:n, size = n) # choose ordering of obs.
  
  # Stopping flags and outputs
  stopped_bounded   <- FALSE
  stopped_unbounded <- FALSE
  Nstop_bounded   <- NA_integer_
  Nstop_unbounded <- NA_integer_
  
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
    # data_i = matrix(data[idx_species_i,], nrow = ni, ncol = Kn) # get data
    data_i = data[idx_species_i,] 
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
    
    # Early exit if all four rules have stopped
    if (stopped_bounded && stopped_unbounded) break
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
  
  return(c(Nstop_bounded,Nstop_unbounded))
}


SR_grid = function( eps_grid, data, Nrep, num_cores, seed0,
                    alpha = 0.05, beta = 1e-5)
{
  Lgrid = length(eps_grid) # grid length
  res_list = vector("list",Lgrid)
  res_list = lapply(res_list, function(x) {
    y = matrix(nrow = Nrep, ncol = 2)
    colnames(y) = c("Bounded","Unbounded")
    y
    }  )


  ## Parallel run (no prints allowed)
  cluster <- makeCluster(num_cores, type = "SOCK")
  doSNOW::registerDoSNOW(cluster)
  clusterExport(cluster, list("SR_grid_single_run"),
                envir = environment())
  res_list = parLapply( cl = cluster, x = eps_grid,
                        fun = SR_grid_multiple_run,
                        data = data, alpha = alpha, beta = beta,
                        seed0 = seed0, Nrep = Nrep)
  stopCluster(cluster)

  return(res_list)
}


## Coverages stopping rules
SRcov_grid_multiple_run <- function(cov, data, seed0, Nrep)
{
  ## Functions
  suppressWarnings(suppressPackageStartupMessages(library(tibble)))
  suppressWarnings(suppressPackageStartupMessages(library(parallel)))
  suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
  suppressWarnings(suppressPackageStartupMessages(library(progress)))
  suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
  Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
  source("../../R/Rfunctions.R")
  
  set.seed(seed0)
  seeds = sample(1:999999, size = Nrep)
  res = matrix(nrow = Nrep, ncol = 2)
  colnames(res) = c("Coverage","Chao2009")
  
  for(ii in 1:Nrep){
    seed = seeds[ii]
    res[ii,] = SRcov_grid_single_run(cov, data, seed)
  }
  return(res)
}

SRcov_grid_single_run <- function(cov, data, seed)
{
  ## Functions
  # suppressWarnings(suppressPackageStartupMessages(library(tibble)))
  # suppressWarnings(suppressPackageStartupMessages(library(parallel)))
  # suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
  # suppressWarnings(suppressPackageStartupMessages(library(progress)))
  # suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
  # Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
  # source("../../R/Rfunctions.R")
  
  
  set.seed(seed) # set seed
  
  n = nrow(data) # total num. obs.
  Kn = ncol(data) # total num. distinct
  ordered_idx = sample(1:n, size = n) # choose ordering of obs.
  
  # Stopping flags and outputs
  stopped_cov       <- FALSE
  stopped_Chao2009  <- FALSE
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
    # data_i = matrix(data[idx_species_i,], nrow = ni, ncol = Kn) # get data
    data_i = data[idx_species_i,] 
    Nj_i = colSums(data_i) # compute frequencies
    Nj_i = Nj_i[Nj_i > 0]
    Kobs_i = length(which(Nj_i > 0))
    if( Kobs_i == 0L) next   # nothing observed yet
    
    ## ---- 5.3 Chao–Jost coverage-based rule (same counts_obs) ----
    if (!stopped_cov) {
      C_hat <- coverage_ChaoJost(ni, Nj_i)
      C_hat = max(0,C_hat)
      if (!is.na(C_hat) && C_hat >= cov) {
        stopped_cov <- TRUE
        Nstop_cov   <- ni
      }
    }
    
    ## ---- 5.4 Chao2009 Eq.15 stopping rule (same counts_obs) ----
    if (!stopped_Chao2009) {
      mg <- Nsample_Chao2009(ni, Nj_i, cov)
      if (!is.na(mg) && mg <= 1) {
        stopped_Chao2009 <- TRUE
        Nstop_Chao2009   <- ni
      }
    }
    
    # Early exit if all four rules have stopped
    if (stopped_cov && stopped_Chao2009) break
  }
  
  ## ------------------------------------------------------------
  ## Post-loop: handle rules that *never* stopped by n_max
  ## ------------------------------------------------------------
  if (!stopped_cov) {
    Nstop_cov <- n_max
  }
  if (!stopped_Chao2009) {
    Nstop_Chao2009 <- n_max
  }
  
  return(c(Nstop_cov,Nstop_Chao2009))
}


SRcov_grid = function( cov_grid, data, Nrep, num_cores, seed0)
{
  Lgrid = length(cov_grid) # grid length
  res_list = vector("list",Lgrid)
  res_list = lapply(res_list, function(x) {
    y = matrix(nrow = Nrep, ncol = 2)
    colnames(y) = c("Coverage","Chao2009")
    y
  }  )
  
  
  ## Parallel run (no prints allowed)
  cluster <- makeCluster(num_cores, type = "SOCK")
  doSNOW::registerDoSNOW(cluster)
  clusterExport(cluster, list("SRcov_grid_single_run"),
                envir = environment())
  res_list = parLapply( cl = cluster, x = cov_grid,
                        fun = SRcov_grid_multiple_run,
                        data = data, seed0 = seed0, Nrep = Nrep)
  stopCluster(cluster)
  
  return(res_list)
}


## ------------------------------------------------------------
## 2. Run - Bounded and Unbounded
## ------------------------------------------------------------
load(paste0("data/Data2019_allRoutes.Rdat"))
data = incidence_matrix
n = nrow(data)

eps_grid = seq(0.001,0.1,length.out = 34*2 )
alpha = 0.05
beta = 1e-5

Nrep = 10
seed0 = 4224
num_cores = 34

# res = SR_grid( eps_grid, data, Nrep, num_cores, seed0, alpha, beta)
# save(res, file = "save/Data2019_allroutes_epsgrid.Rdat")

## ------------------------------------------------------------
## 3. Run - Coverage
## ------------------------------------------------------------
load(paste0("data/Data2019_allRoutes.Rdat"))
data = incidence_matrix
n = nrow(data)

cov_grid = 1 - eps_grid

Nrep = 20
seed0 = 4224
num_cores = 34

# res_cov = SRcov_grid( cov_grid, data, Nrep, num_cores, seed0)
# save(res_cov, file = "save/Data2019_allroutes_covgrid.Rdat")

## ------------------------------------------------------------
## 3. Read and plot
## ------------------------------------------------------------
stop_here = TRUE
ltype = c(1,1,2,2)
mycol = c("darkgreen","darkorange","deeppink","lightblue")
ygrids = vector("list",4)
ygrids[[1]]<-ygrids[[2]]<-eps_grid
ygrids[[3]]<-ygrids[[4]]<-(1-cov_grid)

if(!stop_here){
  load("save/Data2019_allroutes_epsgrid.Rdat")
  load("save/Data2019_allroutes_covgrid.Rdat")
  
  res_list = lapply(res, function(x) apply(x,2,quantile,probs = 0.5));res_med = do.call(rbind,res_list)
  res_cov_list = lapply(res_cov, function(x) apply(x,2,quantile,probs = 0.5));res_cov_med = do.call(rbind,res_cov_list)
  res_all = cbind(res_med,res_cov_med)
  
  par(mfrow = c(1,1),bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,1,0))
  plot(0,0,type = "n", main = "", ylab = "Nstop",
       xlim = range(eps_grid), ylim = c(0,n), 
       xlab = paste0(expression(epsilon)," / 1 - coverage") )
  for(i in 1:4){
    points(y = res_all[,i], x = ygrids[[i]], 
           type = "l", lty = ltype[i], 
           lwd = 3, col = mycol[i], pch = 16)
  }
  legend("topright", c("Bounded","Unbounded","Coverage","Chao2009"), 
         col = mycol, lty = ltype, lwd = 3)
  
}
