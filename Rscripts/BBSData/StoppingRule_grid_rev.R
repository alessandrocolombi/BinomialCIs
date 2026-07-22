wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_bocconi = "/home/colombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100,wd_bocconi)
choose_wd = wd_vec[4] # <--- modify here
wd = paste0(choose_wd,"BinomialCIs/Rscripts/BBSData")
setwd(wd)

# Timestamp iniziale: utile per verificare dal log della VM che lo script sia
# partito nella working directory attesa.
script_start_time <- Sys.time()
cat(sprintf(
  "[%s] SCRIPT START | working directory: %s\n",
  format(script_start_time, "%Y-%m-%d %H:%M:%S"),
  getwd()
))
flush.console()

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
SR_grid_multiple_run <- function(eps, data, seed0, Nrep, alpha, beta, n_max)
{
  # Ogni worker gestisce un valore di eps e tutte le relative repliche.
  # Il PID permette di distinguere i messaggi prodotti in parallelo.
  worker_start_time <- Sys.time()
  worker_pid <- Sys.getpid()
  report_every <- max(1L, ceiling(Nrep / 5L))  # circa 5 update per eps
  cat(sprintf(
    "[%s] EPS START | pid=%d | eps=%.6f | repliche=%d\n",
    format(worker_start_time, "%Y-%m-%d %H:%M:%S"),
    worker_pid, eps, Nrep
  ))
  flush.console()

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
  names_methods = c("Bounded","Unbounded","Bounded_rev","Unbounded_rev")
  res = matrix(nrow = Nrep, ncol = length(names_methods))
  colnames(res) = names_methods
  
  ii = 1
  for(ii in 1:Nrep){
    seed = seeds[ii]
    res[ii,] = SR_grid_single_run(eps, data, seed, alpha, beta, n_max)

    # Aggiornamento periodico senza stampare ogni singola replica.
    if (ii %% report_every == 0L || ii == Nrep) {
      elapsed_seconds <- as.numeric(difftime(
        Sys.time(), worker_start_time, units = "secs"
      ))
      cat(sprintf(
        paste0(
          "[%s] EPS PROGRESS | pid=%d | eps=%.6f | ",
          "replica=%d/%d (%.0f%%) | elapsed=%.1f min\n"
        ),
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        worker_pid, eps, ii, Nrep, 100 * ii / Nrep,
        elapsed_seconds / 60
      ))
      flush.console()
    }
  }

  worker_elapsed <- as.numeric(difftime(
    Sys.time(), worker_start_time, units = "secs"
  ))
  cat(sprintf(
    "[%s] EPS DONE | pid=%d | eps=%.6f | elapsed=%.1f min\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    worker_pid, eps, worker_elapsed / 60
  ))
  flush.console()
  return(res)
}

SR_grid_single_run <- function(eps, data, seed, alpha, beta, Nmax)
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
  stopped_bounded_rev   <- FALSE
  stopped_unbounded_rev <- FALSE
  Nstop_bounded   <- NA_integer_
  Nstop_unbounded <- NA_integer_
  Nstop_bounded_rev   <- NA_integer_
  Nstop_unbounded_rev <- NA_integer_
  
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
    
    ## ---- 5.3 Bounded-alphabet CI rule - revision (on observed data) ----
    if (!stopped_bounded_rev) {
      b_n <- log(ni)
      Mguess = 10 * Kobs_i
      Nj_guess = c(Nj_i, rep(0,Mguess - length(Nj_i) ))
      alpha1 = 0.99*alpha
      delta = 0.01*alpha
      log_m_bar = compute_lmbar_bounded(ni, Nj_guess, Mguess, b_n, alpha, n_max)
      lhs <- log_m_bar + (ni - b_n)*log(1-eps) 
      rhs <- log(alpha - delta)
      if (!is.na(lhs) && lhs < rhs) {
        stopped_bounded_rev <- TRUE
        Nstop_bounded_rev   <- ni
      }
    }
    
    ## ---- 5.4 Unbounded-alphabet CI rule - revision (on observed data) ----
    if (!stopped_unbounded_rev) {
      Shat  <- sum(Nj_i) / ni
      Sbar <- ( sqrt( log(Nmax)-log(beta) / (2 * ni) ) + sqrt( Shat + (log(Nmax)-log(beta) / (2 * ni)) ) )^2
      lhs <- log(Sbar) + ni * log(1 - eps) - log(eps)
      rhs <- log(alpha - beta)
      
      if (!is.na(lhs) && lhs <= rhs) {
        stopped_unbounded_rev <- TRUE
        Nstop_unbounded_rev   <- ni
      }
    }
    
    # Early exit if all four rules have stopped
    if (stopped_bounded && stopped_unbounded && stopped_bounded_rev && stopped_unbounded_rev) break
  }
  
  ## ------------------------------------------------------------
  ## Post-loop: handle rules that *never* stopped by n_max
  ## ------------------------------------------------------------
  if (!stopped_bounded) {
    Nstop_bounded <- n_max
  }
  if (!stopped_unbounded) {
    Nstop_unbounded <- n_max
  }
  if (!stopped_bounded_rev) {
    Nstop_bounded_rev <- n_max
  }
  if (!stopped_unbounded_rev) {
    Nstop_unbounded_rev <- n_max
  }
  
  return(c(Nstop_bounded,Nstop_unbounded,Nstop_bounded_rev,Nstop_unbounded_rev))
}


SR_grid = function( eps_grid, data, Nrep, num_cores, seed0,
                    n_max, alpha = 0.05, beta = 1e-5)
{
  grid_start_time <- Sys.time()
  Lgrid = length(eps_grid) # grid length
  res_list = vector("list",Lgrid)
  res_list = lapply(res_list, function(x) {
    names_methods = c("Bounded","Unbounded","Bounded_rev","Unbounded_rev")
    y = matrix( nrow = Nrep, ncol = length(names_methods) )
    colnames(y) = names_methods
    y
  }  )
  
  
  ## Parallel run. outfile="" inoltra al log principale i messaggi dei worker.
  cat(sprintf(
    paste0(
      "[%s] EPS GRID START | grid points=%d | repliche/punto=%d | ",
      "total runs=%d | cores=%d\n"
    ),
    format(grid_start_time, "%Y-%m-%d %H:%M:%S"),
    Lgrid, Nrep, Lgrid * Nrep, num_cores
  ))
  flush.console()

  cluster <- makeCluster(num_cores, type = "SOCK", outfile = "")
  doSNOW::registerDoSNOW(cluster)
  clusterExport(cluster, list("SR_grid_single_run"),
                envir = environment())
  res_list = parLapply( cl = cluster, x = eps_grid,
                        fun = SR_grid_multiple_run,
                        data = data, alpha = alpha, beta = beta,
                        seed0 = seed0, Nrep = Nrep, n_max = n_max)
  stopCluster(cluster)

  grid_elapsed <- as.numeric(difftime(
    Sys.time(), grid_start_time, units = "secs"
  ))
  cat(sprintf(
    "[%s] EPS GRID DONE | grid points=%d | elapsed=%.1f min\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    Lgrid, grid_elapsed / 60
  ))
  flush.console()
  
  return(res_list)
}


## Coverages stopping rules
SRcov_grid_multiple_run <- function(cov, data, seed0, Nrep)
{
  # Monitoraggio analogo alla griglia epsilon, pronto anche se questa analisi
  # viene riattivata in futuro.
  worker_start_time <- Sys.time()
  worker_pid <- Sys.getpid()
  report_every <- max(1L, ceiling(Nrep / 5L))
  cat(sprintf(
    "[%s] COVERAGE START | pid=%d | target=%.6f | repliche=%d\n",
    format(worker_start_time, "%Y-%m-%d %H:%M:%S"),
    worker_pid, cov, Nrep
  ))
  flush.console()

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

    if (ii %% report_every == 0L || ii == Nrep) {
      elapsed_seconds <- as.numeric(difftime(
        Sys.time(), worker_start_time, units = "secs"
      ))
      cat(sprintf(
        paste0(
          "[%s] COVERAGE PROGRESS | pid=%d | target=%.6f | ",
          "replica=%d/%d (%.0f%%) | elapsed=%.1f min\n"
        ),
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        worker_pid, cov, ii, Nrep, 100 * ii / Nrep,
        elapsed_seconds / 60
      ))
      flush.console()
    }
  }

  worker_elapsed <- as.numeric(difftime(
    Sys.time(), worker_start_time, units = "secs"
  ))
  cat(sprintf(
    "[%s] COVERAGE DONE | pid=%d | target=%.6f | elapsed=%.1f min\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    worker_pid, cov, worker_elapsed / 60
  ))
  flush.console()
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
  grid_start_time <- Sys.time()
  Lgrid = length(cov_grid) # grid length
  res_list = vector("list",Lgrid)
  res_list = lapply(res_list, function(x) {
    y = matrix(nrow = Nrep, ncol = 2)
    colnames(y) = c("Coverage","Chao2009")
    y
  }  )
  
  
  ## Parallel run con inoltro dei messaggi dei worker al log principale.
  cat(sprintf(
    paste0(
      "[%s] COVERAGE GRID START | grid points=%d | repliche/punto=%d | ",
      "total runs=%d | cores=%d\n"
    ),
    format(grid_start_time, "%Y-%m-%d %H:%M:%S"),
    Lgrid, Nrep, Lgrid * Nrep, num_cores
  ))
  flush.console()

  cluster <- makeCluster(num_cores, type = "SOCK", outfile = "")
  doSNOW::registerDoSNOW(cluster)
  clusterExport(cluster, list("SRcov_grid_single_run"),
                envir = environment())
  res_list = parLapply( cl = cluster, x = cov_grid,
                        fun = SRcov_grid_multiple_run,
                        data = data, seed0 = seed0, Nrep = Nrep)
  stopCluster(cluster)

  grid_elapsed <- as.numeric(difftime(
    Sys.time(), grid_start_time, units = "secs"
  ))
  cat(sprintf(
    "[%s] COVERAGE GRID DONE | grid points=%d | elapsed=%.1f min\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    Lgrid, grid_elapsed / 60
  ))
  flush.console()
  
  return(res_list)
}


## ------------------------------------------------------------
## 2. Run - Bounded and Unbounded
## ------------------------------------------------------------
load(paste0("data/Data2019_allRoutes.Rdat"))
data = incidence_matrix
n = nrow(data)

cat(sprintf(
  "[%s] DATA LOADED | observations=%d | species=%d\n",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  nrow(data), ncol(data)
))
flush.console()

eps_grid = seq(0.001,0.1,length.out = 34*2 )
alpha = 0.05
beta = 1e-5
n_max = n

Nrep = 5
seed0 = 4224
num_cores = 34

cat(sprintf(
  paste0(
    "[%s] LAUNCH EPS ANALYSIS | eps range=[%.6f, %.6f] | ",
    "grid points=%d | repliche=%d | cores=%d\n"
  ),
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  min(eps_grid), max(eps_grid), length(eps_grid), Nrep, num_cores
))
flush.console()

#############
# Run
#############
res = SR_grid( eps_grid, data, Nrep, num_cores, seed0, n_max, alpha, beta)
save(res, file = "save/Data2019_allroutes_epsgrid_rev.Rdat")


cat(sprintf(
  "[%s] RESULTS SAVED | file=save/Data2019_allroutes_epsgrid_rev.Rdat\n",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
))
flush.console()

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

cat(sprintf(
  "[%s] COVERAGE ANALYSIS SKIPPED | calls are currently commented out\n",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
))
flush.console()

## ------------------------------------------------------------
## 4. Read and plot
## ------------------------------------------------------------
save_img = FALSE
width = 12; height = 6
cex.lab = 2
cex.axis = 2

stop_here = TRUE
ltype = c(1,1,2,2,2,2)
mycol = c("darkgreen","darkorange","darkgreen","darkorange","deeppink","lightblue")
ygrids = vector("list",length(mycol) )
ygrids[[1]]<-ygrids[[2]]<-ygrids[[3]]<-ygrids[[4]]<-eps_grid
ygrids[[5]]<-ygrids[[6]]<-(1-cov_grid)

if(!stop_here){
  load("save/Data2019_allroutes_epsgrid_rev.Rdat")
  load("save/Data2019_allroutes_covgrid.Rdat")
  
  res_list = lapply(res, function(x) apply(x,2,quantile,probs = 0.5));res_med = do.call(rbind,res_list)
  res_cov_list = lapply(res_cov, function(x) apply(x,2,quantile,probs = 0.5));res_cov_med = do.call(rbind,res_cov_list)
  res_all = cbind(res_med,res_cov_med)
  
  img_name = "img/Data2019_allroutes_SRgrid_rev.pdf"
  if(save_img)
    pdf(img_name, width = width, height = height)
  par(mfrow = c(1,1), bty = "l",  mar = c(3.75,6.5,1,1), mgp=c(5,1,0), las = 1, cex.lab = cex.lab)
  plot(0,0,type = "n", main = "", ylab = "Nstop",
       xaxt = "n", yaxt = "n",
       xlim = range(eps_grid), ylim = c(0,n), 
       xlab = "" )
  axis(2, cex.axis = cex.axis)
  axis(1, cex.axis = cex.axis)
  mtext(paste0(expression(epsilon)," / 1 - coverage"), side = 1, line = 2.5, cex = cex.axis)
  for(i in 1:6){
    points(y = res_all[,i], x = ygrids[[i]], 
           type = "l", lty = ltype[i], 
           lwd = 5, col = mycol[i], pch = 16)
  }
  legend("topright", c("Bounded","Unbounded","Coverage","Chao2009"), 
         col = mycol, lty = ltype, lwd = 5, cex = 1.5)
  
}

script_end_time <- Sys.time()
script_elapsed <- as.numeric(difftime(
  script_end_time, script_start_time, units = "secs"
))
cat(sprintf(
  "[%s] SCRIPT DONE | elapsed=%.1f min\n",
  format(script_end_time, "%Y-%m-%d %H:%M:%S"),
  script_elapsed / 60
))
flush.console()
