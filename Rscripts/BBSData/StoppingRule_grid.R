wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100)
choose_wd = wd_vec[3] # <--- modify here
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
## ------------------------------------------------------------

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


## ------------------------------------------------------------
## 2. Run
## ------------------------------------------------------------
load(paste0("data/Data2019_allRoutes.Rdat"))
data = incidence_matrix
n = nrow(data)

eps_grid = seq(0.001,0.05,length.out = 3)
alpha = 0.05
beta = 1e-5

Nrep = 30
seed0 = 4224
num_cores = 30

res = SR_grid( eps_grid, data, Nrep, num_cores, seed0, alpha, beta)
save(res, file = "save/Data2019_allroutes_epsgrid.Rdat")

## ------------------------------------------------------------
## 3. Read and plot
## ------------------------------------------------------------
stop_here = TRUE
if(!stop_here){
  load("save/Data2019_allroutes_epsgrid.Rdat")
  colMeans(res)
  apply(res, 2, quantile, probs = c(0.25,0.5,0.75))
  
  
  mycol = c("darkgreen","darkred","deeppink","lightblue")
  
  par(mfrow = c(1,1),bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,1,0))
  plot(0,0,type = "n", main = "eps = 0.05", ylab = "Nstop",
       xlim = c(0.5,4.5), ylim = c(150,3200), xlab = "")
  for(i in 1:4){
    boxplot(res[,i], at = i, add = T, 
            col = mycol[i], pch = 16, yaxt = "n")
  }
  
}
