
setwd("C:/Users/colom/BinomialCIs/Rscripts/BoundedAlphabet")

# Librerie ----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
suppressWarnings(suppressPackageStartupMessages(library(progress)))
suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../../R/Rfunctions.R")


# Functions ---------------------------------------------------------------

# Build p_j ∝ j^{-gamma} on {1, ..., M}
build_zipf_p <- function(M = 1500L, gamma = 0.75) {
  j <- seq_len(M) + 1
  p <- j^(-gamma)
  p
}

build_uniform_p <- function(M = 1500L, c=0.05) {
  p <- rep(c, M)
  p
}

build_geom_p <- function(M = 1500L, a = 0.1) {
  if (a >= 1)
    stop("a must be strictly less than 1")
  j <- 2:(M + 1L)
  # This matches your sim_TruncatedGeom_features:
  # w_j = (1 - a)^(j-1), j = 2,...,M+1  -> (1-a)^1,...,(1-a)^M
  w <- (1 - a)^(j - 1L)
  as.numeric(w)
}


build_generic_p = function(name,M,param){
  if(name == "Zipfs"){
    build_zipf_p(M,param)
  }else if( name == "Geom" ){
    build_geom_p(M,param)
  }else if( name == "Uniform"){
    build_uniform_p(M, param)
  }
  else
    stop("Invalid name")
}



# Options -----------------------------------------------------------------
alfa = 0.05
beta = 1e-5

experiments = list(
  "Zipfs" = c(1.05),
  "Geom"  = c(0.05),
  "Uniform" = c(0.006,0.05)
)

M = 1500
n = 2500
eps = 0.005
ngrid = ceiling(seq(10,n, length.out = 100))
Nexp = length(ngrid)
BB <- Brep <- 5


# Run ---------------------------------------------------------------------


save_plot_all = FALSE

par(mfrow = c(2,2),bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,1,0))
i = 1
for(i in 1:length(experiments)){
  name_exp = names(experiments)[i]
  cat("\n Start: ",name_exp,"\n")
  for(p in  experiments[[i]]){
    # trimmed_param = get_first3digits(p,4)
    cat("\n p = ",p,"\n")
    
    ub_Bd_mat <- ub_Unb_mat <- matrix(-Inf,Nexp,Brep) 
    
    prob_true_mat = matrix(0,nrow = BB, ncol = M) # (BB x M) matrix
    prob_true_mat = t(apply(prob_true_mat, 1, function(x) { build_generic_p(name_exp,M,p) } ))
    for(jj in seq_along(ngrid)){
      n = ngrid[jj]
      data = apply(prob_true_mat, 2, function(pj) {rbinom(n = BB, size = n, prob = pj)} ) # (BB x M) matrix
      # data[b,j]: number of obs. of j-th features in b-th repetition
      for(b in 1:BB){
        n_i = data[b,] # number of obs. of all features in b-th repetition
        idx_obs = which(n_i > 0) # index of observed features
        Kn = length(idx_obs) # number of obs. features 
        data_obs = n_i[idx_obs] # get observed features
        
        ## Bounded alphabet
        bn = log(n)
        ub_Bd_mat[jj,b] = compute_UB_analytical(n,n_i,M,bn,alfa,FALSE)
        
        # Unbounded alphabet:
        Shat = sum( n_i )/n
        Sstar = ( sqrt( -log(beta)/(2*n) ) + sqrt( Shat + (-log(beta)/(2*n)) ) )^2
        rn = log( Sstar / (-log(1-alfa+beta)) ) + log(n) - log(log(n))
        ub_Unb_mat[jj,b] = compute_UB_rnorm(n, alfa, beta, rn, Shat)
      }
      
    }
    
    ub_Bd  = apply(ub_Bd_mat, 1, quantile, probs = c(0.025,0.5,0.975))
    ub_Unb = apply(ub_Unb_mat, 1, quantile,    probs = c(0.025,0.5,0.975))
    
    ymax = 0.1; ymin = 0
    ylabs = round(seq(ymin,ymax,by = 0.05),2)
    main_plt = paste0(name_exp,"- ",p)
    plot(0,0,  yaxt = "n",
         xlab = "n", ylab = "bound",
         xlim = c(min(ngrid),max(ngrid) ) , ylim = c(0,0.2), 
         main = main_plt,
         type = "n")
    grid(lty = 1,lwd = 1, col = "gray90" )
    axis(side = 2, at = ylabs, 
         labels = ylabs, las = 1, 
         cex.axis = 1 )
    points( x = ngrid, y = ub_Bd[2,], 
            type = "l", 
            lwd = 3, pch = 16, lty = 1,
            col = "darkgreen" ) 
    points( x = ngrid, y = ub_Unb[2,], 
            type = "l", 
            lwd = 3, pch = 16, lty = 1,
            col = "darkorange" ) 
    abline(h = eps, lty = 2, col = "red")
    legend("topright",c("Bounded","Unbounded"), 
           lwd = 3, col = c("darkgreen","darkorange"))
    

  }
}



