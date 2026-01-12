
setwd("C:/Users/colom/BinomialCIs/Rscripts/SimulationStudy_firstSub")

# Librerie ----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
suppressWarnings(suppressPackageStartupMessages(library(progress)))
suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../../R/Rfunctions.R")

# Sim specific param ----------------------------------------------------------------
RunParallel = TRUE
if(!RunParallel){
  idx = 1 
  param = 1
}

names = c("Zipfs","Geom","Constant")#,"Uniform","TrUnif")
name = names[3] # choose the name here


# Common parameters -------------------------------------------------------
alfa = 0.05
beta = 1e-5
Brep = 100; Bor = 100
seed = 42
set.seed(seed)
trimmed_param = get_first3digits(param,4)

# n fix, M varies ---------------------------------------------------------
save_cov = TRUE
save_exp = TRUE
save_img = TRUE

n = 2000
Mgrid_min =  2000 
Mgrid_max =  10000
Mgrid_step = 250

Mgrid = seq(Mgrid_min,Mgrid_max,by = Mgrid_step)
Nexp = length(Mgrid)


exp_name = paste0("Exp1_nfix_",name,"_",trimmed_param)
file_name = paste0("save/",exp_name,".Rdat")
img_name = paste0("img/",name,"/",exp_name,".pdf")
cov_name = paste0("save/Coverage/",name,"/",exp_name,".Rdat")

## Run  --------------------------------------------------------------------


run_n_fix = TRUE
if(run_n_fix){
  
  # Upper/Lower bound
  oracle_mat <- matrix(-Inf,Nexp,Bor) 
  ub_Bnf_mat <- ub_Bdd_mat <- ub_Ubd_mat <- matrix(-Inf,Nexp,Brep) 
  oracle <- rep(NA,Nexp)
  
  # Coverage
  cov_Bnf_mat <- cov_Bdd_mat <- cov_Ubd_mat <- matrix(0,Nexp,Brep) 
  
  
  pb <- progress_bar$new(total = Nexp); ii = 1
  for(ii in 1:Nexp){
    cat("\n #Exp. = ",ii," \n")
    M = Mgrid[ii]
    BB = max(Bor,Brep)
    prob_true_mat = matrix(0,nrow = BB, ncol = M) # (BB x M) matrix
    prob_true_mat = t(apply(prob_true_mat, 1, function(x) { sim_generic(name,M,param) } ))
    data = apply(prob_true_mat, 2, function(pj) {rbinom(n = BB, size = n, prob = pj)} ) # (BB x M) matrix
    # data[b,j]: number of obs. of j-th features in b-th repetition
    
    # Oracle
    Mmax = rep(NA,Bor) # needed to compute the oracle
    b = 1
    for(b in 1:Bor){
      n_i = data[b,] # number of obs. of all features in b-th repetition
      idx_obs = which(n_i > 0) # index of observed features
      Kn = length(idx_obs) # number of obs. features 
      data_obs = n_i[idx_obs] # get observed features
      if(Kn == M){
        # if here, the whole alphabet has been observed
        Mmax[b] = 0
      }else{
        # if here, some symbols have not been observed
        idx_unobs = which(n_i == 0)
        Mmax[b] = max(prob_true_mat[b,idx_unobs])
      }
      
      ## Get pmax = max{pj : Nj = 0} in b-th repetition
      oracle_mat[ii,b] = Mmax[b]  
      
    }
    
    ## Oracle 
    oracle[ii] = quantile(Mmax, probs = 1-alfa)
    
    # Estimators
    b = 1
    for(b in 1:Brep){
      n_i = data[b,] # number of obs. of all features in b-th repetition
      idx_obs = which(n_i > 0) # index of observed features
      Kn = length(idx_obs) # number of obs. features 
      data_obs = n_i[idx_obs] # get observed features
      pmax = -Inf 
      if(any(n_i == 0)){
        pmax = max( prob_true_mat[b,n_i == 0] ) # get pmax true
      }
      
      ## Benchmark
      ub_Bnf_mat[ii,b] = log(M/alfa)/n
      if(pmax <= ub_Bnf_mat[ii,b]){
        cov_Bnf_mat[ii,b] = 1
      }
      
      ## Bounded
      bn = log(n)
      ub_Bdd_mat[ii,b] = compute_UB_analytical(n,n_i,M,bn,alfa,FALSE)
      if(pmax <= ub_Bdd_mat[ii,b]){
        cov_Bdd_mat[ii,b] = 1
      }
      
      # Unbounded alphabet:
      Shat = sum( n_i )/n
      Sstar = ( sqrt( -log(beta)/(2*n) ) + sqrt( Shat + (-log(beta)/(2*n)) ) )^2
      rn = log( Sstar / (-log(1-alfa+beta)) ) + log(n) - log(log(n))
      ub_Ubd_mat[ii,b] = compute_UB_rnorm(n, alfa, beta, rn, Shat)
      if(pmax <= ub_Ubd_mat[ii,b]){
        cov_Ubd_mat[ii,b] = 1
      }
    }
    
    pb$tick()
  }
  
  
  save_res = list( "ub_Bnf_mat" = ub_Bnf_mat,
                   "oracle_mat"   = oracle_mat,
                   "ub_Bdd_mat" = ub_Bdd_mat,
                   "ub_Ubd_mat" = ub_Ubd_mat,
                   "oracle" = oracle )
  
  if(save_exp)
    save(save_res, file = file_name)
  
}


## Sample coverage ---------------------------------------------------------
cov_mat = Mgrid
cov_mat = cbind(cov_mat, colMeans( t(cov_Bnf_mat) ))
cov_mat = cbind(cov_mat, colMeans( t(cov_Bdd_mat)))
cov_mat = cbind(cov_mat, colMeans( t(cov_Ubd_mat)))
colnames(cov_mat) = c("M","Bonferroni","Bounded","Unbounded")

if(save_cov)
  save(cov_mat, file = cov_name)

## Final summary and plot ---------------------------------------------------

ub_Bnf = apply(ub_Bnf_mat, 1, quantile, probs = c(0.025,0.5,0.975))
ub_Bdd = apply(ub_Bdd_mat, 1, quantile,    probs = c(0.025,0.5,0.975))
ub_Ubd = apply(ub_Ubd_mat, 1, quantile,    probs = c(0.025,0.5,0.975))

ymax = max(ub_Bnf[3,],ub_Bdd[3,],ub_Ubd[3,],oracle) * 1.05 #18*1e-3
ymin = min(ub_Bnf[1,],ub_Bdd[1,],ub_Ubd[1,],oracle) #5*1e-3
if(ymin < 0)
  ymin = 0
ylabs = round(seq(ymin*1e3,ymax*1e3,by = 1),1)


if(save_img)
  pdf(img_name)
par( mfrow = c(1,1), mar = c(3,3,1,0.5), mgp=c(2,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "M", ylab = "1000 * bound",
     xlim = c(min(Mgrid),max(Mgrid) ) , ylim = c(ymin,ymax), 
     main = paste0(" "),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs*1e-3, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
points( x = Mgrid, y = ub_Bnf[2,], 
        type = "l", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkred" ) 
points( x = Mgrid, y = ub_Bdd[2,], 
        type = "l", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkgreen" ) 
points( x = Mgrid, y = ub_Ubd[2,],
        type = "l",
        lwd = 3, pch = 16, lty = 1,
        col = "darkorange" )
legend("bottomright",c("Bonferroni","Bounded","Unbounded"), 
       lwd = 3, col = c("darkred","darkgreen","darkorange"))
if(save_img)
  dev.off()


# M fix, n varies ---------------------------------------------------------
save_cov = TRUE
save_exp = TRUE
save_img = TRUE

M = 5000
Ngrid_max =  5000
Ngrid_min =  1000
Ngrid_step = 250

Ngrid = seq(Ngrid_min,Ngrid_max,by = Ngrid_step)
Nexp = length(Ngrid)

exp_name = paste0("Exp1_Mfix_",name,"_",trimmed_param)
file_name = paste0("save/",exp_name,".Rdat")
img_name = paste0("img/",name,"/",exp_name,".pdf")
cov_name = paste0("save/Coverage/",name,"/",exp_name,".Rdat")

## Run  --------------------------------------------------------------------


run_M_fix = TRUE
if(run_M_fix){
  
  oracle_mat <- matrix(-Inf,Nexp,Bor) 
  ub_Bnf_mat <- ub_Bdd_mat <- ub_Ubd_mat <- matrix(-Inf,Nexp,Brep) 
  oracle <- rep(NA,Nexp)
  
  # Coverage
  cov_Bnf_mat <- cov_Bdd_mat <- cov_Ubd_mat <- matrix(0,Nexp,Brep) 
  
  
  # Generate true distribution
  BB = max(Bor,Brep)
  prob_true_mat = matrix(0,nrow = BB, ncol = M) # (BB x M) matrix
  prob_true_mat = t(apply(prob_true_mat, 1, function(x) { sim_generic(name,M,param) } ))
  
  pb <- progress_bar$new(total = Nexp)
  for(ii in 1:Nexp){
    cat("\n #Exp. = ",ii," \n")
    n = Ngrid[ii]
    data = apply(prob_true_mat, 2, function(pj) {rbinom(n = BB, size = n, prob = pj)} ) # (BB x M) matrix
    # data[b,j]: number of obs. of j-th features in b-th repetition
    
    # Oracle
    Mmax = rep(NA,Bor) # needed to compute the oracle
    for(b in 1:Bor){
      n_i = data[b,] # number of obs. of all features in b-th repetition
      idx_obs = which(n_i > 0) # index of observed features
      Kn = length(idx_obs) # number of obs. features 
      data_obs = n_i[idx_obs] # get observed features
      if(Kn == M){
        # if here, the whole alphabet has been observed
        Mmax[b] = 0
      }else{
        # if here, some symbols have not been observed
        idx_unobs = which(n_i == 0)
        Mmax[b] = max(prob_true_mat[b,idx_unobs])
      }
      
      ## Get pmax = max{pj : Nj = 0} in b-th repetition
      oracle_mat[ii,b] = Mmax[b]  
      
    }
    
    ## Oracle 
    oracle[ii] = quantile(Mmax, probs = 1-alfa)
    
    # Estimators
    for(b in 1:Brep){
      n_i = data[b,] # number of obs. of all features in b-th repetition
      idx_obs = which(n_i > 0) # index of observed features
      Kn = length(idx_obs) # number of obs. features 
      data_obs = n_i[idx_obs] # get observed features
      pmax = -Inf    
      if(any(n_i == 0)){
        pmax = max( prob_true_mat[b,n_i == 0] ) # get pmax true
      }
      
      ## Benchmark
      ub_Bnf_mat[ii,b] = log(M/alfa)/n
      if(pmax <= ub_Bnf_mat[ii,b]){
        cov_Bnf_mat[ii,b] = 1
      }
      
      
      ## Bounded
      bn = log(n)
      ub_Bdd_mat[ii,b] = compute_UB_analytical(n,n_i,M,bn,alfa,FALSE)
      if(pmax <= ub_Bdd_mat[ii,b]){
        cov_Bdd_mat[ii,b] = 1
      }
      
      ## Unbounded alphabet:
      Shat = sum( n_i )/n
      Sstar = ( sqrt( -log(beta)/(2*n) ) + sqrt( Shat + (-log(beta)/(2*n)) ) )^2
      rn = log( Sstar / (-log(1-alfa+beta)) ) + log(n) - log(log(n))
      ub_Ubd_mat[ii,b] = compute_UB_rnorm(n, alfa, beta, rn, Shat)
      if(pmax <= ub_Ubd_mat[ii,b]){
        cov_Ubd_mat[ii,b] = 1
      }
    }
    
    pb$tick()
  }
  
  
  save_res = list( "ub_Bnf_mat" = ub_Bnf_mat,
                   "oracle_mat"   = oracle_mat,
                   "ub_Bdd_mat" = ub_Bdd_mat,
                   "ub_Ubd_mat" = ub_Ubd_mat,
                   "oracle" = oracle )
  
  if(save_exp)
    save(save_res, file = file_name)
  
}

## Sample coverage ---------------------------------------------------------
cov_mat = Ngrid
cov_mat = cbind(cov_mat, colMeans( t(cov_Bnf_mat) ))
cov_mat = cbind(cov_mat, colMeans( t(cov_Bdd_mat)))
cov_mat = cbind(cov_mat, colMeans( t(cov_Ubd_mat)))
colnames(cov_mat) = c("M","Bonferroni","Bounded","Unbounded")

if(save_cov)
  save(cov_mat, file = cov_name)


## Final summary and plot ---------------------------------------------------

ub_Bnf = apply(ub_Bnf_mat, 1, quantile, probs = c(0.025,0.5,0.975))
ub_Bdd = apply(ub_Bdd_mat, 1, quantile, probs = c(0.025,0.5,0.975))
ub_Ubd = apply(ub_Ubd_mat, 1, quantile, probs = c(0.025,0.5,0.975))

ymax = max(ub_Bnf[3,],ub_Bdd[3,],ub_Ubd[3,],oracle) * 1.05 #18*1e-3
ymin = min(ub_Bnf[1,],ub_Bdd[1,],ub_Ubd[1,],oracle) #5*1e-3
if(ymin < 0)
  ymin = 0
ylabs = round(seq(ymin*1e3,ymax*1e3,by = 1),1)


if(save_img)
  pdf(img_name)
par( mfrow = c(1,1), mar = c(4,4,1,0.5), mgp=c(2.5,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "n", ylab = "1000*bound",
     xlim = c(min(Ngrid),max(Ngrid) ) , ylim = c(ymin,ymax), 
     main = paste0(" "),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs*1e-3, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
points( x = Ngrid, y = ub_Bnf[2,], 
        type = "l", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkred" ) 
points( x = Ngrid, y = ub_Bdd[2,], 
        type = "l", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkgreen" ) 
points( x = Ngrid, y = ub_Ubd[2,],
        type = "l",
        lwd = 3, pch = 16, lty = 1,
        col = "darkorange" )
legend("topright",c("Bonferroni","Bounded","Unbounded"), 
       lwd = 3, col = c("darkred","darkgreen","darkorange"))
if(save_img)
  dev.off()