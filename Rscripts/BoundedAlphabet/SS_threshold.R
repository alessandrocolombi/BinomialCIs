
setwd("C:/Users/colom/BinomialCIs/Rscripts/BoundedAlphabet")

# Librerie ----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
suppressWarnings(suppressPackageStartupMessages(library(progress)))
suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../../R/Rfunctions.R")

# Sim specific param ----------------------------------------------------------------
RunParallel = FALSE
if(!RunParallel){
  idx = 1 
  param = 1
}

experiments = list(
  "Constant_new" = c(0.001,0.005,0.01,0.015,0.02),
  "Zipfs" = c(0.75,0.8,0.825,0.85,0.9,0.95,1,1.02,1.05),
  "gamma" = c(0.7,0.8,0.9,1,1.1,1.2),
  "Geom" = c(0.25,0.1,0.08,0.05,0.02),
  "Constant" = c(1)
)

# Common parameters -------------------------------------------------------
alfa = 0.05
beta = 1e-5
Brep = 2; Bor = 2; BB = max(Bor,Brep)
seed = 42
set.seed(seed)

save_img = FALSE
  

# Run Constant - new ---------------------------------------------------------------
n = 1000
M = 1500 
soglia = M/n * log(20*M)
soglia


j = 1;
name_exp = names(experiments)[j]
cat("\n Start: ",name_exp,"\n")

# Names
exp_name = paste0("SS_soglia_",name_exp)
file_name = paste0("save/",exp_name,".Rdat")
img_name = paste0("img/",exp_name,".pdf")
cov_name = paste0("save/Coverage/",exp_name,".Rdat")

params = experiments[[j]]
Nparams = length(params)
ub_Bdd_all <- ub_Unbdd_all <- ub_BddOracle_all <- ub_UnbddOracle_all <- ub_bnc_all <- rep(0,Nparams)
avg_S <- rep(0,Nparams)

ii = 1
for(ii in 1:Nparams){
  param = params[ii]
  cat("\n #Exp. = ",ii,"; s = ",param," \n")
  
  ub_bnc_mat <- ub_An_mat <- ub_Ubb_rnorm_mat <- rep(-Inf,Brep) 
  ub_BddOracle_mat <- ub_UbbOracle_mat <- rep(-Inf,Brep) 
  Strue_rep <- rep(0,Brep)
  
  prob_true_mat = matrix(0,nrow = BB, ncol = M) # (BB x M) matrix
  prob_true_mat = t(apply(prob_true_mat, 1, function(x) { sim_Constant_features(M,param) } ))
  data = apply(prob_true_mat, 2, function(pj) {rbinom(n = BB, size = n, prob = pj)} ) # (BB x M) matrix
  # data[b,j]: number of obs. of j-th features in b-th repetition
  
  # Estimators
  b = 1
  for(b in 1:Brep){
    n_i = data[b,] # number of obs. of all features in b-th repetition
    idx_obs = which(n_i > 0) # index of observed features
    Kn = length(idx_obs) # number of obs. features 
    data_obs = n_i[idx_obs] # get observed features
    
    ## Benchmark
    ub_bnc_mat[b] = (log(M/alfa))/n
    
    bn = log(n)
    ## UB Bdd Oracle
    ub_BddOracle_mat[b] = compute_UB_bdd_oracle(n,prob_true_mat[b,],bn,alfa)
    ## UB Analytical
    ub_An_mat[b] = compute_UB_analytical(n,n_i,M,bn,alfa,FALSE)
    
    # Unbounded alphabet:
    S = sum(prob_true_mat[b,])
    Strue_rep[b] = S
    Shat = sum( n_i )/n
    
    ## Oracle
    ub_UbbOracle_mat[b] = compute_UB_Unbdd_oracle(n, alfa, S)
    ## Data dependent
    Sstar = ( sqrt( -log(beta)/(2*n) ) + sqrt( Shat + (-log(beta)/(2*n)) ) )^2
    rn = log( Sstar / (-log(1-alfa+beta)) ) + log(n) - log(log(n))
    ub_Ubb_rnorm_mat[b] = compute_UB_rnorm(n, alfa, beta, rn, Shat)
  }
  
  # Read results
  ub_Bdd_all[ii]   = quantile(ub_An_mat, probs = c(0.5))
  ub_Unbdd_all[ii] = quantile(ub_Ubb_rnorm_mat, probs = c(0.5)) 
  ub_BddOracle_all[ii] = quantile(ub_BddOracle_mat, probs = c(0.5)) 
  ub_UnbddOracle_all[ii] = quantile(ub_UbbOracle_mat, probs = c(0.5)) 
  ub_bnc_all[ii] = quantile(ub_bnc_mat, probs = c(0.5)) 
  avg_S[ii] = quantile(Strue_rep, probs = c(0.5) )
}


res_mat = cbind(avg_S, ub_Bdd_all, ub_Unbdd_all, ub_BddOracle_all, ub_UnbddOracle_all, ub_bnc_all)
res_mat[,2:ncol(res_mat)] = res_mat[,2:ncol(res_mat)] * 1000 
res_mat[,1] = round(res_mat[,1], 2)
res_mat[,2:ncol(res_mat)] = round(res_mat[,2:ncol(res_mat)], 4)
colnames(res_mat) = c("<S>", "Bdd", "Ubd", "BdO", "UbO", "Bnc")
rownames(res_mat) = as.character(experiments[[j]])
res_mat


ymax = max(ub_Bdd_all,ub_Unbdd_all,ub_BddOracle_all,ub_UnbddOracle_all,ub_bnc_all) * 1.05 #18*1e-3
ymin = 7*1e-3
ylabs = round(seq(ymin*1e3,ymax*1e3,by = 1),1)


if(save_img)
  pdf(img_name)
par( mfrow = c(1,1), mar = c(3,3,1,0.5), mgp=c(2,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "S", ylab = "1000 * bound",
     xlim = range(avg_S) , ylim = c(ymin,ymax), 
     main = paste0(" "),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs*1e-3, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
abline(v = soglia, lty = 2, col = "grey45" )
# points( x = avg_S, y = ub_bnc_all, 
#         type = "b", 
#         lwd = 3, pch = 16, lty = 1,
#         col = "black" ) 
points( x = avg_S, y = ub_Bdd_all, 
        type = "b", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkgreen" ) 
points( x = avg_S, y = ub_Unbdd_all, 
        type = "b", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkorange" ) 
# points( x = avg_S, y = ub_BddOracle_all, 
#         type = "b", 
#         lwd = 3, pch = 16, lty = 1,
#         col = "green" ) 
# points( x = avg_S, y = ub_UnbddOracle_all, 
#         type = "b", 
#         lwd = 3, pch = 16, lty = 1,
#         col = "red" ) 
legend("bottomright",c("Bounded","Unbounded"), 
       lwd = 3, col = c("darkgreen","darkorange"))

if(save_img)
  dev.off()


# Run Zipfs ---------------------------------------------------------------
n = 1000
M = 1500 
soglia = M/n * log(20*M)
soglia


j = 2;
name_exp = names(experiments)[j]
cat("\n Start: ",name_exp,"\n")

# Names
exp_name = paste0("SS_soglia_",name_exp)
file_name = paste0("save/",exp_name,".Rdat")
img_name = paste0("img/",exp_name,".pdf")
cov_name = paste0("save/Coverage/",exp_name,".Rdat")

params = experiments[[j]]
Nparams = length(params)
ub_Bdd_all <- ub_Unbdd_all <- ub_BddOracle_all <- ub_UnbddOracle_all <- ub_bnc_all <- rep(0,Nparams)
avg_S <- rep(0,Nparams)

ii = 5
for(ii in 1:Nparams){
  param = params[ii]
  cat("\n #Exp. = ",ii,"; s = ",param," \n")
  
  ub_bnc_mat <- ub_An_mat <- ub_Ubb_rnorm_mat <- rep(-Inf,Brep) 
  ub_BddOracle_mat <- ub_UbbOracle_mat <- rep(-Inf,Brep) 
  Strue_rep <- rep(0,Brep)
  
  prob_true_mat = matrix(0,nrow = BB, ncol = M) # (BB x M) matrix
  prob_true_mat = t(apply(prob_true_mat, 1, function(x) { sim_TruncatedZipfs_features(M,param) } ))
  data = apply(prob_true_mat, 2, function(pj) {rbinom(n = BB, size = n, prob = pj)} ) # (BB x M) matrix
  # data[b,j]: number of obs. of j-th features in b-th repetition
    
  # Estimators
  b = 1
  for(b in 1:Brep){
    n_i = data[b,] # number of obs. of all features in b-th repetition
    idx_obs = which(n_i > 0) # index of observed features
    Kn = length(idx_obs) # number of obs. features 
    data_obs = n_i[idx_obs] # get observed features
      
    ## Benchmark
    ub_bnc_mat[b] = (log(M/alfa))/n
      
    bn = log(n)
    ## UB Bdd Oracle
    ub_BddOracle_mat[b] = compute_UB_bdd_oracle(n,prob_true_mat[b,],bn,alfa)
    ## UB Analytical
    ub_An_mat[b] = compute_UB_analytical(n,n_i,M,bn,alfa,FALSE)
      
    # Unbounded alphabet:
    S = sum(prob_true_mat[b,])
    Strue_rep[b] = S
    Shat = sum( n_i )/n
      
    ## Oracle
    ub_UbbOracle_mat[b] = compute_UB_Unbdd_oracle(n, alfa, S)
    ## Data dependent
    Sstar = ( sqrt( -log(beta)/(2*n) ) + sqrt( Shat + (-log(beta)/(2*n)) ) )^2
    rn = log( Sstar / (-log(1-alfa+beta)) ) + log(n) - log(log(n))
    ub_Ubb_rnorm_mat[b] = compute_UB_rnorm(n, alfa, beta, rn, Shat)
  }
  
  # Read results
  ub_Bdd_all[ii]   = quantile(ub_An_mat, probs = c(0.5))
  ub_Unbdd_all[ii] = quantile(ub_Ubb_rnorm_mat, probs = c(0.5)) 
  ub_BddOracle_all[ii] = quantile(ub_BddOracle_mat, probs = c(0.5)) 
  ub_UnbddOracle_all[ii] = quantile(ub_UbbOracle_mat, probs = c(0.5)) 
  ub_bnc_all[ii] = quantile(ub_bnc_mat, probs = c(0.5)) 
  avg_S[ii] = quantile(Strue_rep, probs = c(0.5) )
}


res_mat = cbind(avg_S, ub_Bdd_all, ub_Unbdd_all, ub_BddOracle_all, ub_UnbddOracle_all, ub_bnc_all)
res_mat[,2:ncol(res_mat)] = res_mat[,2:ncol(res_mat)] * 1000 
res_mat[,1] = round(res_mat[,1], 2)
res_mat[,2:ncol(res_mat)] = round(res_mat[,2:ncol(res_mat)], 4)
colnames(res_mat) = c("<S>", "Bdd", "Ubd", "BdO", "UbO", "Bnc")
rownames(res_mat) = as.character(experiments[[j]])
res_mat


ymax = max(ub_Bdd_all,ub_Unbdd_all,ub_BddOracle_all,ub_UnbddOracle_all,ub_bnc_all) * 1.05 #18*1e-3
ymin = 7*1e-3
ylabs = round(seq(ymin*1e3,ymax*1e3,by = 1),1)


if(save_img)
  pdf(img_name)
par( mfrow = c(1,1), mar = c(3,3,1,0.5), mgp=c(2,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "S", ylab = "1000 * bound",
     xlim = range(avg_S) , ylim = c(ymin,ymax), 
     main = paste0(" "),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs*1e-3, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
abline(v = soglia, lty = 2, col = "grey45" )
# points( x = avg_S, y = ub_bnc_all, 
#         type = "b", 
#         lwd = 3, pch = 16, lty = 1,
#         col = "black" ) 
points( x = avg_S, y = ub_Bdd_all, 
        type = "b", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkgreen" ) 
points( x = avg_S, y = ub_Unbdd_all, 
        type = "b", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkorange" ) 
# points( x = avg_S, y = ub_BddOracle_all, 
#         type = "b", 
#         lwd = 3, pch = 16, lty = 1,
#         col = "green" ) 
# points( x = avg_S, y = ub_UnbddOracle_all, 
#         type = "b", 
#         lwd = 3, pch = 16, lty = 1,
#         col = "red" ) 
legend("bottomright",c("Bounded","Unbounded"), 
       lwd = 3, col = c("darkgreen","darkorange"))

if(save_img)
  dev.off()



# Run Geom ---------------------------------------------------------------
n = 1000
M = 1500 
soglia = M/n * log(20*M)
soglia


j = 4;
name_exp = names(experiments)[j]
cat("\n Start: ",name_exp,"\n")

# Names
exp_name = paste0("SS_soglia_",name_exp)
file_name = paste0("save/",exp_name,".Rdat")
img_name = paste0("img/",exp_name,".pdf")
cov_name = paste0("save/Coverage/",exp_name,".Rdat")

params = experiments[[j]]
Nparams = length(params)
ub_Bdd_all <- ub_Unbdd_all <- ub_BddOracle_all <- ub_UnbddOracle_all <- ub_bnc_all <- rep(0,Nparams)
avg_S <- rep(0,Nparams)

ii = 2
for(ii in 1:Nparams){
  param = params[ii]
  cat("\n #Exp. = ",ii,"; s = ",param," \n")
  
  ub_bnc_mat <- ub_An_mat <- ub_Ubb_rnorm_mat <- rep(-Inf,Brep) 
  ub_BddOracle_mat <- ub_UbbOracle_mat <- rep(-Inf,Brep) 
  Strue_rep <- rep(0,Brep)
  
  prob_true_mat = matrix(0,nrow = BB, ncol = M) # (BB x M) matrix
  prob_true_mat = t(apply(prob_true_mat, 1, function(x) { sim_TruncatedGeom_features(M,param) } ))
  data = apply(prob_true_mat, 2, function(pj) {rbinom(n = BB, size = n, prob = pj)} ) # (BB x M) matrix
  # data[b,j]: number of obs. of j-th features in b-th repetition
  
  # Estimators
  b = 1
  for(b in 1:Brep){
    n_i = data[b,] # number of obs. of all features in b-th repetition
    idx_obs = which(n_i > 0) # index of observed features
    Kn = length(idx_obs) # number of obs. features 
    data_obs = n_i[idx_obs] # get observed features
    
    ## Benchmark
    ub_bnc_mat[b] = (log(M/alfa))/n
    
    bn = log(n)
    ## UB Bdd Oracle
    ub_BddOracle_mat[b] = compute_UB_bdd_oracle(n,prob_true_mat[b,],bn,alfa)
    ## UB Analytical
    ub_An_mat[b] = compute_UB_analytical(n,n_i,M,bn,alfa,FALSE)
    
    # Unbounded alphabet:
    S = sum(prob_true_mat[b,])
    Strue_rep[b] = S
    Shat = sum( n_i )/n
    
    ## Oracle
    ub_UbbOracle_mat[b] = compute_UB_Unbdd_oracle(n, alfa, S)
    ## Data dependent
    Sstar = ( sqrt( -log(beta)/(2*n) ) + sqrt( Shat + (-log(beta)/(2*n)) ) )^2
    rn = log( Sstar / (-log(1-alfa+beta)) ) + log(n) - log(log(n))
    ub_Ubb_rnorm_mat[b] = compute_UB_rnorm(n, alfa, beta, rn, Shat)
  }
  
  # Read results
  ub_Bdd_all[ii]   = quantile(ub_An_mat, probs = c(0.5))
  ub_Unbdd_all[ii] = quantile(ub_Ubb_rnorm_mat, probs = c(0.5)) 
  ub_BddOracle_all[ii] = quantile(ub_BddOracle_mat, probs = c(0.5)) 
  ub_UnbddOracle_all[ii] = quantile(ub_UbbOracle_mat, probs = c(0.5)) 
  ub_bnc_all[ii] = quantile(ub_bnc_mat, probs = c(0.5)) 
  avg_S[ii] = quantile(Strue_rep, probs = c(0.5) )
}


res_mat = cbind(avg_S, ub_Bdd_all, ub_Unbdd_all, ub_BddOracle_all, ub_UnbddOracle_all, ub_bnc_all)
res_mat[,2:ncol(res_mat)] = res_mat[,2:ncol(res_mat)] * 1000 
res_mat[,1] = round(res_mat[,1], 2)
res_mat[,2:ncol(res_mat)] = round(res_mat[,2:ncol(res_mat)], 4)
colnames(res_mat) = c("<S>", "Bdd", "Ubd", "BdO", "UbO", "Bnc")
rownames(res_mat) = as.character(experiments[[j]])
res_mat


ymax = max(ub_Bdd_all,ub_Unbdd_all,ub_BddOracle_all,ub_UnbddOracle_all,ub_bnc_all) * 1.05 #18*1e-3
ymin = 7*1e-3
ylabs = round(seq(ymin*1e3,ymax*1e3,by = 1),1)


if(save_img)
  pdf(img_name)
par( mfrow = c(1,1), mar = c(3,3,1,0.5), mgp=c(2,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "S", ylab = "1000 * bound",
     xlim = range(avg_S) , ylim = c(ymin,ymax), 
     main = paste0(" "),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs*1e-3, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
abline(v = soglia, lty = 2, col = "grey45" )
# points( x = avg_S, y = ub_bnc_all, 
#         type = "b", 
#         lwd = 3, pch = 16, lty = 1,
#         col = "black" ) 
points( x = avg_S, y = ub_Bdd_all, 
        type = "b", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkgreen" ) 
points( x = avg_S, y = ub_Unbdd_all, 
        type = "b", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkorange" ) 
# points( x = avg_S, y = ub_BddOracle_all, 
#         type = "b", 
#         lwd = 3, pch = 16, lty = 1,
#         col = "green" ) 
# points( x = avg_S, y = ub_UnbddOracle_all, 
#         type = "b", 
#         lwd = 3, pch = 16, lty = 1,
#         col = "red" ) 
legend("bottomright",c("Bounded","Unbounded"), 
       lwd = 3, col = c("darkgreen","darkorange"))

if(save_img)
  dev.off()








# Run gamma ---------------------------------------------------------------
n = 1000
M = 1500 
soglia = M/n * log(20*M)
soglia


j = 3;
name_exp = names(experiments)[j]
cat("\n Start: ",name_exp,"\n")

# Names
exp_name = paste0("SS_soglia_",name_exp)
file_name = paste0("save/",exp_name,".Rdat")
img_name = paste0("img/",exp_name,".pdf")
cov_name = paste0("save/Coverage/",exp_name,".Rdat")

params = experiments[[j]]
Nparams = length(params)
ub_Bdd_all <- ub_Unbdd_all <- ub_BddOracle_all <- ub_UnbddOracle_all <- ub_bnc_all <- rep(0,Nparams)
avg_S <- rep(0,Nparams)

ii = 2
for(ii in 1:Nparams){
  param = params[ii]
  cat("\n #Exp. = ",ii,"; s = ",param," \n")
  
  ub_bnc_mat <- ub_An_mat <- ub_Ubb_rnorm_mat <- rep(-Inf,Brep) 
  ub_BddOracle_mat <- ub_UbbOracle_mat <- rep(-Inf,Brep) 
  Strue_rep <- rep(0,Brep)
  
  prob_true_mat = matrix(0,nrow = BB, ncol = M) # (BB x M) matrix
  prob_true_mat = t(apply(prob_true_mat, 1, function(x) { sim_TruncatedZipfs_features(M,param) } ))
  data = apply(prob_true_mat, 2, function(pj) {rbinom(n = BB, size = n, prob = pj)} ) # (BB x M) matrix
  # data[b,j]: number of obs. of j-th features in b-th repetition
  
  # Estimators
  b = 1
  for(b in 1:Brep){
    n_i = data[b,] # number of obs. of all features in b-th repetition
    idx_obs = which(n_i > 0) # index of observed features
    Kn = length(idx_obs) # number of obs. features 
    data_obs = n_i[idx_obs] # get observed features
    
    ## Benchmark
    ub_bnc_mat[b] = (log(M/alfa))/n
    
    bn = log(n)
    ## UB Bdd Oracle
    ub_BddOracle_mat[b] = compute_UB_bdd_oracle(n,prob_true_mat[b,],bn,alfa)
    ## UB Analytical
    ub_An_mat[b] = compute_UB_analytical(n,n_i,M,bn,alfa,FALSE)
    
    # Unbounded alphabet:
    S = sum(prob_true_mat[b,])
    Strue_rep[b] = S
    Shat = sum( n_i )/n
    
    ## Oracle
    ub_UbbOracle_mat[b] = compute_UB_Unbdd_oracle(n, alfa, S)
    ## Data dependent
    Sstar = ( sqrt( -log(beta)/(2*n) ) + sqrt( Shat + (-log(beta)/(2*n)) ) )^2
    rn = log( Sstar / (-log(1-alfa+beta)) ) + log(n) - log(log(n))
    ub_Ubb_rnorm_mat[b] = compute_UB_rnorm(n, alfa, beta, rn, Shat)
  }
  
  # Read results
  ub_Bdd_all[ii]   = quantile(ub_An_mat, probs = c(0.5))
  ub_Unbdd_all[ii] = quantile(ub_Ubb_rnorm_mat, probs = c(0.5)) 
  ub_BddOracle_all[ii] = quantile(ub_BddOracle_mat, probs = c(0.5)) 
  ub_UnbddOracle_all[ii] = quantile(ub_UbbOracle_mat, probs = c(0.5)) 
  ub_bnc_all[ii] = quantile(ub_bnc_mat, probs = c(0.5)) 
  avg_S[ii] = quantile(Strue_rep, probs = c(0.5) )
}


res_mat = cbind(avg_S, ub_Bdd_all, ub_Unbdd_all, ub_BddOracle_all, ub_UnbddOracle_all, ub_bnc_all)
res_mat[,2:ncol(res_mat)] = res_mat[,2:ncol(res_mat)] * 1000 
res_mat[,1] = round(res_mat[,1], 2)
res_mat[,2:ncol(res_mat)] = round(res_mat[,2:ncol(res_mat)], 4)
colnames(res_mat) = c("<S>", "Bdd", "Ubd", "BdO", "UbO", "Bnc")
rownames(res_mat) = as.character(experiments[[j]])
res_mat


ymax = max(ub_Bdd_all,ub_Unbdd_all,ub_BddOracle_all,ub_UnbddOracle_all,ub_bnc_all) * 1.05 #18*1e-3
ymin = 7*1e-3
ylabs = round(seq(ymin*1e3,ymax*1e3,by = 1),1)


if(save_img)
  pdf(img_name)
par( mfrow = c(1,1), mar = c(3,3,1,0.5), mgp=c(2,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "S", ylab = "1000 * bound",
     xlim = range(avg_S) , ylim = c(ymin,ymax), 
     main = paste0(" "),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs*1e-3, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
abline(v = soglia, lty = 2, col = "grey45" )
# points( x = avg_S, y = ub_bnc_all, 
#         type = "b", 
#         lwd = 3, pch = 16, lty = 1,
#         col = "black" ) 
points( x = avg_S, y = ub_Bdd_all, 
        type = "b", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkgreen" ) 
points( x = avg_S, y = ub_Unbdd_all, 
        type = "b", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkorange" ) 
# points( x = avg_S, y = ub_BddOracle_all, 
#         type = "b", 
#         lwd = 3, pch = 16, lty = 1,
#         col = "green" ) 
# points( x = avg_S, y = ub_UnbddOracle_all, 
#         type = "b", 
#         lwd = 3, pch = 16, lty = 1,
#         col = "red" ) 
legend("bottomright",c("Bounded","Unbounded"), 
       lwd = 3, col = c("darkgreen","darkorange"))

if(save_img)
  dev.off()








# Run - Constant ---------------------------------------------------------------------
n = 1000
M = 5000 # this is Mtrue, the true size of the alphabet
Smin = 1
Smax = 200
Sgrid = seq(Smin,Smax,length.out = 10)
Nexp = length(Sgrid)


j = 5;
name_exp = names(experiments)[j]
cat("\n Start: ",name_exp,"\n")

# Names
exp_name = paste0("SS_soglia_",name_exp)
file_name = paste0("save/",exp_name,".Rdat")
img_name = paste0("img/",exp_name,".pdf")
cov_name = paste0("save/Coverage/",exp_name,".Rdat")

params = experiments[[j]]
Nparams = length(params)
ub_Bdd_all <- ub_Unbdd_all <- ub_BddOracle_all <- ub_UnbddOracle_all <- ub_bnc_all <- matrix(0, nrow = Nparams, ncol = Nexp)
for(ii in 1:Nparams){
  cat("\n #Exp. = ",ii," \n")
  param = params[ii]
  
  ub_bnc_mat <- ub_An_mat <- ub_Ubb_rnorm_mat <- matrix(-Inf,Nexp,Brep) 
  ub_BddOracle_mat <- ub_UbbOracle_mat <- matrix(-Inf,Nexp,Brep) 
  
  for(ss in 1:Nexp){
    S = Sgrid[ss]
    cat("\n #S = ",S," \n")
    BB = max(Bor,Brep)
    prob_true_mat = matrix(0,nrow = BB, ncol = M) # (BB x M) matrix
    prob_true_mat = t(apply(prob_true_mat, 1, function(x) { sim_CntS_Constant(S,M) } ))
    data = apply(prob_true_mat, 2, function(pj) {rbinom(n = BB, size = n, prob = pj)} ) # (BB x M) matrix
    # data[b,j]: number of obs. of j-th features in b-th repetition
    
    # Estimators
    b = 1
    for(b in 1:Brep){
      n_i = data[b,] # number of obs. of all features in b-th repetition
      idx_obs = which(n_i > 0) # index of observed features
      Kn = length(idx_obs) # number of obs. features 
      data_obs = n_i[idx_obs] # get observed features
      
      ## Benchmark
      ub_bnc_mat[ss,b] = (log(M/alfa))/n
      
      bn = log(n)
      ## UB Bdd Oracle
      ub_BddOracle_mat[ss,b] = compute_UB_bdd_oracle(n,prob_true_mat[b,],bn,alfa)
      ## UB Analytical
      ub_An_mat[ss,b] = compute_UB_analytical(n,n_i,M,bn,alfa,FALSE)
      
      # Unbounded alphabet:
      Shat = sum( n_i )/n
      
      ## Oracle
      ub_UbbOracle_mat[ss,b] = compute_UB_Unbdd_oracle(n, alfa, S)
      ## Data dependent
      Sstar = ( sqrt( -log(beta)/(2*n) ) + sqrt( Shat + (-log(beta)/(2*n)) ) )^2
      rn = log( Sstar / (-log(1-alfa+beta)) ) + log(n) - log(log(n))
      ub_Ubb_rnorm_mat[ss,b] = compute_UB_rnorm(n, alfa, beta, rn, Shat)
    }
  }
  
  # Read results
  ub_Bdd_all[ii,]   = apply(ub_An_mat, 1, quantile, probs = c(0.5))
  ub_Unbdd_all[ii,] = apply(ub_Ubb_rnorm_mat, 1, quantile, probs = c(0.5))
  ub_BddOracle_all[ii,] = apply(ub_BddOracle_mat, 1, quantile, probs = c(0.5))
  ub_UnbddOracle_all[ii,] = apply(ub_UbbOracle_mat, 1, quantile, probs = c(0.5))
  ub_bnc_all[ii,] = apply(ub_bnc_mat, 1, quantile, probs = c(0.5))
  
}

ymax = max(ub_Bdd_all,ub_Unbdd_all,ub_BddOracle_all,ub_UnbddOracle_all,ub_bnc_all) * 1.05 #18*1e-3
ymin = 7*1e-3
ylabs = round(seq(ymin*1e3,ymax*1e3,by = 1),1)


if(save_img)
  pdf(img_name)
par( mfrow = c(1,1), mar = c(3,3,1,0.5), mgp=c(2,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "S", ylab = "1000 * bound",
     xlim = range(Sgrid) , ylim = c(ymin,ymax), 
     main = paste0(" "),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs*1e-3, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
abline(v = (M/n)*log(20*M), lty = 2, col = "grey45" )
for(ii in 1:Nparams){
  # points( x = Sgrid, y = ub_bnc_all[ii,], 
  #         type = "b", 
  #         lwd = 3, pch = 16, lty = 1,
  #         col = "black" ) 
  points( x = Sgrid, y = ub_Bdd_all[ii,], 
          type = "b", 
          lwd = 3, pch = 16, lty = 1,
          col = "darkgreen" ) 
  points( x = Sgrid, y = ub_Unbdd_all[ii,], 
          type = "b", 
          lwd = 3, pch = 16, lty = 1,
          col = "darkorange" ) 
  # points( x = Sgrid, y = ub_BddOracle_all[ii,], 
  #         type = "b", 
  #         lwd = 3, pch = 16, lty = 1,
  #         col = "green" ) 
  # points( x = Sgrid, y = ub_UnbddOracle_all[ii,], 
  #         type = "b", 
  #         lwd = 3, pch = 16, lty = 1,
  #         col = "red" ) 
  legend("bottomright",c("Bounded","Unbounded"), 
         lwd = 3, col = c("darkgreen","darkorange"))
  if(save_img)
    dev.off()
}

