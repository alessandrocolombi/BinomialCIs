
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

names = c("Zipfs","Constant")
name = names[1] # choose the name here

# Common parameters -------------------------------------------------------
alfa = 0.05
beta = 1e-5
Brep = 100; Bor = 100
seed = 42
set.seed(seed)

save_cov = TRUE
save_img = TRUE
save_exp = TRUE

n = 1000
M = 5000 # this is Mtrue, the true size of the alphabet
Madd_grid = c(1e0,1e1,1e2,1e3,1e4,1e5,1e6)
logMadd_grid = log10(Madd_grid)
Mguess_grid = M + Madd_grid

Nexp = length(Madd_grid)

trimmed_param = get_first3digits(param,4)
exp_name = paste0("Exp2_",name,"_",trimmed_param)
file_name = paste0("save/",exp_name,".Rdat")
img_name = paste0("img/",name,"/",exp_name,".pdf")
cov_name = paste0("save/Coverage/",name,"/",exp_name,".Rdat")

## Run  --------------------------------------------------------------------


run_n_fix = TRUE
if(run_n_fix){
  
  ub_Bdd_mat <- ub_Ubd_mat <- matrix(-Inf,Nexp,Brep) 
  
  pb <- progress_bar$new(total = Nexp); ii = 1
  for(ii in 1:Nexp){
    cat("\n #Exp. = ",ii," \n")
    BB = max(Bor,Brep)
    prob_true_mat = matrix(0,nrow = BB, ncol = M) # (BB x M) matrix
    prob_true_mat = t(apply(prob_true_mat, 1, function(x) { sim_generic(name,M,param) } ))
    data = apply(prob_true_mat, 2, function(pj) {rbinom(n = BB, size = n, prob = pj)} ) # (BB x M) matrix
    # data[b,j]: number of obs. of j-th features in b-th repetition
    
    # Estimators
    b = 1
    for(b in 1:Brep){
      n_i = data[b,] # number of obs. of all features in b-th repetition
      idx_obs = which(n_i > 0) # index of observed features
      Kn = length(idx_obs) # number of obs. features 
      data_obs = n_i[idx_obs] # get observed features
      
      Mguess = Mguess_grid[ii]
      ## UB Analytical
      bn = log(n)
      n_i_guess = c( n_i, rep(0,Madd_grid[ii])  )
      ub_Bdd_mat[ii,b] = compute_UB_analytical(n,n_i_guess,Mguess,bn,alfa,FALSE)
      
      # Unbounded alphabet:
      Shat = sum( n_i )/n
      
      ## r-norm
      Sstar = ( sqrt( -log(beta)/(2*n) ) + sqrt( Shat + (-log(beta)/(2*n)) ) )^2
      rn = log( Sstar / (-log(1-alfa+beta)) ) + log(n) - log(log(n))
      ub_Ubd_mat[ii,b] = compute_UB_rnorm(n, alfa, beta, rn, Shat)
    }
    
    pb$tick()
  }
  
  # save 
  save_res = list( "ub_Bdd_mat" = ub_Bdd_mat,
                   "ub_Ubd_mat" = ub_Ubd_mat)
  
  if(save_exp)
    save(save_res, file = file_name)
  
}



## Final summary and plot ---------------------------------------------------

# load("save/Exp2_Zipfs_05.Rdat")
# ub_Bdd_mat = save_res$ub_Bdd_mat
# ub_Ubd_mat = save_res$ub_Ubd_mat

ub_Bdd  = apply(ub_Bdd_mat, 1, quantile, probs = c(0.025,0.5,0.975))
ub_Ubd  = apply(ub_Ubd_mat, 1, quantile, probs = c(0.025,0.5,0.975))

ymax = max(ub_Bdd[3,],ub_Ubd[3,]) * 1.05 #18*1e-3
# ymin = min(ub_Bdd[1,],ub_Ubd[1,]) * 1/12 #5*1e-3
# if(ymin < 0)
#   ymin = 0
ymin = 0
ylabs = round(seq(ymin*1e3,ymax*1e3,by = 1),1)




if(save_img)
  pdf(img_name)
par( mfrow = c(1,1), mar = c(3,3,1,0.5), mgp=c(2,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = expression(log[10](M[add])), ylab = "1000 * bound",
     xlim = range(logMadd_grid) , ylim = c(ymin,ymax), 
     main = paste0(" "),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs*1e-3, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
points( x = logMadd_grid, y = ub_Bdd[2,], 
        type = "l", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkgreen" ) 
points( x = logMadd_grid, y = ub_Ubd[2,], 
        type = "l", 
        lwd = 3, pch = 16, lty = 1,
        col = "darkorange" ) 
legend("bottomright",c("Bounded","Unbounded"), 
       lwd = 3, col = c("darkgreen","darkorange"))
if(save_img)
  dev.off()



