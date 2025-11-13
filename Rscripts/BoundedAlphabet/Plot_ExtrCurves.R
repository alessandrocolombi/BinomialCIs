
setwd("C:/Users/colom/BinomialCIs/Rscripts/BoundedAlphabet")

# Librerie ----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
suppressWarnings(suppressPackageStartupMessages(library(progress)))
suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../../R/Rfunctions.R")


# Plot extrapolation curves -----------------------------------------------

experiments = list(
  "Zipfs" = c(0.25,0.5,1.02,1.5),
  "Geom" = c(0.005,0.1,0.25),
  "Constant" = c(0.5,0.05,0.0001),
  "Uniform" = c(NA)
)

M = 5000
n = 2500
ngrid = ceiling(seq(1,n, length.out = 100))
BB = 50


save_plot_all = TRUE

i = 1
for(i in 1:length(experiments)){
  name_exp = names(experiments)[i]
  cat("\n Start: ",name_exp,"\n")
  if(save_plot_all)
    pdf(paste0("img/",name_exp,"/","ExtCurve_",name_exp,".pdf"))
  par(mfrow = c(1,1),bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,1,0))
  for(p in  experiments[[i]]){
    trimmed_param = get_first3digits(p,4)
    cat("\n p = ",p,"\n")
    prob_true_mat = matrix(0,nrow = BB, ncol = M) # (BB x M) matrix
    prob_true_mat = t(apply(prob_true_mat, 1, function(x) { sim_generic(name_exp,M,p) } ))
    
    NumObs = matrix(0,nrow = length(ngrid), ncol = BB)
    for(nn in 1:length(ngrid)){
      data = apply(prob_true_mat, 2, function(pj) {rbinom(n = BB, size = ngrid[nn], prob = pj)} ) # (BB x M) matrix
      NumObs[nn, ] = rowSums(data > 0)
    }
    
    ExtrCrv = apply(NumObs, 1, quantile, probs = c(0.025,0.5,0.975))
    plot(x = 0, y = 0, type = "n",
         main = paste0(name_exp," - ", trimmed_param), xlab = "#obs.", ylab = "#distincts",
         ylim = c(0,max(ExtrCrv)+1),
         xlim = c(0,n+1),
         pch = 1) # init plot
    polygon( c(ngrid, rev(ngrid)),
             c(ExtrCrv[1,], rev(ExtrCrv[3,])),
             col = "grey75",
             border = NA) # plot in-sample bands
    points(x = ngrid, y = ExtrCrv[2,], type = "l", lwd = 3) # plot mean obs
  }
  if(save_plot_all)
    dev.off()
}

  













