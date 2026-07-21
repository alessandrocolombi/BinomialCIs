
setwd("C:/Users/colom/BinomialCIs/Rscripts/SimulationStudy_firstSub")

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
  "Zipfs" = c(1.05,0.85,0.75)
)

M = 1500
n = 1000
ngrid = ceiling(seq(1,n, length.out = 100))
BB = 50
save_img = TRUE
width = 8; height = 6
cex.lab = 2
cex.axis = 2

seed = 231231
set.seed(seed)


i = 1
for(i in 1:length(experiments)){
  name_exp = names(experiments)[i]
  cat("\n Start: ",name_exp,"\n")
  ii = 1
  for(ii in  seq_along(experiments[[i]])){
    p = experiments[[i]][ii]
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
    
    img_name = paste0("img/",name_exp,"/AccCurve_",trimmed_param,".pdf")
    
    if(save_img)
      pdf(img_name, width = width, height = height)
    if(ii == 1){
      par( mfrow = c(1,1), mar = c(3,6,1,1), mgp=c(4,0.75,0), bty = "l", las = 1, cex.lab = cex.lab )
      plot(x = 0, y = 0, 
           type = "n", xaxt = "n", yaxt = "n",
           main = "", 
           xlab = "", ylab = "#symbols",
           ylim = c(0,max(ExtrCrv)+1),
           xlim = c(0,n+1),
           pch = 1) # init plot
    } else{
      par( mfrow = c(1,1), mar = c(3,4.5,1,1), mgp=c(3,0.75,0), bty = "l", las = 1, cex.lab = cex.lab )
      plot(x = 0, y = 0, 
           type = "n", xaxt = "n", yaxt = "n",
           main = "", 
           xlab = "", ylab = "",
           ylim = c(0,max(ExtrCrv)+1),
           xlim = c(0,n+1),
           pch = 1) # init plot
    }
    axis(2, cex.axis = cex.axis)
    axis(1, cex.axis = cex.axis);mtext("n", side = 1, line = 2, cex = cex.axis)
    polygon( c(ngrid, rev(ngrid)),
             c(ExtrCrv[1,], rev(ExtrCrv[3,])),
             col = "grey75",
             border = NA) # plot in-sample bands
    points(x = ngrid, y = ExtrCrv[2,], type = "l", lwd = 5) # plot mean obs
    if(save_img)
      dev.off()
    
  }

}















