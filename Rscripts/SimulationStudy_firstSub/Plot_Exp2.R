
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
name = "Zipfs"
params = c(0.25,0.5,1.02)
# Common parameters -------------------------------------------------------
alfa = 0.05
beta = 1e-5
Brep = 100; Bor = 100
seed = 42
set.seed(seed)

save_img = TRUE
width = 8; height = 6
cex.lab = 2
cex.axis = 2


n = 1000
M = 5000 # this is Mtrue, the true size of the alphabet
Madd_grid = c(1e0,1e1,1e2,1e3,1e4,1e5,1e6)
logMadd_grid = log10(Madd_grid)
Mguess_grid = M + Madd_grid
Nexp = length(Madd_grid)


ii = 1
for(ii in seq_along(params)){
  
  ## Get param
  param = params[ii]
  trimmed_param = get_first3digits(param,4)
  exp_name = paste0("Exp2_",name,"_",trimmed_param)
  file_name = paste0("save/",exp_name,".Rdat")
  img_name = paste0("img/",name,"/",exp_name,".pdf")
  
  ## Load data
  load(file_name)
  ub_Bdd_mat = save_res$ub_Bdd_mat
  ub_Ubd_mat = save_res$ub_Ubd_mat
  
  ## Read and plot
  ub_Bdd  = apply(ub_Bdd_mat, 1, quantile, probs = c(0.025,0.5,0.975))
  ub_Ubd  = apply(ub_Ubd_mat, 1, quantile, probs = c(0.025,0.5,0.975))
  
  ymax = 17*1e-3; ymin = 0*1e-3
  ylabs = round(seq(ymin*1e3,ymax*1e3,by = 1),1)
  
  
  if(save_img)
    pdf(img_name, width = width, height = height)
  if(ii == 1){
    par( mfrow = c(1,1), mar = c(4,4.5,1,1), mgp=c(3,1,0), bty = "l", las = 1, cex.lab = cex.lab )
    plot(0,0,  yaxt = "n", xaxt = "n",
         xlab = expression(log[10](M[add])), ylab = "1000 * bound",
         xlim = range(logMadd_grid) , ylim = c(ymin,ymax), 
         main = paste0(" "),
         type = "n")
  }else{
    par( mfrow = c(1,1), mar = c(4,3,1,1),  mgp=c(3,1,0), bty = "l", las = 1, cex.lab = cex.lab )
    plot(0,0,  yaxt = "n", xaxt = "n",
         xlab = expression(log[10](M[add])), ylab = " ",
         xlim = range(logMadd_grid) , ylim = c(ymin,ymax), 
         main = paste0(" "),
         type = "n")
  }
  grid(lty = 1,lwd = 1, col = "gray90" )
  axis(side = 2, at = ylabs*1e-3, 
       labels = ylabs, las = 1, 
       cex.axis = cex.axis )
  axis(1, cex.axis = cex.axis)
  points( x = logMadd_grid, y = ub_Bdd[2,], 
          type = "l", 
          lwd = 5, pch = 16, lty = 1,
          col = "darkgreen" ) 
  points( x = logMadd_grid, y = ub_Ubd[2,], 
          type = "l", 
          lwd = 5, pch = 16, lty = 1,
          col = "darkorange" ) 
  legend("bottomright",c("Bounded","Unbounded"), 
         lwd = 5, col = c("darkgreen","darkorange"),cex = 1.5)
  if(save_img)
    dev.off()
}




