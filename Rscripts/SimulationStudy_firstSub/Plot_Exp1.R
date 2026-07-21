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
name = "Zipfs" # choose the name here
params = c(0.25,0.5,1.02)

# Common parameters -------------------------------------------------------
alfa = 0.05
beta = 1e-5
Brep = 100; Bor = 100
seed = 42
set.seed(seed)
width = 8; height = 6
cex.lab = 2
cex.axis = 2

# n fix, M varies ---------------------------------------------------------
save_img = TRUE

n = 2000
Mgrid_min =  2000 
Mgrid_max =  10000
Mgrid_step = 250

Mgrid = seq(Mgrid_min,Mgrid_max,by = Mgrid_step)
Nexp = length(Mgrid)

ii = 3
for(ii in seq_along(params)){
  
  ## Get param
  param = params[ii]
  trimmed_param = get_first3digits(param,4)
  exp_name = paste0("Exp1_nfix_",name,"_",trimmed_param)
  file_name = paste0("save/",exp_name,".Rdat")
  img_name = paste0("img/",name,"/",exp_name,".pdf")
  cov_name = paste0("save/Coverage/",name,"/",exp_name,".Rdat")
  
  ## Load data
  load(file_name)
  ub_Bnf_mat = save_res$ub_Bnf_mat
  ub_Bdd_mat = save_res$ub_Bdd_mat
  ub_Ubd_mat = save_res$ub_Ubd_mat
  
  ## Read and plot
  ub_Bnf = apply(ub_Bnf_mat, 1, quantile, probs = c(0.025,0.5,0.975))
  ub_Bdd = apply(ub_Bdd_mat, 1, quantile,    probs = c(0.025,0.5,0.975))
  ub_Ubd = apply(ub_Ubd_mat, 1, quantile,    probs = c(0.025,0.5,0.975))
  
  ymax = 8*1e-3;ymin = 0
  ylabs = round(seq(ymin*1e3,ymax*1e3,by = 1),1)
  
  
  if(save_img)
    pdf(img_name, width = width, height = height)
  if(ii == 1){
    par( mfrow = c(1,1), mar = c(3.5,4,1,1), mgp=c(2.5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
    plot(0,0,  yaxt = "n", xaxt = "n",
         xlab = "", ylab = "1000 * bound",
         xlim = c(min(Mgrid),max(Mgrid) ) , ylim = c(ymin,ymax), 
         main = paste0(" "),
         type = "n")
  }
  else{
    par( mfrow = c(1,1), mar = c(3.5,2,1,1), mgp=c(2,1,0), bty = "l", las = 1, cex.lab = cex.lab )
    plot(0,0,  yaxt = "n", xaxt = "n",
         xlab = "", ylab = " ",
         xlim = c(min(Mgrid),max(Mgrid) ) , ylim = c(ymin,ymax), 
         main = paste0(" "),
         type = "n")
  }
  grid(lty = 1,lwd = 1, col = "gray90" )
  axis(side = 2, at = ylabs*1e-3, 
       labels = ylabs, las = 1, 
       cex.axis = cex.axis )
  axis(1, cex.axis = cex.axis);mtext("M", side = 1, line = 2.5, cex = cex.axis)
  points( x = Mgrid, y = ub_Bnf[2,], 
          type = "l", 
          lwd = 5, pch = 16, lty = 1,
          col = "darkred" ) 
  points( x = Mgrid, y = ub_Bdd[2,], 
          type = "l", 
          lwd = 5, pch = 16, lty = 1,
          col = "darkgreen" ) 
  points( x = Mgrid, y = ub_Ubd[2,],
          type = "l",
          lwd = 5, pch = 16, lty = 1,
          col = "darkorange" )
  legend("bottomright",c("Bonferroni","Bounded","Unbounded"), 
         lwd = 5, col = c("darkred","darkgreen","darkorange"), cex = 1.5)
  if(save_img)
    dev.off()
}

# M fix, n varies ---------------------------------------------------------
M = 5000
Ngrid_max =  5000
Ngrid_min =  1000
Ngrid_step = 250

Ngrid = seq(Ngrid_min,Ngrid_max,by = Ngrid_step)
Nexp = length(Ngrid)

ii = 1
for(ii in seq_along(params)){
  
  ## Get param
  param = params[ii]
  trimmed_param = get_first3digits(param,4)
  exp_name = paste0("Exp1_Mfix_",name,"_",trimmed_param)
  file_name = paste0("save/",exp_name,".Rdat")
  img_name = paste0("img/",name,"/",exp_name,".pdf")
  cov_name = paste0("save/Coverage/",name,"/",exp_name,".Rdat")
  
  ## Load data
  load(file_name)
  ub_Bnf_mat = save_res$ub_Bnf_mat
  ub_Bdd_mat = save_res$ub_Bdd_mat
  ub_Ubd_mat = save_res$ub_Ubd_mat
  
  ## Read and plot
  ub_Bnf = apply(ub_Bnf_mat, 1, quantile, probs = c(0.025,0.5,0.975))
  ub_Bdd = apply(ub_Bdd_mat, 1, quantile,    probs = c(0.025,0.5,0.975))
  ub_Ubd = apply(ub_Ubd_mat, 1, quantile,    probs = c(0.025,0.5,0.975))
  
  ymax = 14*1e-3;ymin = 0
  ylabs = round(seq(ymin*1e3,ymax*1e3,by = 1),1)
  
  
  if(save_img)
    pdf(img_name, width = width, height = height)
  if(ii == 1){
    par( mfrow = c(1,1), mar = c(3.5,4.5,1,1), mgp=c(2.75,1,0), bty = "l", las = 1, cex.lab = cex.lab )
    plot(0,0,  yaxt = "n", xaxt = "n",
         xlab = "", ylab = "1000 * bound",
         xlim = c(min(Ngrid),max(Ngrid) ) , ylim = c(ymin,ymax), 
         main = paste0(" "),
         type = "n")
  }
  else{
    par( mfrow = c(1,1), mar = c(3.5,2.75,1,1), mgp=c(2,1,0), bty = "l", las = 1, cex.lab = cex.lab )
    plot(0,0,  yaxt = "n", xaxt = "n",
         xlab = "", ylab = " ",
         xlim = c(min(Ngrid),max(Ngrid) ) , ylim = c(ymin,ymax), 
         main = paste0(" "),
         type = "n")
  }
  grid(lty = 1,lwd = 1, col = "gray90" )
  axis(side = 2, at = ylabs*1e-3, 
       labels = ylabs, las = 1, 
       cex.axis = cex.axis )
  axis(1, cex.axis = cex.axis);mtext("n", side = 1, line = 2.5, cex = cex.axis)
  points( x = Ngrid, y = ub_Bnf[2,], 
          type = "l", 
          lwd = 5, pch = 16, lty = 1,
          col = "darkred" ) 
  points( x = Ngrid, y = ub_Bdd[2,], 
          type = "l", 
          lwd = 5, pch = 16, lty = 1,
          col = "darkgreen" ) 
  points( x = Ngrid, y = ub_Ubd[2,],
          type = "l",
          lwd = 5, pch = 16, lty = 1,
          col = "darkorange" )
  legend("topright",c("Bonferroni","Bounded","Unbounded"), 
         lwd = 5, col = c("darkred","darkgreen","darkorange"), cex = 1.5)
  if(save_img)
    dev.off()
}


