setwd("C:/Users/colom/bnp_upperbounds/Rscripts/features")

# Librerie ----------------------------------------------------------------

suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
suppressWarnings(suppressPackageStartupMessages(library(progress)))
suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../../R/Rfunctions.R")

# Funzioni ----------------------------------------------------------------

# Common parameters -------------------------------------------------------
alfa = 0.05
Bor  = 500 # nel suo codice, painsky ripete 10k volte, noi 5k?
Bnoi <- Bfr <- 500
n = 1000
Rmax = 100; 
Mgrid = seq(10,1010,by = 50)
Nexp = length(Mgrid)
M_max = 500


# Zipfs: load data -------------------------------------------------------------------

s = 1.01
var_gamma = 10
var_nb    = 10
Nrep = 5

lub_PP3_mat <- lub_MixPois_mat <- ub_Freq_mat <- lub_MixBin_mat <- oracle_mat <- c()
for(ii in 1:Nrep){
  file_name = paste0("../save/SimStudyFeatures_zipfs_",ii,".Rdat")  
  load(file_name)
  lub_PP3_mat = cbind(lub_PP3_mat, save_res$lub_PP3_mat)
  lub_MixPois_mat = cbind(lub_MixPois_mat, save_res$lub_MixPois_mat)
  ub_Freq_mat = cbind(ub_Freq_mat,save_res$ub_Freq_mat)
  lub_MixBin_mat = cbind(lub_MixBin_mat,save_res$lub_MixBin_mat)
  oracle_mat = cbind(oracle_mat,save_res$oracle_mat)
}



ub_PP3 = exp(apply(lub_PP3_mat, 1, quantile, probs = c(0.025,0.5,0.975)))
ub_MixPois = exp(apply(lub_MixPois_mat, 1, quantile, probs = c(0.025,0.5,0.975)))
ub_Freq = apply(ub_Freq_mat, 1, quantile, probs = c(0.025,0.5,0.975))
ub_MixBin  = exp(apply(lub_MixBin_mat, 1, quantile, probs = c(0.025,0.5,0.975)))
ub_oracle = apply(oracle_mat, 1, quantile, probs = c(0.025,0.5,0.975))

# Zipfs: Plot curves -------------------------------------------------------------------
ymax = (11/10) * (max(ub_oracle[2,],ub_PP3[2,],ub_MixPois[2,],ub_MixBin[2,],ub_Freq[2,]))
ymin = 0
ylabs = round(seq(0,ymax,length.out = 5),3)


save_img = FALSE
img_name = paste0("../img/SSFeatures.pdf")

if(save_img)
  pdf(img_name)
par( mfrow = c(1,1), mar = c(3,3,1,0.5), mgp=c(1.5,0.5,0), bty = "l" )
plot(0,0,  yaxt = "n",
     xlab = "M", ylab = " ",
     xlim = c(0,max(Mgrid) ) , ylim = c(ymin,ymax), 
     main = paste0("Zipfs"),
     type = "n")
grid(lty = 1,lwd = 1, col = "gray90" )
axis(side = 2, at = ylabs, 
     labels = ylabs, las = 1, 
     cex.axis = 1 )
points( x = Mgrid, y = ub_oracle[2,], 
        type = "b", lwd = 1, pch = 16, lty = 2,
        col = "black" ) 
points( x = Mgrid, y = ub_PP3[2,], 
        type = "b", lwd = 1, pch = 16, lty = 2, 
        col = "darkred" ) 
points( x = Mgrid, y = ub_MixBin[2,],
        type = "b", lwd = 1, pch = 16, lty = 2,
        col = "darkblue" )
points( x = Mgrid, y = ub_MixPois[2,], 
        type = "b", lwd = 1, pch = 16, lty = 2,
        col = "darkgreen" ) 
points( x = Mgrid, y = ub_Freq[2,], 
        type = "b", lwd = 1, pch = 16, lty = 2,
        col = "darkorange" ) 
## Add conf. intervals
# segments(x0 = Mgrid, x1 = Mgrid,
#          y0 = ub_oracle[1,], y1 = ub_oracle[3,],
#          col = "black", lwd = 1)
# segments(x0 = Mgrid, x1 = Mgrid,
#          y0 = ub_PP3[1,], y1 = ub_PP3[3,],
#          col = "darkred", lwd = 1)
# segments(x0 = Mgrid, x1 = Mgrid,
#          y0 = ub_MixBin[1,], y1 = ub_MixBin[3,],
#          col = "darkblue", lwd = 1)
# segments(x0 = Mgrid, x1 = Mgrid,
#          y0 = ub_MixPois[1,], y1 = ub_MixPois[3,],
#          col = "darkgreen", lwd = 1)
# segments(x0 = Mgrid, x1 = Mgrid,
#          y0 = ub_Freq[1,], y1 = ub_Freq[3,],
#          col = "darkorange", lwd = 1)
legend("bottomright",c("Oracle","PP","MixedPoisson","MixedBinomial","Freq."), 
       lwd = 3, col = c("black","darkred","darkgreen","darkblue","darkorange"))
if(save_img)
  dev.off()

# Zipfs: Coverage -------------------------------------------------------------------
Nrep_tot = ncol(oracle_mat)
PP3_cov     = rowSums( exp(lub_PP3_mat) > oracle_mat)/Nrep_tot
MixPois_cov = rowSums(exp(lub_MixPois_mat) > oracle_mat)/Nrep_tot
MixBin_cov  = rowSums(exp(lub_MixBin_mat) > oracle_mat)/Nrep_tot
Freq_cov    = rowSums(ub_Freq_mat > oracle_mat)/Nrep_tot

coverages = rbind(PP3_cov,MixPois_cov,MixBin_cov,Freq_cov)
coverages

