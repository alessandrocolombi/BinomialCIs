# setwd("C:/Users/colom/bnp_upperbounds/Rscripts/species")
setwd("/home/lucia.paci/Lucia/Ale/bnp_upperbounds/Rscripts/species")

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
Bor  = 3 # nel suo codice, painsky ripete 10k volte, noi 5k?
Bnoi = Bor 
n = 1000
Rmax = 100; RmaxFD = 50
Mgrid = seq(10,1010,by = 50)
Nexp = length(Mgrid)
M_max = 500


# Zipfs: load data (from different repetitions) -------------------------------------------------------------------
# stop("Data have all been collected in SimStudyFeatures_zipfs_all.Rdat ")
s = 1.01
Nsim = 10
pain_mat <- lub_Mrk_mat <- lub_FD_mat <- oracle_mat <- c()
for(ii in 1:Nsim){
  file_name = paste0("../save/SimStudySpecies_zipfs_",ii,".Rdat")  
  load(file_name)
  pain_mat = cbind(pain_mat, save_res$pain_mat)
  lub_Mrk_mat = cbind(lub_Mrk_mat, save_res$lub_Mrk_mat)
  lub_FD_mat = cbind(lub_FD_mat,save_res$lub_FD_mat)
  oracle_mat = cbind(oracle_mat,save_res$oracle_mat)
}

# Save collected dataset
# save_res = list("lub_PP3_mat" = lub_PP3_mat,
#                 "lub_MixPois_mat" = lub_MixPois_mat,
#                 "lub_MixBin_mat"  = lub_MixBin_mat,
#                 "ub_Freq_mat"  = ub_Freq_mat,
#                 "oracle_mat"  = oracle_mat)
# file_name = paste0("../save/SimStudyFeatures_zipfs_all.Rdat")
# save(save_res, file = file_name)

# Zipfs: load data (all at once) -------------------------------------------------------------------

s = 1.01
var_gamma = 10
var_nb    = 10
Nrep = 5

load("../save/SimStudyFeatures_zipfs_all.Rdat")
lub_PP3_mat = save_res$lub_PP3_mat
lub_MixPois_mat = save_res$lub_MixPois_mat
ub_Freq_mat = save_res$ub_Freq_mat
lub_MixBin_mat = save_res$lub_MixBin_mat
oracle_mat = save_res$oracle_mat

# Compute final summary ---------------------------------------------------

ub_PYP = exp(apply(lub_Mrk_mat, 1, quantile, probs = c(0.025,0.5,0.975)))
ub_FD = exp(apply(lub_FD_mat, 1, quantile, probs = c(0.025,0.5,0.975)))
ub_pain = apply(pain_mat, 1, quantile, probs = c(0.025,0.5,0.975))
ub_oracle = apply(oracle_mat, 1, quantile, probs = c(0.025,0.5,0.975))

# Zipfs: Plot curves -------------------------------------------------------------------
ymax = (11/10) * (max(ub_oracle[2,],ub_PYP[2,],ub_FD[2,],ub_pain[2,]))
ymin = 0
ylabs = round(seq(0,ymax,length.out = 5),3)


save_img = FALSE
img_name = paste0("../img/SSSpecies.pdf")

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
points( x = Mgrid, y = ub_PYP[2,], 
        type = "b", lwd = 1, pch = 16, lty = 2, 
        col = "darkred" ) 
points( x = Mgrid, y = ub_FD[2,],
        type = "b", lwd = 1, pch = 16, lty = 2,
        col = "darkblue" )
points( x = Mgrid, y = ub_pain[2,], 
        type = "b", lwd = 1, pch = 16, lty = 2,
        col = "darkorange" ) 
## Add conf. intervals
# segments(x0 = Mgrid, x1 = Mgrid,
#          y0 = ub_oracle[1,], y1 = ub_oracle[3,],
#          col = "black", lwd = 1)
legend("bottomright",c("Oracle","PD","FDP","Freq."), 
       lwd = 3, col = c("black","darkred","darkblue","darkorange"))
if(save_img)
  dev.off()

# Zipfs: Coverage -------------------------------------------------------------------
Nrep_tot = ncol(oracle_mat)
PYP_cov     = rowSums( exp(lub_Mrk_mat) > oracle_mat)/Nrep_tot
FDP_cov     = rowSums(exp(lub_FD_mat) > oracle_mat)/Nrep_tot
Freq_cov    = rowSums(pain_mat > oracle_mat)/Nrep_tot

coverages = rbind(PYP_cov,FDP_cov,Freq_cov)
coverages = as.data.frame(coverages)
names(coverages) = paste0("M=",as.character(Mgrid))
coverages

