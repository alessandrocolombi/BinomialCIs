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
# Bor  = 3 # nel suo codice, painsky ripete 10k volte, noi 5k?
# Bnoi = Bor 
n = 1000
Rmax = 100; RmaxFD = 50
Mgrid = seq(10,1010,by = 50)
Nexp = length(Mgrid)
M_max = 500


# Zipfs: load data (from different repetitions) -------------------------------------------------------------------
stop("Data have all been collected in SimStudyFeatures_zipfs_all.Rdat ")
s = 1.01
Nsim = 10
sim_idx = 11:30
sim_idx = sim_idx[-11]

pain_mat <- lub_Mrk_mat <- lub_FD_mat <- oracle_mat <- c()
for(ii in sim_idx){
  file_name = paste0("../save/SimStudySpecies_zipfs_",ii,".Rdat")  
  load(file_name)
  pain_mat = cbind(pain_mat, save_res$pain_mat)
  lub_Mrk_mat = cbind(lub_Mrk_mat, save_res$lub_Mrk_mat)
  lub_FD_mat = cbind(lub_FD_mat,save_res$lub_FD_mat)
  oracle_mat = cbind(oracle_mat,save_res$oracle_mat)
}

# # Save collected dataset
# save_res = list("lub_Mrk_mat" = lub_Mrk_mat,
#                 "lub_FD_mat" = lub_FD_mat,
#                 "pain_mat"  = pain_mat,
#                 "oracle_mat"  = oracle_mat)
# file_name = paste0("../save/SimStudySpecies_zipfs_all.dat")
# save(save_res, file = file_name)

# Zipfs: load data (all at once) -------------------------------------------------------------------

s = 1.01
var_gamma = 10
var_nb    = 10

load("../save/SimStudySpecies_zipfs_all.dat")
lub_Mrk_mat = save_res$lub_Mrk_mat
lub_FD_mat = save_res$lub_FD_mat
pain_mat = save_res$pain_mat
oracle_mat = save_res$oracle_mat

# Compute final summary ---------------------------------------------------

ub_PYP = exp(apply(lub_Mrk_mat, 1, quantile, probs = c(0.025,0.5,0.975)))
ub_FD = exp(apply(lub_FD_mat, 1, quantile, probs = c(0.025,0.5,0.975)))
ub_pain = apply(pain_mat, 1, quantile, probs = c(0.025,0.5,0.975))
ub_oracle = apply(oracle_mat, 1, quantile, probs = 1-alfa)


# Zipfs: Plot curves -------------------------------------------------------------------
ymax = (11/10) * (max(ub_oracle,ub_PYP[2,],ub_FD[2,],ub_pain[2,]))
ymin = 0
ylabs = round(seq(0,ymax,length.out = 5),3)


save_img = TRUE
img_name = paste0("../img/SSSpecies_zipfs.pdf")

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
points( x = Mgrid, y = ub_oracle, 
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

file_name = paste0("../save/SimStudySpecies_coverages_zipfs_all.Rdat")
save(coverages, file = file_name)






# Neg. Bin.: load data (from different repetitions) -------------------------------------------------------------------
stop("Data have all been collected in SimStudyFeatures_NegBin_all.Rdat ")
l = 1; r = 0.003
Nsim = 10
sim_idx = 1:20
# sim_idx = sim_idx[-11]

pain_mat <- lub_Mrk_mat <- lub_FD_mat <- oracle_mat <- c()
for(ii in sim_idx){
  file_name = paste0("../save/SimStudySpecies_NegBib_",ii,".Rdat")  
  load(file_name)
  pain_mat = cbind(pain_mat, save_res$pain_mat)
  lub_Mrk_mat = cbind(lub_Mrk_mat, save_res$lub_Mrk_mat)
  lub_FD_mat = cbind(lub_FD_mat,save_res$lub_FD_mat)
  oracle_mat = cbind(oracle_mat,save_res$oracle_mat)
}

# Save collected dataset
# save_res = list("lub_Mrk_mat" = lub_Mrk_mat,
#                 "lub_FD_mat" = lub_FD_mat,
#                 "pain_mat"  = pain_mat,
#                 "oracle_mat"  = oracle_mat)
# file_name = paste0("../save/SimStudySpecies_NegBin_all.dat")
# save(save_res, file = file_name)

# NegBin: load data (all at once) -------------------------------------------------------------------

load("../save/SimStudySpecies_NegBin_all.Rdat")
lub_Mrk_mat = save_res$lub_Mrk_mat
lub_FD_mat = save_res$lub_FD_mat
pain_mat = save_res$pain_mat
oracle_mat = save_res$oracle_mat

# Compute final summary ---------------------------------------------------

ub_PYP = exp(apply(lub_Mrk_mat, 1, quantile, probs = c(0.025,0.5,0.975)))
ub_FD = exp(apply(lub_FD_mat, 1, quantile, probs = c(0.025,0.5,0.975)))
ub_pain = apply(pain_mat, 1, quantile, probs = c(0.025,0.5,0.975))
ub_oracle = apply(oracle_mat, 1, quantile, probs = 1-alfa)


# NegBin: Plot curves -------------------------------------------------------------------
ymax = (11/10) * (max(ub_oracle,ub_PYP[2,],ub_FD[2,],ub_pain[2,]))
ymin = 0
ylabs = round(seq(0,ymax,length.out = 5),3)


save_img = TRUE
img_name = paste0("../img/SSSpecies_NegBin.pdf")

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
points( x = Mgrid, y = ub_oracle, 
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

# NegBin: Coverage -------------------------------------------------------------------
Nrep_tot = ncol(oracle_mat)
PYP_cov     = rowSums( exp(lub_Mrk_mat) > oracle_mat)/Nrep_tot
FDP_cov     = rowSums(exp(lub_FD_mat) > oracle_mat)/Nrep_tot
Freq_cov    = rowSums(pain_mat > oracle_mat)/Nrep_tot

coverages = rbind(PYP_cov,FDP_cov,Freq_cov)
coverages = as.data.frame(coverages)
names(coverages) = paste0("M=",as.character(Mgrid))
coverages

file_name = paste0("../save/SimStudySpecies_coverages_NegBin_all.dat")
save(coverages, file = file_name)

