# Load and read all --------------------------------------------------------------------
wd = "C:/Users/colom/BinomialCIs/Rscripts/"
setwd(wd)

Rcpp::sourceCpp("../src/RcppFunctions.cpp")
source("../R/Rfunctions.R")

# 1) Beetles ---------------------------------------------------------------------
library(SpadeR)
data(DiversityData)
Nj = DiversityData$Inci[-1,]
names = rownames(DiversityData$Inci); names = names[-1]
n = c(DiversityData$Inci[1,1])
Kn = length(which(Nj > 0))

# Bounded
b_n <- log(n)
Mguess = 10 * Kn
Nj_guess = c(Nj, rep(0,Mguess - length(Nj) ))
alpha = 0.05
U_bounded <- compute_UB_analytical(n, Nj_guess, Mguess, b_n, alpha, FALSE)
U_bounded

# Unbounded
beta = 1e-5
Shat  <- sum(Nj) / n
Sstar <- ( sqrt( -log(beta) / (2 * n) ) +
             sqrt( Shat + (-log(beta) / (2 * n)) ) )^2
r_n   <- log( Sstar / (-log(1 - alpha + beta)) ) +
         log(n) - log(log(n))
U_unbounded <- compute_UB_rnorm(n, alpha, beta, r_n, Shat)
U_unbounded

# Coverage
C_hat <- coverage_ChaoJost(n, Nj)
C_hat



# 2) Soil extraction+ -----------------------------------------------------
library(SpadeR)
data(DiversityData)
DiversityData$Inci_count 
dim(DiversityData$Inci_count)
n = DiversityData$Inci_count[1,1]
pari = DiversityData$Inci_count[seq(2,36,by = 2),1]
dispari = DiversityData$Inci_count[seq(3,37,by = 2),1]
n_pairs = length(pari)
Nj = c()
for(j in 1:n_pairs){
  Nj = c(Nj, rep(pari[j], dispari[j] ) )
}
Kn = length(which(Nj > 0))

# Bounded
b_n <- log(n)
Mguess = 10 * Kn
Nj_guess = c(Nj, rep(0,Mguess - length(Nj) ))
alpha = 0.05
U_bounded <- compute_UB_analytical(n, Nj_guess, Mguess, b_n, alpha, FALSE)
U_bounded

# Unbounded
beta = 1e-5
Shat  <- sum(Nj) / n
Sstar <- ( sqrt( -log(beta) / (2 * n) ) +
             sqrt( Shat + (-log(beta) / (2 * n)) ) )^2
r_n   <- log( Sstar / (-log(1 - alpha + beta)) ) +
         log(n) - log(log(n))
U_unbounded <- compute_UB_rnorm(n, alpha, beta, r_n, Shat)
U_unbounded

# Coverage
C_hat <- coverage_ChaoJost(n, Nj)
C_hat


# 3) seed-bank records ---------------------------------------------
data(DiversityData)
dim(DiversityData$Inci_raw)
data = t(DiversityData$Inci_raw)
n = nrow(data)
Nj = c(apply(data, 2, sum))
Kobs = length( which(Nj > 0) )
TabNj = tabulate(Nj, nbins = n) 


# Bounded
b_n <- log(n)
Mguess = Kobs
Nj_guess = c(Nj, rep(0,Mguess - length(Nj) ))
alpha = 0.05
U_bounded <- compute_UB_analytical(n, Nj_guess, Mguess, b_n, alpha, FALSE)
U_bounded

# Unbounded
beta = 1e-5
Shat  <- sum(Nj) / n
Sstar <- ( sqrt( -log(beta) / (2 * n) ) +
             sqrt( Shat + (-log(beta) / (2 * n)) ) )^2
r_n   <- log( Sstar / (-log(1 - alpha + beta)) ) +
         log(n) - log(log(n))
U_unbounded <- compute_UB_rnorm(n, alpha, beta, r_n, Shat)
U_unbounded

# Coverage
C_hat <- coverage_ChaoJost(n, Nj)
C_hat


# Run
ExtrCurve = rarefaction.array(object = data, n_reorderings = 100, seed = 1234)

# Plot
save_plot_all = FALSE
if(save_plot_all)
  pdf("img/Seeds_ExCurves.pdf")
par(mfrow = c(1,1),bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,1,0))
plot(x = 0, y = 0, type = "n",
     main = "Seeds", xlab = "#obs.", ylab = "#species",
     ylim = c(0,Kobs+1),
     xlim = c(0,n+1),
     pch = 1) # init plot
polygon( c(1:n, rev(1:n)),
         c(ExtrCurve[1,], rev(ExtrCurve[3,])),
         col = "grey75",
         border = NA) # plot in-sample bands
points(x = 1:n, y = ExtrCurve[2,], type = "l", lwd = 3) # plot mean obs
if(save_plot_all)
  dev.off()




temp = sort(Nj,decreasing = TRUE, index.retur = TRUE )
ordered_Nj = temp$x
bp1 <- barplot(height = ordered_Nj)

ordering_index = temp$ix
ordered_types = plants_names[ordering_index]
sapply(ordered_types, function(x){trimm})

save_plot = FALSE
if(save_plot)
  pdf("img/Plants_Frequencies_all.pdf")

par(mfrow = c(1,1),bty = "l",  mgp=c(2.5,0.5,0), mar = c(2.5,3.5,2,0))
bars_names = as.character(0:4)
barplot( height = ordered_Nj, 
         las = 2, 
         col = "darkorange", border = NA,
         cex.names = 1,
         main = "", ylab = "#obs.", yaxt = "n", las = 1 )
axis( side = 2, at = round(seq(0, max(x), length.out = 5)), las = 1)
# text( x = bp1, y = par("usr")[3] - 0.02*max(Nj), 
#       labels = substr(ordered_types, 1, 3), 
#       srt = 45, adj = 1, xpd = TRUE, cex = 0.5 )
if(save_plot)
  dev.off()


# 4) soil ciliate species -------------------------------------------------
data(SimilarityMultData)
data = t(SimilarityMultData$Inci_raw)
n = nrow(data)
Nj = c(apply(data, 2, sum))
Kobs = length( which(Nj > 0) )
TabNj = tabulate(Nj, nbins = n) 


# Bounded
b_n <- log(n)
Mguess = 10 * Kobs
Nj_guess = c(Nj, rep(0,Mguess - length(Nj) ))
alpha = 0.05
U_bounded <- compute_UB_analytical(n, Nj_guess, Mguess, b_n, alpha, FALSE)
U_bounded

# Unbounded
beta = 1e-5
Shat  <- sum(Nj) / n
Sstar <- ( sqrt( -log(beta) / (2 * n) ) +
             sqrt( Shat + (-log(beta) / (2 * n)) ) )^2
r_n   <- log( Sstar / (-log(1 - alpha + beta)) ) +
         log(n) - log(log(n))
U_unbounded <- compute_UB_rnorm(n, alpha, beta, r_n, Shat)
U_unbounded

# Coverage
C_hat <- coverage_ChaoJost(n, Nj)
C_hat


# Run
ExtrCurve = rarefaction.array(object = data, n_reorderings = 100, seed = 1234)

# Plot
save_plot_all = FALSE
if(save_plot_all)
  pdf("img/Soil_ExCurves.pdf")
par(mfrow = c(1,1),bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,1,0))
plot(x = 0, y = 0, type = "n",
     main = "Soil", xlab = "#obs.", ylab = "#species",
     ylim = c(0,Kobs+1),
     xlim = c(0,n+1),
     pch = 1) # init plot
polygon( c(1:n, rev(1:n)),
         c(ExtrCurve[1,], rev(ExtrCurve[3,])),
         col = "grey75",
         border = NA) # plot in-sample bands
points(x = 1:n, y = ExtrCurve[2,], type = "l", lwd = 3) # plot mean obs
if(save_plot_all)
  dev.off()




temp = sort(Nj,decreasing = TRUE, index.retur = TRUE )
ordered_Nj = temp$x
bp1 <- barplot(height = ordered_Nj)

ordering_index = temp$ix
ordered_types = plants_names[ordering_index]
sapply(ordered_types, function(x){trimm})

save_plot = FALSE
if(save_plot)
  pdf("img/Plants_Frequencies_all.pdf")

par(mfrow = c(1,1),bty = "l",  mgp=c(2.5,0.5,0), mar = c(2.5,3.5,2,0))
bars_names = as.character(0:4)
barplot( height = ordered_Nj, 
         las = 2, 
         col = "darkorange", border = NA,
         cex.names = 1,
         main = "", ylab = "#obs.", yaxt = "n", las = 1 )
axis( side = 2, at = round(seq(0, max(x), length.out = 5)), las = 1)
# text( x = bp1, y = par("usr")[3] - 0.02*max(Nj), 
#       labels = substr(ordered_types, 1, 3), 
#       srt = 45, adj = 1, xpd = TRUE, cex = 0.5 )
if(save_plot)
  dev.off()


# 5) Capture-recapture ----------------------------------------------------

library(SpadeR)
data(ChaoSpeciesData)
?ChaoSpeciesData
data = t(ChaoSpeciesData$Inci_raw)
n = nrow(data)
Nj = c(apply(data, 2, sum))
Kobs = length( which(Nj > 0) )
TabNj = tabulate(Nj, nbins = n) 


# Bounded
b_n <- log(n)
Mguess = 10 * Kobs
Nj_guess = c(Nj, rep(0,Mguess - length(Nj) ))
alpha = 0.05
U_bounded <- compute_UB_analytical(n, Nj_guess, Mguess, b_n, alpha, FALSE)
U_bounded

# Unbounded
beta = 1e-5
Shat  <- sum(Nj) / n
Sstar <- ( sqrt( -log(beta) / (2 * n) ) +
             sqrt( Shat + (-log(beta) / (2 * n)) ) )^2
r_n   <- log( Sstar / (-log(1 - alpha + beta)) ) +
  log(n) - log(log(n))
U_unbounded <- compute_UB_rnorm(n, alpha, beta, r_n, Shat)
U_unbounded

# Coverage
C_hat <- coverage_ChaoJost(n, Nj)
C_hat


# Run
ExtrCurve = rarefaction.array(object = data, n_reorderings = 100, seed = 1234)

# Plot
save_plot_all = FALSE
if(save_plot_all)
  pdf("img/Recapture_ExCurves.pdf")
par(mfrow = c(1,1),bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,1,0))
plot(x = 0, y = 0, type = "n",
     main = "Recapture", xlab = "#obs.", ylab = "#species",
     ylim = c(0,Kobs+1),
     xlim = c(0,n+1),
     pch = 1) # init plot
polygon( c(1:n, rev(1:n)),
         c(ExtrCurve[1,], rev(ExtrCurve[3,])),
         col = "grey75",
         border = NA) # plot in-sample bands
points(x = 1:n, y = ExtrCurve[2,], type = "l", lwd = 3) # plot mean obs
if(save_plot_all)
  dev.off()




temp = sort(Nj,decreasing = TRUE, index.retur = TRUE )
ordered_Nj = temp$x
bp1 <- barplot(height = ordered_Nj)

ordering_index = temp$ix
ordered_types = plants_names[ordering_index]
sapply(ordered_types, function(x){trimm})

save_plot = FALSE
if(save_plot)
  pdf("img/Plants_Frequencies_all.pdf")

par(mfrow = c(1,1),bty = "l",  mgp=c(2.5,0.5,0), mar = c(2.5,3.5,2,0))
bars_names = as.character(0:4)
barplot( height = ordered_Nj, 
         las = 2, 
         col = "darkorange", border = NA,
         cex.names = 1,
         main = "", ylab = "#obs.", yaxt = "n", las = 1 )
axis( side = 2, at = round(seq(0, max(x), length.out = 5)), las = 1)
# text( x = bp1, y = par("usr")[3] - 0.02*max(Nj), 
#       labels = substr(ordered_types, 1, 3), 
#       srt = 45, adj = 1, xpd = TRUE, cex = 0.5 )
if(save_plot)
  dev.off()


# 6) Hong Kong Birds ----------------------------------------------------

library(SpadeR)
data(ChaoSharedData)
data = t(ChaoSharedData$Inci_raw)
dim(data)
n = nrow(data)
Nj = c(apply(data, 2, sum))
Kobs = length( which(Nj > 0) )
TabNj = tabulate(Nj, nbins = n) 


# Bounded
b_n <- log(n)
Mguess = 10 * Kobs
Nj_guess = c(Nj, rep(0,Mguess - length(Nj) ))
alpha = 0.05
U_bounded <- compute_UB_analytical(n, Nj_guess, Mguess, b_n, alpha, FALSE)
U_bounded

# Unbounded
beta = 1e-5
Shat  <- sum(Nj) / n
Sstar <- ( sqrt( -log(beta) / (2 * n) ) +
             sqrt( Shat + (-log(beta) / (2 * n)) ) )^2
r_n   <- log( Sstar / (-log(1 - alpha + beta)) ) +
         log(n) - log(log(n))
U_unbounded <- compute_UB_rnorm(n, alpha, beta, r_n, Shat)
U_unbounded

# Coverage
C_hat <- coverage_ChaoJost(n, Nj)
C_hat


# Run
ExtrCurve = rarefaction.array(object = data, n_reorderings = 100, seed = 1234)

# Plot
par(mfrow = c(1,1),bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,1,0))
plot(x = 0, y = 0, type = "n",
     main = "Birds", xlab = "#obs.", ylab = "#species",
     ylim = c(0,Kobs+1),
     xlim = c(0,n+1),
     pch = 1) # init plot
polygon( c(1:n, rev(1:n)),
         c(ExtrCurve[1,], rev(ExtrCurve[3,])),
         col = "grey75",
         border = NA) # plot in-sample bands
points(x = 1:n, y = ExtrCurve[2,], type = "l", lwd = 3) # plot mean obs




temp = sort(Nj,decreasing = TRUE, index.retur = TRUE )
ordered_Nj = temp$x
bp1 <- barplot(height = ordered_Nj)

ordering_index = temp$ix
ordered_types = plants_names[ordering_index]
sapply(ordered_types, function(x){trimm})

save_plot = FALSE
if(save_plot)
  pdf("img/Plants_Frequencies_all.pdf")

par(mfrow = c(1,1),bty = "l",  mgp=c(2.5,0.5,0), mar = c(2.5,3.5,2,0))
bars_names = as.character(0:4)
barplot( height = ordered_Nj, 
         las = 2, 
         col = "darkorange", border = NA,
         cex.names = 1,
         main = "", ylab = "#obs.", yaxt = "n", las = 1 )
axis( side = 2, at = round(seq(0, max(x), length.out = 5)), las = 1)
# text( x = bp1, y = par("usr")[3] - 0.02*max(Nj), 
#       labels = substr(ordered_types, 1, 3), 
#       srt = 45, adj = 1, xpd = TRUE, cex = 0.5 )
if(save_plot)
  dev.off()


