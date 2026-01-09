# Load and read all --------------------------------------------------------------------
wd = "C:/Users/colom/BinomialCIs/Rscripts/BBSData"
setwd(wd)

Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../../R/Rfunctions.R")

load(paste0("data/Data2019_allRoutes.Rdat"))
data = incidence_matrix
n = nrow(data)
Nj = c(apply(data, 2, sum))
Kobs = length( which(Nj > 0) )
TabNj = tabulate(Nj, nbins = n) 

# Run
ExtrCurve = rarefaction.array(object = data, n_reorderings = 20, seed = 1234)

# Plot
par(mfrow = c(1,1),bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,1,0))
plot(x = 0, y = 0, type = "n",
     main = "", xlab = "#obs.", ylab = "#species",
     ylim = c(0,Kobs+1),
     xlim = c(0,n+1),
     pch = 1) # init plot
polygon( c(1:n, rev(1:n)),
         c(ExtrCurve[1,], rev(ExtrCurve[3,])),
         col = "grey75",
         border = NA) # plot in-sample bands
points(x = 1:n, y = ExtrCurve[2,], type = "l", lwd = 3) # plot mean obs


