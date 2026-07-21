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

save_img = TRUE
width = 12; height = 6
cex.lab = 2
cex.axis = 2

img_name = "img/Data2019_allRoutes_AccCurve.pdf"
# Plot
if(save_img)
  pdf(img_name, width = width, height = height)
par(mfrow = c(1,1), bty = "l",  mar = c(3.5,5.8,1,1), mgp=c(4,1,0), las = 1, cex.lab = cex.lab)
plot(x = 0, y = 0, type = "n",
     xaxt = "n", yaxt = "n",
     main = "", xlab = "", ylab = "#species",
     ylim = c(0,Kobs+1),
     xlim = c(0,n+1),
     pch = 1) # init plot
axis(2, cex.axis = cex.axis)
axis(1, cex.axis = cex.axis);mtext("n", side = 1, line = 2, cex = cex.axis)
polygon( c(1:n, rev(1:n)),
         c(ExtrCurve[1,], rev(ExtrCurve[3,])),
         col = "grey75",
         border = NA) # plot in-sample bands
points(x = 1:n, y = ExtrCurve[2,], type = "l", lwd = 5) # plot mean obs
if(save_img)
  dev.off()

