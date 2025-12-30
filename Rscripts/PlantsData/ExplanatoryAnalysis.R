source("../../R/Rfunctions.R")
# Load and read all --------------------------------------------------------------------
wd = "C:/Users/colom/BinomialCIs/Rscripts/PlantsData/"
setwd(wd)

mazz2016_Plants <- read.csv("C:/Users/colom/BinomialCIs/Rscripts/PlantsData/mazz2016_Plants.csv")
dim(mazz2016_Plants)
plants_names = mazz2016_Plants$X
data = t(mazz2016_Plants[,2:103])
n = nrow(data)


save_plot_all = TRUE


  
Nj = c(apply(data, 2, sum))
Kobs = length( which(Nj > 0) )
TabNj = tabulate(Nj, nbins = n) 

  
# Run
ExtrCurve = rarefaction.array(object = data, n_reorderings = 100, seed = 1234)

# Plot
if(save_plot_all)
  pdf("img/Plants_ExCurves.pdf")
par(mfrow = c(1,1),bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,1,0))
plot(x = 0, y = 0, type = "n",
     main = "Plants", xlab = "#obs.", ylab = "#species",
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

