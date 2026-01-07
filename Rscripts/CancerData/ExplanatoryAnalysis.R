
# Load and read all --------------------------------------------------------------------
wd = "C:/Users/colom/BinomialCIs/Rscripts/CancerData/"
setwd(wd)
load("cancer_types.Rdat")
load("cancer_names_easy.Rdat")


# Accumulation curves -----------------------------------------------------

TabNj_list = vector("list", length = length(cancer_names))
ExtrCurve_list = vector("list", length = length(cancer_names))
N = 0

idx = 1
save_plot_all = FALSE
save_plot_individual = FALSE

if(save_plot_all)
  pdf("img/ExtrapolationCurvs_all.pdf")

par(mfrow = c(2,2),bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,1,0))
for(idx in 1:length(cancer_names)){
  cat("\n idx = ",idx,"\n")
  cancer_name = cancer_names[idx]
  filename = paste0(wd,"TCGA/",cancer_types[idx],"_targeted.RData")
  load(filename)
  Nj = c(apply(Z, 2, sum))
  n = nrow(Z)
  N = N + n
  Kobs = length( which(Nj > 0) )
  TabNj = c(length(which(Nj == 0)), tabulate(Nj, nbins = n) )
  TabNj_list[[idx]] = TabNj
  
  # Run
  ExtrCurve = rarefaction.array(object = Z, n_reorderings = 20, seed = 1234)
  ExtrCurve_list[[idx]] = ExtrCurve

  # Plot
  if(save_plot_individual)
    pdf(paste0("img/ExtrapolationCurvs_",cancer_types[idx],".pdf"))
  par(mfrow = c(1,1),bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,1,0))
  plot(x = 0, y = 0, type = "n",
       main = cancer_name, xlab = "#obs.", ylab = "#variants",
       ylim = c(0,Kobs+1),
       xlim = c(0,n+1),
       pch = 1) # init plot
  polygon( c(1:n, rev(1:n)),
           c(ExtrCurve[1,], rev(ExtrCurve[3,])),
           col = "grey75",
           border = NA) # plot in-sample bands
  points(x = 1:n, y = ExtrCurve[2,], type = "l", lwd = 3) # plot mean obs
  if(save_plot_individual)
    dev.off()
}
if(save_plot_all)
  dev.off()




# Exploratory plots (#obs & #variants) -------------------------------------------------------
nj = sapply(ExtrCurve_list, function(x){ncol(x)})
Ncancers = length(nj)
Kobsj = sapply(TabNj_list, function(x) sum( x[-1] ) )

ordered_nj = sort(nj,decreasing = TRUE, index.retur = TRUE )
nj = ordered_nj$x
bp1 <- barplot(height = nj)

ordering_index = ordered_nj$ix
ordered_types = cancer_types[ordering_index]
ordered_Kj = Kobsj[ordering_index]
bp2 <- barplot(height = ordered_Kj)

# Saved in custom size: 6x14 in
par(mfrow = c(1,2), mgp=c(2.5,0.5,0), mar = c(2.5,3.5,1,0))
barplot( height = nj, 
         names.arg = "", las = 2, col = "darkred", border = NA,
         main = " ", ylab = "#obs.", yaxt = "n" )
axis( side = 2, at = seq(0, max(nj), by=100), las = 1)
text( x = bp1, y = par("usr")[3] - 0.02*max(nj), 
      labels = ordered_types, 
      srt = 45, adj = 1, xpd = TRUE, cex = 0.5 )

barplot( height = ordered_Kj, names.arg = "",
         las = 2, cex.names = 0.8, col = "darkblue", border = NA,                
         main = " ", ylab = "#variants", yaxt = "n" )
axis( side = 2, at = seq(0, max(ordered_Kj), by=1000),las = 1 )
text( x = bp2, y = par("usr")[3] - 0.02*max(Kobsj),
      labels = ordered_types, 
      srt = 45, adj = 1, xpd = TRUE, cex = 0.5 )


# Frequencies -------------------------------------------------------------

save_plot = TRUE
if(save_plot)
  pdf("img/Frequencies_all.pdf")

par(mfrow = c(2,2),bty = "l",  mgp=c(2.5,0.5,0), mar = c(2.5,3.5,2,0))
for(idx in 1:Ncancers){
  x = TabNj_list[[idx]]
  cancer_name = cancer_names[idx]
  # positive_indices <- which(x > 0)
  # max_index <- max(positive_indices)
  # x = x[1:max_index]
  x = x[1:5]
  
  bars_names = as.character(0:4)
  barplot( height = x, 
           names.arg = bars_names, las = 2, 
           col = "darkorange", border = NA,
           cex.names = 1,
           main = cancer_name, ylab = "#obs.", yaxt = "n", las = 1 )
  axis( side = 2, at = round(seq(0, max(x), length.out = 5)), las = 1)

}

if(save_plot)
  dev.off()


# Bounded or not? ---------------------------------------------------------

# Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
# source("../../R/Rfunctions.R")

d = length(cancer_names)
idx = 5;
Mstar <- rep(-1,d)
for(idx in 1:d){
  cancer_name = cancer_names[idx]
  cat("\n idx = ",idx,"; cancer name: ",cancer_name," ")
  
  exp_name = paste0("Cancer_",cancer_types[idx])
  
  filename = paste0(wd,"TCGA/",cancer_types[idx],"_targeted.RData")
  load(filename)
  Nj = c(apply(Z, 2, sum))
  n = nrow(Z)
  Kobs = length( which(Nj > 0) )
  TabNj = c(length(which(Nj == 0)), tabulate(Nj, nbins = n) )
  S<-Shat<- sum( Nj )/n
  library(VGAM)
  z <- S*n - log(20)
  W <- ifelse(z < 700, VGAM::lambertW(exp(z)), z - log(z))
  cat("...",W,"\n")
  Mstar[idx] = W
}


par(mfrow = c(1,1), mgp=c(2.5,0.5,0), mar = c(2.5,3.5,1,0), bty = "l")
plot(Mstar)




