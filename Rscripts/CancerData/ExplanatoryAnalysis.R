# Funzioni Lorenzo --------------------------------------------------------
convert_features_list <- function(feature_matrix){
  
  feat_list <- vector("list", nrow(feature_matrix))
  
  for (i in 1:nrow(feature_matrix)){
    feat_list[[i]] <- which(feature_matrix[i,]==1, arr.ind = TRUE)
  }
  
  return (feat_list)
}
rarefaction.array <- function(object, n_reorderings = 1, seed = 1234) {
  
  feature_list <- convert_features_list(object)
  n <- nrow(object)
  
  if (n_reorderings == 1){
    
    rare_curve <- sapply(1:n, function(i) length(unique(unlist(feature_list[1:i]))) )
    
  } else {
    
    rare_curve <- matrix(NA, nrow = n_reorderings, ncol = n)
    
    for (j in 1:n_reorderings){
      
      f_list <- sample(feature_list)
      
      rare_curve[j, ] <- sapply(1:n, function(i) length(unique(unlist(f_list[1:i]))) )
    }
    
    # rare_curve <- colMeans(rare_curve)
    rare_curve_qnt = apply(rare_curve,2,quantile, prob = c(0.025,0.5,0.975))
  }
  
  return(rare_curve_qnt)
  
}



# Load and read all --------------------------------------------------------------------


load("cancer_types.Rdat")
wd = "C:/Users/colom/BinomialCIs/Rscripts/CancerData/"
load("cancer_names_easy.Rdat")

TabNj_list = vector("list", length = length(cancer_names))
ExtrCurve_list = vector("list", length = length(cancer_names))

idx = 1
save_plot = FALSE

if(save_plot)
  pdf("img/ExtrapolationCurvs_all.pdf")

par(mfrow = c(2,2),bty = "l",  mgp=c(1.5,0.5,0), mar = c(2.5,2.5,1,0))
for(idx in 1:length(cancer_names)){
  cat("\n idx = ",idx,"\n")
  cancer_name = cancer_names[idx]
  filename = paste0(wd,"TCGA/",cancer_types[idx],"_targeted.RData")
  load(filename)
  Nj = c(apply(Z, 2, sum))
  n = nrow(Z)
  Kobs = length( which(Nj > 0) )
  TabNj = c(length(which(Nj == 0)), tabulate(Nj, nbins = n) )
  TabNj_list[[idx]] = TabNj
  
  # Run
  ExtrCurve = rarefaction.array(object = Z, n_reorderings = 20, seed = 1234)
  ExtrCurve_list[[idx]] = ExtrCurve
  
  # Plot
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
}
if(save_plot)
  dev.off()




# Exploratory plots (#obs & #variants) -------------------------------------------------------
Ncancers = length(nj)
nj = sapply(ExtrCurve_list, function(x){ncol(x)})
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

