# Wd and functions ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100)
choose_wd = wd_vec[1] # <--- modify here
wd = paste0(choose_wd,"BinomialCIs/Rscripts/CancerData/")
setwd(wd)


suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))

load("cancer_types.Rdat")
load("cancer_names_easy.Rdat")
d = length(cancer_names)

num_cores = d
idxs = 1:d


# Plot in linear scale ----------------------------------------------------

save_img = FALSE
width = 8; height = 6
cex.lab = 2
cex.axis = 2

idx = 15
for(idx in idxs) {
  cancer_name = cancer_names[idx]
  filename = paste0(wd,"TCGA/",cancer_types[idx],"_targeted.RData")
  load(filename)
  cat("\n idx = ",idx,"; cancer name: ",cancer_name,"\n")
  
  seed = 42
  set.seed(seed)
  
  exp_name = paste0("Cancer_",cancer_types[idx])
  file_name = paste0("save/",exp_name,".Rdat")
  img_name_1 = paste0("img/",exp_name,"_len.pdf")
  
  Nj = c(apply(Z, 2, sum))
  n = nrow(Z)
  Kobs = length( which(Nj > 0) )
  TabNj = c(length(which(Nj == 0)), tabulate(Nj, nbins = n) )
  train_prop = seq(0.05,1,by = 0.05)
  Nexp = length(train_prop)
  Nrep_max = 200
  Nrep = ceiling(seq(Nrep_max,1,length.out = Nexp))
  alfa = 0.05
  beta = 1e-5
  
  load(file_name)  
  ub_unb_list = res$ub_unb_list
  ub_bdd_list = res$ub_bdd_list
  ub_unb  = lapply(ub_unb_list, quantile, probs = c(0.025,0.5,0.975)); ub_unb = do.call(cbind, ub_unb)
  ub_bdd  = lapply(ub_bdd_list, quantile, probs = c(0.025,0.5,0.975)); ub_bdd = do.call(cbind, ub_bdd)
  
  ymax = 1#max( c( ub_unb[3,],ub_bdd[3,])) * 1.05
  ymin = 0 #min( c( nub_unb[1,],nub_bdd[1,])) 
  ylabs = round(seq(ymin,ymax,length.out = 5),1)
  
  n_train_vec = ceiling( n*train_prop )
  
  
  ub_unb[which(ub_unb > 1)] = 1
  ub_bdd[which(ub_bdd > 1)] = 1
  
  if(save_img)
    pdf(img_name_1, width = width, height = height )
  
  if(idx == 15){
    par( mfrow = c(1,1), mar = c(4,5,1,0.5), mgp=c(3.5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
    plot(0,0,  yaxt = "n", xaxt = "n",
         xlab = "", ylab = "1000 * Bound",
         xlim = c(0,max(n_train_vec) ) , ylim = c(ymin,ymax), 
         main = paste0(" "),
         type = "n")
  } else{
    par( mfrow = c(1,1), mar = c(4,3.5,1,0.5), mgp=c(2,1,0), bty = "l", las = 1, cex.lab = cex.lab )
    plot(0,0,  yaxt = "n", xaxt = "n",
         xlab = "", ylab = "",
         xlim = c(0,max(n_train_vec) ) , ylim = c(ymin,ymax), 
         main = paste0(" "),
         type = "n")
  }
  grid(lty = 1,lwd = 1, col = "gray90" )
  axis(side = 2, at = ylabs, 
       labels = ylabs, las = 1, 
       cex.axis = cex.axis )
  axis(1, cex.axis = cex.axis);mtext("n", side = 1, line = 2, cex = cex.axis)
  points( x = n_train_vec, y = ub_unb[2,], 
          type = "l", 
          lwd = 5, pch = 16, lty = 1,
          col = "darkorange" )
  points( x = n_train_vec, y = ub_bdd[2,], 
          type = "l", 
          lwd = 5, pch = 16, lty = 1,
          col = "darkgreen" ) 
  abline( h = c(0.05,0.01), col = "grey50", lty = c(2,3), lwd = 3)
  legend("topright",c("Unbounded","Bounded"), 
         lwd = 3, col = c("darkorange","darkgreen"), cex = 1.5)
  if(save_img)
    dev.off()
  
}


# Plot in log scale -------------------------------------------------------

save_img = TRUE
width = 8; height = 6
cex.lab = 2
cex.axis = 2

idx = 15
for(idx in idxs) {
  cancer_name = cancer_names[idx]
  filename = paste0(wd,"TCGA/",cancer_types[idx],"_targeted.RData")
  load(filename)
  cat("\n idx = ",idx,"; cancer name: ",cancer_name,"\n")
  
  seed = 42
  set.seed(seed)
  
  exp_name = paste0("Cancer_",cancer_types[idx])
  file_name = paste0("save/",exp_name,".Rdat")
  img_name_2 = paste0("img/",exp_name,"_log.pdf")
  
  Nj = c(apply(Z, 2, sum))
  n = nrow(Z)
  Kobs = length( which(Nj > 0) )
  TabNj = c(length(which(Nj == 0)), tabulate(Nj, nbins = n) )
  train_prop = seq(0.05,1,by = 0.05)
  Nexp = length(train_prop)
  Nrep_max = 200
  Nrep = ceiling(seq(Nrep_max,1,length.out = Nexp))
  alfa = 0.05
  beta = 1e-5
  
  load(file_name)  
  ub_unb_list = res$ub_unb_list
  ub_bdd_list = res$ub_bdd_list
  ub_unb  = lapply(ub_unb_list, quantile, probs = c(0.025,0.5,0.975)); ub_unb = do.call(cbind, ub_unb)
  ub_bdd  = lapply(ub_bdd_list, quantile, probs = c(0.025,0.5,0.975)); ub_bdd = do.call(cbind, ub_bdd)
  
  ymin = 1e-2; ymax = 1
  # ylabs = round(seq(ymin,ymax,length.out = 5),1)
  
  n_train_vec = ceiling( n*train_prop )
  
  
  ub_unb[which(ub_unb > 1)] = 1
  ub_bdd[which(ub_bdd > 1)] = 1
  
  if(save_img)
    pdf(img_name_2, width = width, height = height )
  
  if(idx == 15){
    par( mfrow = c(1,1), mar = c(3,6,1,0.5), mgp=c(4.5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
    plot(1,1,  yaxt = "n", xaxt = "n",
         xlab = "", ylab = "Bound",
         xlim = c(1,max(n_train_vec) ) , ylim = c(ymin,ymax), 
         main = paste0(" "),
         log = "y",
         type = "n")
  } else{
    par( mfrow = c(1,1), mar = c(3,4.5,1,0.5), mgp=c(4.5,1,0), bty = "l", las = 1, cex.lab = cex.lab )
    plot(1,1,  yaxt = "n", xaxt = "n",
         xlab = "", ylab = " ",
         xlim = c(1,max(n_train_vec) ) , ylim = c(ymin,ymax), 
         main = paste0(" "),
         log = "y",
         type = "n")
  }
  grid(lty = 1,lwd = 1, col = "gray90" )
  axis(2, las = 1, cex.axis = cex.axis)
  axis(1, cex.axis = cex.axis);mtext("n", side = 1, line = 2, cex = cex.axis)
  points( x = n_train_vec, y = ub_unb[2,], 
          type = "l", 
          lwd = 5, pch = 16, lty = 1,
          col = "darkorange" )
  points( x = n_train_vec, y = ub_bdd[2,], 
          type = "l", 
          lwd = 5, pch = 16, lty = 1,
          col = "darkgreen" ) 
  abline( h = c(0.05,0.01), col = "grey50", lty = c(2,3), lwd = 3)
  legend("topright",c("Unbounded","Bounded"), 
         lwd = 3, col = c("darkorange","darkgreen"), cex = 1.5)
  if(save_img)
    dev.off()
  
}

