# Generate Distributions (features) ----------------------------------------------------
sim_Uniform_features = function(M){
  w = runif(M)
  w
}

sim_Constant_features = function(M, c){
  w = rep(c,M)
  w
}

sim_TruncatedZipfs_features = function(M,s){
  w = sapply(2:(M+1),function(j) j^(-s))
  w 
}

sim_TruncatedGeom_features = function(M,a){
  w = sapply(2:(M+1),function(j) (1-a)^(j-1) )
  w 
}


sim_generic = function(name,M,param){
  if(name == "Zipfs"){
    sim_TruncatedZipfs_features(M,param)
  }else if(name == "Constant"){
    sim_Constant_features(M,param)
  }else if( name == "Geom" ){
    sim_TruncatedGeom_features(M,param)
  }else if( name == "Uniform"){
    sim_Uniform_features(M)
  }
  else
    stop("Invalid name")
    
}

# This function generates a zipfs distribution with parameter s.
# Since the distribution is unbounded, it is truncated at level Mstar
# where Mstar is so that of NOT observing any feature after the Mstar
# is larger that 1-eps
# P(\sum_{j=Mstar}^\infty N_j = 0) >= 1-eps
sim_zipfs_features = function(s,n,eps){
  if(s<=1)
    stop("s must be positive")
  logmstar = 1/(s-1) * ( log(n) - log(s-1) - log(-log(1-eps)) )
  mstar = ceiling(logmstar)
  
  w = sim_TruncatedZipfs_features(mstar,s)
  w
}







get_first3digits <- function(x,d=3) {
  # convert to string without scientific notation
  s <- format(x, scientific = FALSE, trim = TRUE)
  # remove decimal point
  s <- gsub("\\.", "", s)
  # take first 3 digits
  s <- substr(s, 1, d)
  return(s)
}



# Lorenzo Ghilotti's functions for extrapolation curves -------------------------------------
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





