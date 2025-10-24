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







