# Generate Distributions (features) ----------------------------------------------------
sim_Uniform_features = function(M){
  w = runif(M)
}


sim_TruncatedZipfs_features = function(M,s){
  w = sapply(2:(M+1),function(j) j^(-s))
  w 
}

sim_TruncatedGeom_features = function(M,a){
  w = sapply(2:(M+1),function(j) (1-a)^(j-1) )
  w 
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















