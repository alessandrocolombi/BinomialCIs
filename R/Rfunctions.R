# Generate Distributions (features) ----------------------------------------------------

sim_TruncatedZipfs_features = function(M,s){
  w = sapply(1:(M),function(j) j^(-s))
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
