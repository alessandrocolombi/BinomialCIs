# Generate Distributions (features) ----------------------------------------------------

sim_TruncatedZipfs_features = function(M,s){
  w = sapply(1:(M),function(j) j^(-s))
  w 
}



