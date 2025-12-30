# Generate Distributions (features) ----------------------------------------------------
sim_Uniform_features = function(M, a){
  w = runif(M,min = 0, max = a)
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
  if(a >= 1)
    stop("a must be strictly less than 1")
  w = sapply(2:(M+1),function(j) (1-a)^(j-1) )
  w 
}

sim_TruncatedGamma_features = function(M,lambda){
  w = sapply(2:(M+1),function(j) exp(-lambda*j) )
  w 
}

sim_TruncatedUnif_features = function(M, Meff = 2000, pmax = 1){
  if(M < Meff)
    stop("Error in sim_TruncatedUnif_features. M must be >= than Meff")
  w = rep(0,M)
  w[1:Meff] = runif(n = Meff, min = 0, max = pmax)
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
    sim_Uniform_features(M, param)
  }
  else if( name == "TrUnif"){
    sim_TruncatedUnif_features(M, Meff = 2000, param)
  }
  else if( name == "gamma"){
    sim_TruncatedGamma_features(M,param)
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


# Generate pj controlling S

sim_CntS_Constant = function(S,M){
  if(S <= 0)
    stop("Error in sim_CntS_Constant. S must be positive")
  w = rep(S/M,M)
  w
}

# sim_CntS_Zipfs = function(S,s){
#   if(S <= 0)
#     stop("Error in sim_CntS_Constant. S must be positive")
#   Mmax = 1000000
#   cumS = 0
#   j = 1 
#   flag = TRUE
#   w = c()
#   while(flag) {
#     new = (j+1)^(-s)
#     w = c(w,new)
#     cumS = cumS + new
#     if(cumS >= S)
#       flag = FALSE
#     if( (j+1) > Mmax  )
#       flag = FALSE
#     
#     j = j + 1
#   }
#   if( sum(w) < S )
#     stop("The sum is smaller than S")
#   w
# }

sim_CntS_gamma = function(M,S,lambda){
  if(S <= 0)
    stop("Error in sim_CntS_Constant. S must be positive")
  w = sapply(1:(M),function(j) exp(-lambda*j)  )
  stop("DEVO GARANTIRE CHE w SIA IN (0,1)")
}


# utilities ---------------------------------------------------------------




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



# Stopping rule functions -------------------------------------------------

# Chao–Jost sample coverage for abundance data
# counts: vector of species abundances (only observed species)
coverage_ChaoJost <- function(n,counts) {
  if (n <= 1) return(NA_real_)
  counts <- counts[counts > 0]
  
  f1 <- sum(counts == 1L)
  f2 <- sum(counts == 2L)
  
  if (f1 == 0L) return(1)  # all species have count >= 2, coverage ~ 1
  
  # Chao & Jost (2012), incidence-based sample coverage
  C_hat <- 1 - (f1 / n) * ( ( (n - 1) * f1 ) / ( (n - 1) * f1 + 2 * f2 ) )
  C_hat
}


Nsample_Chao2009 <- function(n,counts,g_Chao2009){
  
  if (n <= 1) return(NA_real_)
  counts <- counts[counts > 0]
  
  f1 <- sum(counts == 1L) # num. sigleton
  f2 <- sum(counts == 2L) # num. doubleton
  
  Sobs <- length(counts)  # num. distinct observed
  Sest <- Sobs
  if(f2 > 0){
    Sest <- Sest + ( (1 - 1/n)*f1*f1 )/(2*f2) # num. distinct estimated
  } else {
    Sest <- Sest + ( f1*(f1 - 1) )/(2) # num. distinct estimated
  }
  
  if(g_Chao2009*Sest < Sobs)
    return( 0 )
  
  num = log( 1 - (n*2*f2*(g_Chao2009*Sest - Sobs))/( (n-1)*f1*f1 )  ) # numerator of Eq.(15) in Chao et. al. (2009)
  den = log( 1 - (2*f2)/( (n-1)*f1 + 2*f2 )  )                        # denominator of Eq.(15) in Chao et. al. (2009)
  num/den
}




