
# Gen. Dist. (species) ----------------------------------------------------
sim_zipfs = function(M,s){
  w = sapply(1:M,function(j) j^(-s))
  w / sum(w)
}
sim_geom = function(M,a){
  w = sapply(0:(M-1),function(j) (1-a)^(j-1)  )
  w / sum(w)
}
sim_negbin = function(M,l,r){
  w = dnbinom(x = 0:(M-1), size = l, prob = r) # Painsky parametrization
  w / sum(w)
}
sim_betabin = function(M,a,b){
  w = dbetabinom.ab(x = 0:(M-1), size = M-1, shape1 = a, shape2 = b)
  # lw = sapply(1:M,function(j) lchoose(M,j) + lbeta(j+a,M-j+b) - lbeta(a,b) )
  # lw = lchoose(M,1:M) + lbeta(1:M + a, M - 1:M + b )  
  # max_lw = max(lw)
  # w = exp(lw - max_lw)*exp(max_lw)
  w / sum(w)
}
sim_unif = function(M){
  w = 1:M
  w / sum(w)
}

# Upper bound Painsky -----------------------------------------------------
ub_pain = function(n,Rmax,alfa){
  rgrid = 1:Rmax
  qstar_pan = function(n,r){
    (r-1)/(r-1+n)
  }
  pan_ub_grid = sapply(rgrid, function(r){
    qstar = qstar_pan(n,r)
    res = qstar^(r-1) * (1 - qstar )^n * (1/alfa)
    res^(1/r)
  })
  pan_ub = min(pan_ub_grid)
  pan_ub
}


# Worst case dist. --------------------------------------------------------
compute_r_unbounded = function(n,Rmax,alfa){
  rgrid = 1:Rmax
  qstar_pan = function(n,r){
    (r-1)/(r-1+n)
  }
  pan_ub_grid = sapply(rgrid, function(r){
    qstar = qstar_pan(n,r)
    res = qstar^(r-1) * (1 - qstar )^n * (1/alfa)
    res^(1/r)
  })
  which.min(pan_ub_grid)
}


# Objective functions (species) -----------------------------------------------------
llik_pyp <- function(x) {
  alpha <- x[1]
  theta <- x[2]
  -log_eppfPYP( n, Kn, data_obs, alpha, theta )
}
llik_FD <- function(x) {
  gamma <- x[1]
  Lambda <- x[2]
  -log_eppfFD(n,Kn,data_obs,gamma,Lambda,M_max)
}

# Objective functions (features) -----------------------------------------------------
llik_PP <- function(x) {
  alpha <- x[1]
  u <- x[2]
  c <- u - alpha
  -log_efpfBeBePois(n, Kn, data_obs, alpha, c, gamma )
}
llik_PP2 <- function(x) {
  alpha <- 0
  gamma <- x[1]
  c <- x[2]
  -log_efpfBeBePois(n, Kn, data_obs, alpha, c, gamma )
}
llik_PP3Parm <- function(x) {
  alpha <- x[1]
  gamma <- x[2]
  u     <- x[3]
  c     <- u - alpha
  -log_efpfBeBePois(n, Kn, data_obs, alpha, c, gamma )
}
llik_MixPois <- function(x) {
  alpha <- x[1]
  u <- x[2]
  c <- u - alpha
  mu_gamma <- x[3]
  -log_efpfBeBeMixPois( n, Kn, data_obs, alpha, c, mu_gamma, var_gamma )
}
llik_MixBin <- function(x) {
  a <- x[1]
  b <- x[2]
  mu_nb <- x[3]
  -log_efpfBeBeMixNBin( n, Kn, data_obs, a, b, mu_nb, var_nb )
}
f_beta <- function(x){
  if(n <= 0)
    stop("Error in lf_beta: n<=0")
  if(alfa-x <= 0)
    stop("Error in lf_beta: alfa-x<=0")
  if( 1/x <= 0 )
    stop("Error in lf_beta: 1/x<=0")
  n/(alfa-x) * (  sqrt( (log(1/x))/(2*n) ) + sqrt( (log(1/x))/(2*n) + Shat) )
}
lf_beta <- function(x){
  if(n <= 0)
    stop("Error in lf_beta: n<=0")
  if(alfa-x <= 0)
    stop("Error in lf_beta: alfa-x<=0")
  if( 1/x <= 0 )
    stop("Error in lf_beta: 1/x<=0")
  
  log(n) -log(alfa-x) + log(  sqrt( (log(1/x))/(2*n) ) + sqrt( (log(1/x))/(2*n) + Shat) )
}

# Gen. Dist. (features) ----------------------------------------------------
sim_zipfs_features = function(M,s){
  w = sapply(1:(M),function(j) j^(-s))
  w 
}
sim_geom_features = function(M,a){
  w = sapply(0:(M-1),function(j) (1-a)^(j-1)  )
}
sim_ghilo_features = function(M){
  m = floor(M/3)
  w = rep(0,M)
  w[1:m] = 0.015
  w[(m+1):(2*m)] = 0.01
  w[(2*m+1):M] = 0.005
  w 
}


# Utilities ---------------------------------------------------------------
gamma_moments = function(a,b){
  list("mean" = a/b, "var" = a/(b*b))
}
gamma_shape_rate = function(mu,sig2){
  list("shape" = (mu*mu)/sig2, "rate" = mu/(sig2))
}

NegBin_moments = function(r,p){
  list("mean" = r*(1-p)/(p), "var" = r*(1-p)/(p*p))
}
NegBin_params = function(mu,sig2){
  list("p" = (mu)/sig2, "r" = (mu*mu)/(sig2-mu))
}
