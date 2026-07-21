wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_bocconi = "/home/colombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100,wd_bocconi)
choose_wd = wd_vec[1] # <--- modify here
wd = paste0(choose_wd,"BinomialCIs/Rscripts/SecondSub")
setwd(wd)

# Librerie ----------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
suppressWarnings(suppressPackageStartupMessages(library(progress)))
suppressWarnings(suppressPackageStartupMessages(library(VGAM)))
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
source("../../R/Rfunctions.R")

# Sim specific param ----------------------------------------------------------------
RunParallel = FALSE
if(!RunParallel){
  idx = 1 
  param = 1
}

## ------------------------------------------------------------
## Packages
## ------------------------------------------------------------
library(future.apply)

## ------------------------------------------------------------
## 1. Utilities: communities and coverage estimator
## ------------------------------------------------------------

# Build p_j ∝ j^{-gamma} on {1, ..., M}
build_zipf_p <- function(M = 1500L, gamma = 0.75) {
  j <- seq_len(M) + 1
  p <- j^(-gamma)
  p
}

build_uniform_p <- function(M = 1500L, c=0.05) {
  p <- rep(c, M)
  p
}

build_geom_p <- function(M = 1500L, a = 0.1) {
  if (a >= 1)
    stop("a must be strictly less than 1")
  j <- 2:(M + 1L)
  # This matches your sim_TruncatedGeom_features:
  # w_j = (1 - a)^(j-1), j = 2,...,M+1  -> (1-a)^1,...,(1-a)^M
  w <- (1 - a)^(j - 1L)
  as.numeric(w)
}

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


## ============================================================================
## STOPPING RULES: DEFINIZIONE DEI CRITERI DI STOP
## ============================================================================
## Questa e' l'unica sezione da modificare quando cambia la disuguaglianza
## che determina lo stop. Il calcolo delle statistiche U_*, C_hat e mg resta
## nel ciclo di simulazione, subito prima della chiamata a queste funzioni.

STOPPING_RULE_NAMES <- c(
  "Bounded",
  "Unbounded",
  "Chao2009",
  "Coverage",
  "Bounded new",
  "Unbounded new"
)

# 1) Bounded: stop quando l'upper bound e' strettamente minore di eps.
stop_Bounded <- function(U_bounded, eps) {
  !is.na(U_bounded) && U_bounded < eps
}

# 2) Unbounded: stop quando l'upper bound e' minore o uguale a eps.
stop_Unbounded <- function(U_unbounded, eps) {
  !is.na(U_unbounded) && U_unbounded <= eps
}

# 3) Chao2009: stop quando il batch corrente copre il campione addizionale mg.
stop_Chao2009 <- function(mg, batch_size) {
  !is.na(mg) && batch_size >= mg
}

# 4) Coverage: stop quando la coverage stimata raggiunge il target.
stop_Coverage <- function(C_hat, C_target) {
  !is.na(C_hat) && C_hat >= C_target
}

# 5) Unbounded new.
stop_Unbounded_new <- function(log_Sn_bar, eps, n, alpha, beta) {
  lhs <- log_Sn_bar + n * log(1 - eps) - log(eps)
  rhs <- log(alpha - beta)
  !is.na(log_Sn_bar) && lhs < rhs
}

# 6) Bounded new.
stop_Bounded_new <- function(log_m_bar, eps, n, alpha, delta, b) {
  lhs <- log_m_bar + log(n - b) + log(eps) - log(1 - eps)
  rhs <- log(alpha - delta)
  !is.na(log_m_bar) && lhs < rhs
}


## ------------------------------------------------------------
## 2. One simulation run for a fixed community and (q, ε, C_target)
## ------------------------------------------------------------
# Assumes you have:
#   compute_UB_analytical(n, n_i, M, b_n, alpha, log_scale_flag)
#   compute_UB_rnorm(n, alpha, beta, r_n, Shat)
# already defined somewhere in your project.

simulate_one_run_old <- function(p,
                             q_error        = 0,
                             eps            = 0.01,
                             alpha          = 0.05,
                             C_target       = 0.95,
                             g_Chao2009     = 0.95,
                             batch_size     = 50L,
                             n_max          = 5000L,
                             beta           = 1e-5) {
  
  M <- length(p)
  # "Big" true species for correctness check: p_j >= eps
  big_species <- which(p >= eps)
  
  # Counters
  n        <- 0L            # number of sampling units (rows)
  n_error  <- 0L            # number of singleton error species observed so far
  counts_true <- integer(M) # counts of true species across all rows
  
  # Stopping flags and outputs
  stopped_bounded   <- FALSE
  stopped_unbounded <- FALSE
  stopped_cov       <- FALSE
  stopped_Chao2009  <- FALSE
  
  Nstop_bounded   <- NA_integer_
  Nstop_unbounded <- NA_integer_
  Nstop_cov       <- NA_integer_
  Nstop_Chao2009  <- NA_integer_
  
  ok_bounded   <- NA
  ok_unbounded <- NA
  ok_cov       <- NA
  ok_Chao2009  <- NA
  
  # New metrics: missed big species and extra species at stopping
  missed_big_bounded      <- NA_integer_
  missed_big_unbounded    <- NA_integer_
  missed_big_coverage     <- NA_integer_
  missed_big_Chao2009     <- NA_integer_
  extra_species_bounded   <- NA_integer_
  extra_species_unbounded <- NA_integer_
  extra_species_coverage  <- NA_integer_
  extra_species_Chao2009  <- NA_integer_
  
  n_batches <- ceiling(n_max / batch_size)
  
  for (b in seq_len(n_batches)) {
    # cat("\n n = ",n," \n")
    # Allow for a non-multiple n_max if needed
    remaining <- n_max - n
    if (remaining <= 0L) break
    this_b <- min(batch_size, remaining)
    
    ## ---- Bernoulli product: X_{i,j} ~ Bern(p_j), i = 1,...,this_b ----
    # For each species j, number of 1's in this batch: Binomial(this_b, p_j)
    successes <- stats::rbinom(n = M, size = this_b, prob = p)
    
    ## ---- Entry-level contamination ----
    # Each 1 can be an error with prob q_error:
    # if it's an error, we remove it from the true species and
    # create a new singleton error species.
    if (q_error > 0) {
      errors <- stats::rbinom(n = M, size = successes, prob = q_error)
    } else {
      errors <- integer(M)
    }
    
    counts_true <- counts_true + (successes - errors)
    n_error     <- n_error     + sum(errors)
    n           <- n + this_b
    
    ## ---- Observed abundance vector (true + error species) ----
    counts_obs <- c(counts_true[counts_true > 0L],
                    rep(1L, n_error))   # n_error singleton OTUs
    M_obs <- length(counts_obs)         # observed alphabet size
    
    if (length(counts_obs) == 0L) next   # nothing observed yet
    
    ## ---- 5.1 Bounded-alphabet CI rule (on observed data) ----
    if (!stopped_bounded) {
      b_n <- log(n)
      # Now Nj has length M_obs and we pass M_obs to the C++ routine
      nj <- c(counts_obs, rep(0L, 10 * M - M_obs))
      U_bounded <- compute_UB_analytical(n, nj, 10 * M, b_n, alpha, FALSE)
      
      if (stop_Bounded(U_bounded, eps)) {
        stopped_bounded <- TRUE
        Nstop_bounded   <- n
        # Correctness still checked on TRUE big species only
        ok_bounded      <- all(counts_true[big_species] > 0L)
        
        # Missed big species at stopping
        if (length(big_species) > 0L) {
          big_seen <- counts_true[big_species] > 0L
          missed_big_bounded <- sum(!big_seen)
        } else {
          missed_big_bounded <- 0L
        }
        # Extra species at stopping: total observed minus observed big species
        K_true_seen <- sum(counts_true > 0L)
        observed_big <- if (length(big_species) > 0L) sum(counts_true[big_species] > 0L) else 0L
        extra_species_bounded <- K_true_seen + n_error - observed_big
      }
    }
    
    ## ---- 5.2 Unbounded-alphabet CI rule (on observed data) ----
    if (!stopped_unbounded) {
      Shat  <- sum(counts_obs) / n
      Sstar <- ( sqrt( -log(beta) / (2 * n) ) +
                   sqrt( Shat + (-log(beta) / (2 * n)) ) )^2
      r_n   <- log( Sstar / (-log(1 - alpha + beta)) ) +
        log(n) - log(log(n))
      U_unbounded <- compute_UB_rnorm(n, alpha, beta, r_n, Shat)
      
      if (stop_Unbounded(U_unbounded, eps)) {
        stopped_unbounded <- TRUE
        Nstop_unbounded   <- n
        ok_unbounded      <- all(counts_true[big_species] > 0L)
        
        # Missed big species at stopping
        if (length(big_species) > 0L) {
          big_seen <- counts_true[big_species] > 0L
          missed_big_unbounded <- sum(!big_seen)
        } else {
          missed_big_unbounded <- 0L
        }
        # Extra species at stopping
        K_true_seen <- sum(counts_true > 0L)
        observed_big <- if (length(big_species) > 0L) sum(counts_true[big_species] > 0L) else 0L
        extra_species_unbounded <- K_true_seen + n_error - observed_big
      }
    }
    
    ## ---- 5.3 Chao–Jost coverage-based rule (same counts_obs) ----
    if (!stopped_cov) {
      C_hat <- coverage_ChaoJost(n,counts_obs)
      if (stop_Coverage(C_hat, C_target)) {
        # cat(" --> n = ",n," C_hat = ",C_hat," ---> fine Coverage \n")
        stopped_cov <- TRUE
        Nstop_cov   <- n
        ok_cov      <- all(counts_true[big_species] > 0L)
        
        # Missed big species at stopping
        if (length(big_species) > 0L) {
          big_seen <- counts_true[big_species] > 0L
          missed_big_coverage <- sum(!big_seen)
        } else {
          missed_big_coverage <- 0L
        }
        # Extra species at stopping
        K_true_seen <- sum(counts_true > 0L)
        observed_big <- if (length(big_species) > 0L) sum(counts_true[big_species] > 0L) else 0L
        extra_species_coverage <- K_true_seen + n_error - observed_big
      }
    }
    
    ## ---- 5.4 Chao2009 Eq.15 stopping rule (same counts_obs) ----
    if (!stopped_Chao2009) {
      mg <- Nsample_Chao2009(n,counts_obs,g_Chao2009)
      if (stop_Chao2009(mg, this_b)) {
        # cat(" --> n = ",n," mg = ",mg," ---> fine Chao2009 \n")
        stopped_Chao2009 <- TRUE
        Nstop_Chao2009   <- n + mg
        ok_Chao2009      <- all(counts_true[big_species] > 0L)
        
        # Missed big species at stopping
        if (length(big_species) > 0L) {
          big_seen <- counts_true[big_species] > 0L
          missed_big_Chao2009 <- sum(!big_seen)
        } else {
          missed_big_Chao2009 <- 0L
        }
        # Extra species at stopping
        K_true_seen <- sum(counts_true > 0L)
        observed_big <- if (length(big_species) > 0L) sum(counts_true[big_species] > 0L) else 0L
        extra_species_Chao2009 <- K_true_seen + n_error - observed_big
      }
    }
    
    # Early exit if all four rules have stopped
    if (stopped_bounded && stopped_unbounded && stopped_cov && stopped_Chao2009) break
  }
  cat("\n ################## \n")
  
  # If a rule never stopped by n_max, set Nstop_* = n_max
  # (missed/extra and ok_* remain NA to distinguish "never stopped").
  if (!stopped_bounded)   Nstop_bounded   <- n_max
  if (!stopped_unbounded) Nstop_unbounded <- n_max
  if (!stopped_cov)       Nstop_cov       <- n_max
  if (!stopped_Chao2009)  Nstop_Chao2009 <- n_max
  
  list(
    Nstop_bounded          = Nstop_bounded,
    Nstop_unbounded        = Nstop_unbounded,
    Nstop_coverage         = Nstop_cov,
    Nstop_Chao2009         = Nstop_Chao2009,
    ok_bounded             = ok_bounded,
    ok_unbounded           = ok_unbounded,
    ok_coverage            = ok_cov,
    ok_Chao2009            = ok_Chao2009,
    missed_big_bounded     = missed_big_bounded,
    missed_big_unbounded   = missed_big_unbounded,
    missed_big_coverage    = missed_big_coverage,
    missed_big_Chao2009    = missed_big_Chao2009,
    extra_species_bounded  = extra_species_bounded,
    extra_species_unbounded= extra_species_unbounded,
    extra_species_coverage = extra_species_coverage,
    extra_species_Chao2009 = extra_species_Chao2009
  )
}

simulate_one_run <- function(p,
                             q_error        = 0,
                             eps            = 0.01,
                             alpha          = 0.05,
                             C_target       = 0.95,
                             g_Chao2009     = 0.95,
                             batch_size     = 50L,
                             n_max          = 5000L,
                             beta           = 1e-5,
                             delta          = 1e-5,
                             b_new          = 1L) {
  
  M <- length(p)
  # "Big" true species for correctness check: p_j >= eps
  big_species <- which(p >= eps)
  
  # Counters
  n        <- 0L            # number of sampling units (rows)
  n_error  <- 0L            # number of singleton error species observed so far
  counts_true <- integer(M) # counts of true species across all rows
  
  # Stopping flags and outputs
  stopped_bounded   <- FALSE
  stopped_unbounded <- FALSE
  stopped_cov       <- FALSE
  stopped_Chao2009  <- FALSE
  stopped_bounded_new   <- FALSE
  stopped_unbounded_new <- FALSE
  
  Nstop_bounded   <- NA_integer_
  Nstop_unbounded <- NA_integer_
  Nstop_cov       <- NA_integer_
  Nstop_Chao2009  <- NA_integer_
  Nstop_bounded_new   <- NA_integer_
  Nstop_unbounded_new <- NA_integer_
  
  ok_bounded   <- NA
  ok_unbounded <- NA
  ok_cov       <- NA
  ok_Chao2009  <- NA
  ok_bounded_new   <- NA
  ok_unbounded_new <- NA
  
  # New metrics: missed big species and extra species at stopping
  missed_big_bounded      <- NA_integer_
  missed_big_unbounded    <- NA_integer_
  missed_big_coverage     <- NA_integer_
  missed_big_Chao2009     <- NA_integer_
  missed_big_bounded_new  <- NA_integer_
  missed_big_unbounded_new<- NA_integer_
  extra_species_bounded   <- NA_integer_
  extra_species_unbounded <- NA_integer_
  extra_species_coverage  <- NA_integer_
  extra_species_Chao2009  <- NA_integer_
  extra_species_bounded_new   <- NA_integer_
  extra_species_unbounded_new <- NA_integer_
  
  # Helper: compute (# missed big species, # extra species) given current counts
  compute_missed_and_extra <- function() {
    if (length(big_species) > 0L) {
      big_seen    <- counts_true[big_species] > 0L
      missed_big  <- sum(!big_seen)
      observed_big <- sum(big_seen)
    } else {
      missed_big   <- 0L
      observed_big <- 0L
    }
    K_true_seen   <- sum(counts_true > 0L)
    extra_species <- K_true_seen + n_error - observed_big
    list(missed_big = missed_big, extra_species = extra_species)
  }
  
  n_batches <- ceiling(n_max / batch_size)
  
  for (b in seq_len(n_batches)) {
    # Allow for a non-multiple n_max if needed
    remaining <- n_max - n
    if (remaining <= 0L) break
    this_b <- min(batch_size, remaining)
    
    ## ---- Bernoulli product: X_{i,j} ~ Bern(p_j), i = 1,...,this_b ----
    # For each species j, number of 1's in this batch: Binomial(this_b, p_j)
    successes <- stats::rbinom(n = M, size = this_b, prob = p)
    
    ## ---- Entry-level contamination ----
    # Each 1 can be an error with prob q_error:
    # if it's an error, we remove it from the true species and
    # create a new singleton error species.
    if (q_error > 0) {
      errors <- stats::rbinom(n = M, size = successes, prob = q_error)
    } else {
      errors <- integer(M)
    }
    
    counts_true <- counts_true + (successes - errors)
    n_error     <- n_error     + sum(errors)
    n           <- n + this_b
    
    ## ---- Observed abundance vector (true + error species) ----
    counts_obs <- c(counts_true[counts_true > 0L],
                    rep(1L, n_error))   # n_error singleton OTUs
    M_obs <- length(counts_obs)         # observed alphabet size
    
    if (length(counts_obs) == 0L) next   # nothing observed yet
    
    ## ---- 5.1 Bounded-alphabet CI rule (on observed data) ----
    if (!stopped_bounded) {
      b_n <- log(n)
      # Now Nj has length M_obs and we pass 10*M to the C++ routine
      nj <- c(counts_obs, rep(0L, 10 * M - M_obs))
      U_bounded <- compute_UB_analytical(n, nj, 10 * M, b_n, alpha, FALSE)
      
      if (stop_Bounded(U_bounded, eps)) {
        stopped_bounded <- TRUE
        Nstop_bounded   <- n
        ok_bounded      <- all(counts_true[big_species] > 0L)
        
        me <- compute_missed_and_extra()
        missed_big_bounded    <- me$missed_big
        extra_species_bounded <- me$extra_species
      }
    }
    
    ## ---- 5.2 Unbounded-alphabet CI rule (on observed data) ----
    if (!stopped_unbounded) {
      Shat  <- sum(counts_obs) / n
      Sstar <- ( sqrt( -log(beta) / (2 * n) ) +
                   sqrt( Shat + (-log(beta) / (2 * n)) ) )^2
      r_n   <- log( Sstar / (-log(1 - alpha + beta)) ) +
        log(n) - log(log(n))
      U_unbounded <- compute_UB_rnorm(n, alpha, beta, r_n, Shat)
      
      if (stop_Unbounded(U_unbounded, eps)) {
        stopped_unbounded <- TRUE
        Nstop_unbounded   <- n
        ok_unbounded      <- all(counts_true[big_species] > 0L)
        
        me <- compute_missed_and_extra()
        missed_big_unbounded    <- me$missed_big
        extra_species_unbounded <- me$extra_species
      }
    }
    
    ## ---- 5.3 Chao–Jost coverage-based rule (same counts_obs) ----
    if (!stopped_cov) {
      C_hat <- coverage_ChaoJost(n, counts_obs)
      if (stop_Coverage(C_hat, C_target)) {
        stopped_cov <- TRUE
        Nstop_cov   <- n
        ok_cov      <- all(counts_true[big_species] > 0L)
        
        me <- compute_missed_and_extra()
        missed_big_coverage    <- me$missed_big
        extra_species_coverage <- me$extra_species
      }
    }
    
    ## ---- 5.4 Chao2009 Eq.15 stopping rule (same counts_obs) ----
    if (!stopped_Chao2009) {
      mg <- Nsample_Chao2009(n, counts_obs, g_Chao2009)
      if (stop_Chao2009(mg, this_b)) {
        stopped_Chao2009 <- TRUE
        Nstop_Chao2009   <- n + mg
        ok_Chao2009      <- all(counts_true[big_species] > 0L)
        
        me <- compute_missed_and_extra()
        missed_big_Chao2009    <- me$missed_big
        extra_species_Chao2009 <- me$extra_species
      }
    }

    ## ---- 5.5 Bounded new (punto di inserimento) ----
    if (!stopped_bounded_new) {
      b_n <- log(n)
      # Now Nj has length M_obs and we pass 10*M to the C++ routine
      nj <- c(counts_obs, rep(0L, 10 * M - M_obs))
      alpha1 = 0.99*alpha
      delta = 0.01*alpha
      log_m_bar = compute_lmbar_bounded(n, nj, 10*M, b_n, alpha, n_max)

      if (stop_Bounded_new(log_m_bar, eps, n, alpha1, delta, b_n)) {
        stopped_bounded_new <- TRUE
        Nstop_bounded_new   <- n
        ok_bounded_new      <- all(counts_true[big_species] > 0L)

        me <- compute_missed_and_extra()
        missed_big_bounded_new    <- me$missed_big
        extra_species_bounded_new <- me$extra_species
      }
    }

    ## ---- 5.6 Unbounded new (punto di inserimento) ----
    if (!stopped_unbounded_new) {
      Shat  <- sum(counts_obs) / n
      Sstar <- ( sqrt( -log(beta) / (2 * n) ) + sqrt( Shat + (-log(beta) / (2 * n)) ) )^2
      r_n   <- log( Sstar / (-log(1 - alpha + beta)) ) + log(n) - log(log(n))
      
      Sn_bar <- ( sqrt( log(n_max)-log(beta) / (2 * n) ) + sqrt( Shat + (log(n_max)-log(beta) / (2 * n)) ) )^2
      Rn_bar <- log( Sn_bar / (-log(1 - alpha + beta)) ) + log(n) - log(log(n))
      
      if (stop_Unbounded_new(log(Sn_bar), eps, n, alpha, beta)) {
        stopped_unbounded_new <- TRUE
        Nstop_unbounded_new   <- n
        ok_unbounded_new      <- all(counts_true[big_species] > 0L)

        me <- compute_missed_and_extra()
        missed_big_unbounded_new    <- me$missed_big
        extra_species_unbounded_new <- me$extra_species
      }
    }
    
    # Early exit if all six rules have stopped
    if (stopped_bounded && stopped_unbounded && stopped_cov &&
        stopped_Chao2009 && stopped_bounded_new && stopped_unbounded_new) break
  }
  
  ## ------------------------------------------------------------
  ## Post-loop: handle rules that *never* stopped by n_max
  ## ------------------------------------------------------------
  # For those, set Nstop_* = n_max and compute missed/extra at n_max.
  # ok_* stays NA to keep your current type I error definition.
  me_final <- compute_missed_and_extra()
  
  if (!stopped_bounded) {
    Nstop_bounded          <- n_max
    missed_big_bounded     <- me_final$missed_big
    extra_species_bounded  <- me_final$extra_species
  }
  if (!stopped_unbounded) {
    Nstop_unbounded        <- n_max
    missed_big_unbounded   <- me_final$missed_big
    extra_species_unbounded<- me_final$extra_species
  }
  if (!stopped_cov) {
    Nstop_cov              <- n_max
    missed_big_coverage    <- me_final$missed_big
    extra_species_coverage <- me_final$extra_species
  }
  if (!stopped_Chao2009) {
    Nstop_Chao2009         <- n_max
    missed_big_Chao2009    <- me_final$missed_big
    extra_species_Chao2009 <- me_final$extra_species
  }
  if (!stopped_bounded_new) {
    Nstop_bounded_new         <- n_max
    missed_big_bounded_new    <- me_final$missed_big
    extra_species_bounded_new <- me_final$extra_species
  }
  if (!stopped_unbounded_new) {
    Nstop_unbounded_new         <- n_max
    missed_big_unbounded_new    <- me_final$missed_big
    extra_species_unbounded_new <- me_final$extra_species
  }
  
  list(
    Nstop_bounded           = Nstop_bounded,
    Nstop_unbounded         = Nstop_unbounded,
    Nstop_coverage          = Nstop_cov,
    Nstop_Chao2009          = Nstop_Chao2009,
    Nstop_bounded_new       = Nstop_bounded_new,
    Nstop_unbounded_new     = Nstop_unbounded_new,
    ok_bounded              = ok_bounded,
    ok_unbounded            = ok_unbounded,
    ok_coverage             = ok_cov,
    ok_Chao2009             = ok_Chao2009,
    ok_bounded_new          = ok_bounded_new,
    ok_unbounded_new        = ok_unbounded_new,
    missed_big_bounded      = missed_big_bounded,
    missed_big_unbounded    = missed_big_unbounded,
    missed_big_coverage     = missed_big_coverage,
    missed_big_Chao2009     = missed_big_Chao2009,
    missed_big_bounded_new  = missed_big_bounded_new,
    missed_big_unbounded_new= missed_big_unbounded_new,
    extra_species_bounded   = extra_species_bounded,
    extra_species_unbounded = extra_species_unbounded,
    extra_species_coverage  = extra_species_coverage,
    extra_species_Chao2009  = extra_species_Chao2009,
    extra_species_bounded_new   = extra_species_bounded_new,
    extra_species_unbounded_new = extra_species_unbounded_new
  )
}


## ------------------------------------------------------------
## 3'. SEQUENTIAL driver over many synthetic datasets
## ------------------------------------------------------------

# Converte un numero di secondi in una stringa leggibile nei log della VM.
format_elapsed_time <- function(seconds) {
  if (!is.finite(seconds) || seconds < 0) return("n/a")
  seconds <- as.integer(round(seconds))
  hours <- seconds %/% 3600L
  minutes <- (seconds %% 3600L) %/% 60L
  secs <- seconds %% 60L
  sprintf("%02d:%02d:%02d", hours, minutes, secs)
}

run_simulation_sequential <- function(distro         = c("zipf_heavy",
                                                         "unif_small",
                                                         "unif_big",
                                                         "geom_0.05"),
                                      M              = 1500L,
                                      q_vec          = c(0, 0.001, 0.005, 0.01, 0.05),
                                      eps            = 0.01,
                                      alpha          = 0.05,
                                      C_target       = 0.99,
                                      g_Chao2009     = 0.99,
                                      batch_size     = 50L,
                                      n_max          = 5000L,
                                      n_reps         = 200L,
                                      beta           = 1e-5,
                                      delta          = 1e-5,
                                      b_new          = 1L,
                                      progress_every = 50L,
                                      verbose        = TRUE,
                                      seed           = 123) {
  
  distro <- match.arg(distro)
  
  # Build the fixed "true" community p, depending on scenario
  p <- switch(
    distro,
    "zipf_regular"   = build_zipf_p(M = M, gamma = 0.7),
    "zipf_heavy"     = build_zipf_p(M = M, gamma = 1.05),
    "unif_small"     = build_uniform_p(M = M, c = 0.006),
    "unif_big"       = build_uniform_p(M = M, c = 0.05),
    "geom_0.1"       = build_geom_p(M = M, a = 0.1),
    "geom_0.05"      = build_geom_p(M = M, a = 0.05),
    stop("Unknown distro")
  )
  
  # Grid of (q, replicate)
  grid <- expand.grid(
    q_error = q_vec,
    rep_id  = seq_len(n_reps),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Reproducible seeds per row
  set.seed(seed)
  grid$sim_seed <- sample.int(.Machine$integer.max, nrow(grid))

  ## ------------------------------------------------------------------------
  ## Monitoraggio dell'esecuzione
  ## ------------------------------------------------------------------------
  ## Un job corrisponde a una coppia (livello di contaminazione q, replica).
  ## progress_every controlla ogni quanti job viene scritta una riga nel log.
  total_jobs <- nrow(grid)
  progress_every <- max(1L, as.integer(progress_every))
  start_time <- Sys.time()

  if (isTRUE(verbose)) {
    cat(sprintf(
      "[%s] START  scenario=%s | jobs=%d (%d q x %d repliche) | M=%d | n_max=%d\n",
      format(start_time, "%Y-%m-%d %H:%M:%S"),
      distro,
      total_jobs,
      length(q_vec),
      n_reps,
      M,
      n_max
    ))
    flush.console()
  }
  
  # Sequential loop (no future.apply)
  res_list <- lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    set.seed(g$sim_seed)
    
    out <- simulate_one_run(
      p          = p,
      q_error    = g$q_error,
      eps        = eps,
      alpha      = alpha,
      C_target   = C_target,
      g_Chao2009 = g_Chao2009,
      batch_size = batch_size,
      n_max      = n_max,
      beta       = beta,
      delta      = delta,
      b_new      = b_new
    )
    
    result_i <- c(
      list(
        distro  = distro,
        q_error = g$q_error,
        rep_id  = g$rep_id
      ),
      out
    )

    # Scrive un aggiornamento periodico e sempre alla conclusione dello scenario.
    if (isTRUE(verbose) && (i %% progress_every == 0L || i == total_jobs)) {
      now <- Sys.time()
      elapsed_seconds <- as.numeric(difftime(now, start_time, units = "secs"))
      jobs_per_second <- i / max(elapsed_seconds, .Machine$double.eps)
      eta_seconds <- (total_jobs - i) / jobs_per_second

      cat(sprintf(
        paste0(
          "[%s] PROGRESS scenario=%s | %d/%d (%.1f%%) | ",
          "q=%g rep=%d | elapsed=%s | ETA=%s | rate=%.3f job/s\n"
        ),
        format(now, "%Y-%m-%d %H:%M:%S"),
        distro,
        i,
        total_jobs,
        100 * i / total_jobs,
        g$q_error,
        g$rep_id,
        format_elapsed_time(elapsed_seconds),
        format_elapsed_time(eta_seconds),
        jobs_per_second
      ))
      flush.console()
    }

    result_i
  })

  if (isTRUE(verbose)) {
    end_time <- Sys.time()
    total_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
    cat(sprintf(
      "[%s] DONE   scenario=%s | jobs=%d | total elapsed=%s\n",
      format(end_time, "%Y-%m-%d %H:%M:%S"),
      distro,
      total_jobs,
      format_elapsed_time(total_seconds)
    ))
    flush.console()
  }
  
  # Back to a data.frame
  res_df <- do.call(rbind, lapply(res_list, as.data.frame))
  rownames(res_df) <- NULL
  res_df
}

eps <- 0.005

hyperparams <- list(
  distros      = c("zipf_heavy",
                   "unif_small",
                   "unif_big",
                   "geom_0.005"),
  M            = 1500L,
  q_vec        = c(0, 0.0001, 0.0005, 0.001, 0.0025, 0.005),
  eps          = eps,
  alpha        = 0.05,
  batch_size   = 50L,
  n_max        = 10000L,
  C_target     = 0.99,
  g_Chao2009   = 0.99,
  n_reps       = 500L,
  beta         = 1e-5,
  delta        = 1e-5,
  b_new        = 1L,
  # Una riga di avanzamento ogni 50 job; ridurre per log piu' frequenti.
  progress_every = 50L,
  verbose        = TRUE,
  seed_regular = 123L,
  seed_heavy   = 456L  # you can reuse or add more seeds if you want
)

# Riepilogo iniziale dell'intero esperimento. Questo messaggio permette di
# verificare subito dal log della VM che i parametri caricati siano corretti.
experiment_start_time <- Sys.time()
cat(sprintf(
  paste0(
    "[%s] EXPERIMENT START | scenarios=%d | q-levels=%d | ",
    "repliche/scenario/q=%d | total jobs=%d\n"
  ),
  format(experiment_start_time, "%Y-%m-%d %H:%M:%S"),
  length(hyperparams$distros),
  length(hyperparams$q_vec),
  hyperparams$n_reps,
  length(hyperparams$distros) * length(hyperparams$q_vec) * hyperparams$n_reps
))
flush.console()


# Uniform big
res_unif_big <- run_simulation_sequential(
  distro     = "unif_big",
  M          = hyperparams$M,
  q_vec      = hyperparams$q_vec,
  eps        = hyperparams$eps,
  alpha      = hyperparams$alpha,
  C_target   = hyperparams$C_target,
  g_Chao2009 = hyperparams$g_Chao2009,
  batch_size = hyperparams$batch_size,
  n_max      = hyperparams$n_max,
  n_reps     = hyperparams$n_reps,
  beta       = hyperparams$beta,
  delta      = hyperparams$delta,
  b_new      = hyperparams$b_new,
  progress_every = hyperparams$progress_every,
  verbose    = hyperparams$verbose,
  seed       = 789L
)


# Zipf heavy
res_heavy <- run_simulation_sequential(
  distro     = "zipf_heavy",
  M          = hyperparams$M,
  q_vec      = hyperparams$q_vec,
  eps        = hyperparams$eps,
  alpha      = hyperparams$alpha,
  C_target   = hyperparams$C_target,
  g_Chao2009 = hyperparams$g_Chao2009,
  batch_size = hyperparams$batch_size,
  n_max      = hyperparams$n_max,
  n_reps     = hyperparams$n_reps,
  beta       = hyperparams$beta,
  delta      = hyperparams$delta,
  b_new      = hyperparams$b_new,
  progress_every = hyperparams$progress_every,
  verbose    = hyperparams$verbose,
  seed       = hyperparams$seed_heavy
)

# Uniform small
res_unif_small <- run_simulation_sequential(
  distro     = "unif_small",
  M          = hyperparams$M,
  q_vec      = hyperparams$q_vec,
  eps        = hyperparams$eps,
  alpha      = hyperparams$alpha,
  C_target   = hyperparams$C_target,
  g_Chao2009 = hyperparams$g_Chao2009,
  batch_size = hyperparams$batch_size,
  n_max      = hyperparams$n_max,
  n_reps     = hyperparams$n_reps,
  beta       = hyperparams$beta,
  delta      = hyperparams$delta,
  b_new      = hyperparams$b_new,
  progress_every = hyperparams$progress_every,
  verbose    = hyperparams$verbose,
  seed       = 789L
)





# Geometric a = 0.05
res_geom_005 <- run_simulation_sequential(
  distro     = "geom_0.05",
  M          = hyperparams$M,
  q_vec      = hyperparams$q_vec,
  eps        = hyperparams$eps,
  alpha      = hyperparams$alpha,
  C_target   = hyperparams$C_target,
  g_Chao2009 = hyperparams$g_Chao2009,
  batch_size = hyperparams$batch_size,
  n_max      = hyperparams$n_max,
  n_reps     = hyperparams$n_reps,
  beta       = hyperparams$beta,
  delta      = hyperparams$delta,
  b_new      = hyperparams$b_new,
  progress_every = hyperparams$progress_every,
  verbose    = hyperparams$verbose,
  seed       = 792L
)

# Combine all scenarios
res_all <- rbind(
  res_heavy,
  res_unif_small,
  res_unif_big,
  res_geom_005
)

# Da questo punto le simulazioni sono concluse; iniziano aggregazioni, tabelle
# e grafici, che normalmente richiedono molto meno tempo.
simulation_end_time <- Sys.time()
cat(sprintf(
  "[%s] ALL SIMULATIONS DONE | rows=%d | elapsed=%s | starting summaries\n",
  format(simulation_end_time, "%Y-%m-%d %H:%M:%S"),
  nrow(res_all),
  format_elapsed_time(as.numeric(difftime(
    simulation_end_time, experiment_start_time, units = "secs"
  )))
))
flush.console()

## ------------------------------------------------------------
## Aggregation helpers
## ------------------------------------------------------------

agg_mean <- function(x, name) {
  out <- aggregate(x,
                   by = list(distro = res_all$distro,
                             q_error = res_all$q_error),
                   FUN = mean,
                   na.rm = TRUE)
  names(out)[3] <- name
  out
}

agg_median <- function(x, name) {
  out <- aggregate(x,
                   by = list(distro = res_all$distro,
                             q_error = res_all$q_error),
                   FUN = median,
                   na.rm = TRUE)
  names(out)[3] <- name
  out
}

agg_sd <- function(x, name) {
  out <- aggregate(x,
                   by = list(distro = res_all$distro,
                             q_error = res_all$q_error),
                   FUN = sd,
                   na.rm = TRUE)
  names(out)[3] <- name
  out
}

## ------------------------------------------------------------
## Type I errors (prob of stopping while missing some big species)
## ------------------------------------------------------------

type1_bounded   <- agg_mean(!res_all$ok_bounded,   "type1_bounded")
type1_unbounded <- agg_mean(!res_all$ok_unbounded, "type1_unbounded")
type1_coverage  <- agg_mean(!res_all$ok_coverage,  "type1_coverage")
type1_Chao2009  <- agg_mean(!res_all$ok_Chao2009,  "type1_Chao2009")
type1_bounded_new   <- agg_mean(!res_all$ok_bounded_new,   "type1_bounded_new")
type1_unbounded_new <- agg_mean(!res_all$ok_unbounded_new, "type1_unbounded_new")

## ------------------------------------------------------------
## Stopping probabilities
## ------------------------------------------------------------

prop_stop_bounded   <- agg_mean(!is.na(res_all$ok_bounded),   "prop_stopped_bounded")
prop_stop_unbounded <- agg_mean(!is.na(res_all$ok_unbounded), "prop_stopped_unbounded")
prop_stop_coverage  <- agg_mean(!is.na(res_all$ok_coverage),  "prop_stopped_coverage")
prop_stop_Chao2009  <- agg_mean(!is.na(res_all$ok_Chao2009),  "prop_stopped_Chao2009")
prop_stop_bounded_new   <- agg_mean(!is.na(res_all$ok_bounded_new),   "prop_stopped_bounded_new")
prop_stop_unbounded_new <- agg_mean(!is.na(res_all$ok_unbounded_new), "prop_stopped_unbounded_new")

## ------------------------------------------------------------
## N_stop: mean, median, sd
## ------------------------------------------------------------

mean_Nstop_bounded    <- agg_mean(res_all$Nstop_bounded,    "mean_Nstop_bounded")
mean_Nstop_unbounded  <- agg_mean(res_all$Nstop_unbounded,  "mean_Nstop_unbounded")
mean_Nstop_coverage   <- agg_mean(res_all$Nstop_coverage,   "mean_Nstop_coverage")
mean_Nstop_Chao2009   <- agg_mean(res_all$Nstop_Chao2009,   "mean_Nstop_Chao2009")
mean_Nstop_bounded_new   <- agg_mean(res_all$Nstop_bounded_new,   "mean_Nstop_bounded_new")
mean_Nstop_unbounded_new <- agg_mean(res_all$Nstop_unbounded_new, "mean_Nstop_unbounded_new")

median_Nstop_bounded   <- agg_median(res_all$Nstop_bounded,   "median_Nstop_bounded")
median_Nstop_unbounded <- agg_median(res_all$Nstop_unbounded, "median_Nstop_unbounded")
median_Nstop_coverage  <- agg_median(res_all$Nstop_coverage,  "median_Nstop_coverage")
median_Nstop_Chao2009  <- agg_median(res_all$Nstop_Chao2009,  "median_Nstop_Chao2009")
median_Nstop_bounded_new   <- agg_median(res_all$Nstop_bounded_new,   "median_Nstop_bounded_new")
median_Nstop_unbounded_new <- agg_median(res_all$Nstop_unbounded_new, "median_Nstop_unbounded_new")

sd_Nstop_bounded    <- agg_sd(res_all$Nstop_bounded,    "sd_Nstop_bounded")
sd_Nstop_unbounded  <- agg_sd(res_all$Nstop_unbounded,  "sd_Nstop_unbounded")
sd_Nstop_coverage   <- agg_sd(res_all$Nstop_coverage,   "sd_Nstop_coverage")
sd_Nstop_Chao2009   <- agg_sd(res_all$Nstop_Chao2009,   "sd_Nstop_Chao2009")
sd_Nstop_bounded_new   <- agg_sd(res_all$Nstop_bounded_new,   "sd_Nstop_bounded_new")
sd_Nstop_unbounded_new <- agg_sd(res_all$Nstop_unbounded_new, "sd_Nstop_unbounded_new")

## ------------------------------------------------------------
## Missed big species and extra species: mean / sd
## ------------------------------------------------------------

# Missed big species
mean_missed_bounded    <- agg_mean(res_all$missed_big_bounded,    "mean_missed_bounded")
mean_missed_unbounded  <- agg_mean(res_all$missed_big_unbounded,  "mean_missed_unbounded")
mean_missed_coverage   <- agg_mean(res_all$missed_big_coverage,   "mean_missed_coverage")
mean_missed_Chao2009   <- agg_mean(res_all$missed_big_Chao2009,   "mean_missed_Chao2009")
mean_missed_bounded_new   <- agg_mean(res_all$missed_big_bounded_new,   "mean_missed_bounded_new")
mean_missed_unbounded_new <- agg_mean(res_all$missed_big_unbounded_new, "mean_missed_unbounded_new")

sd_missed_bounded      <- agg_sd(res_all$missed_big_bounded,      "sd_missed_bounded")
sd_missed_unbounded    <- agg_sd(res_all$missed_big_unbounded,    "sd_missed_unbounded")
sd_missed_coverage     <- agg_sd(res_all$missed_big_coverage,     "sd_missed_coverage")
sd_missed_Chao2009     <- agg_sd(res_all$missed_big_Chao2009,     "sd_missed_Chao2009")
sd_missed_bounded_new   <- agg_sd(res_all$missed_big_bounded_new,   "sd_missed_bounded_new")
sd_missed_unbounded_new <- agg_sd(res_all$missed_big_unbounded_new, "sd_missed_unbounded_new")

# Extra species
mean_extra_bounded     <- agg_mean(res_all$extra_species_bounded,   "mean_extra_bounded")
mean_extra_unbounded   <- agg_mean(res_all$extra_species_unbounded, "mean_extra_unbounded")
mean_extra_coverage    <- agg_mean(res_all$extra_species_coverage,  "mean_extra_coverage")
mean_extra_Chao2009    <- agg_mean(res_all$extra_species_Chao2009,  "mean_extra_Chao2009")
mean_extra_bounded_new   <- agg_mean(res_all$extra_species_bounded_new,   "mean_extra_bounded_new")
mean_extra_unbounded_new <- agg_mean(res_all$extra_species_unbounded_new, "mean_extra_unbounded_new")

sd_extra_bounded       <- agg_sd(res_all$extra_species_bounded,     "sd_extra_bounded")
sd_extra_unbounded     <- agg_sd(res_all$extra_species_unbounded,   "sd_extra_unbounded")
sd_extra_coverage      <- agg_sd(res_all$extra_species_coverage,    "sd_extra_coverage")
sd_extra_Chao2009      <- agg_sd(res_all$extra_species_Chao2009,    "sd_extra_Chao2009")
sd_extra_bounded_new   <- agg_sd(res_all$extra_species_bounded_new,   "sd_extra_bounded_new")
sd_extra_unbounded_new <- agg_sd(res_all$extra_species_unbounded_new, "sd_extra_unbounded_new")

## ------------------------------------------------------------
## Combine everything into one data frame
## ------------------------------------------------------------

summary_df <- Reduce(
  function(x, y) merge(x, y, by = c("distro", "q_error"), sort = TRUE),
  list(
    type1_bounded,
    type1_unbounded,
    type1_coverage,
    type1_Chao2009,
    type1_bounded_new,
    type1_unbounded_new,
    prop_stop_bounded,
    prop_stop_unbounded,
    prop_stop_coverage,
    prop_stop_Chao2009,
    prop_stop_bounded_new,
    prop_stop_unbounded_new,
    mean_Nstop_bounded,
    mean_Nstop_unbounded,
    mean_Nstop_coverage,
    mean_Nstop_Chao2009,
    mean_Nstop_bounded_new,
    mean_Nstop_unbounded_new,
    median_Nstop_bounded,
    median_Nstop_unbounded,
    median_Nstop_coverage,
    median_Nstop_Chao2009,
    median_Nstop_bounded_new,
    median_Nstop_unbounded_new,
    sd_Nstop_bounded,
    sd_Nstop_unbounded,
    sd_Nstop_coverage,
    sd_Nstop_Chao2009,
    sd_Nstop_bounded_new,
    sd_Nstop_unbounded_new,
    mean_missed_bounded,
    mean_missed_unbounded,
    mean_missed_coverage,
    mean_missed_Chao2009,
    mean_missed_bounded_new,
    mean_missed_unbounded_new,
    sd_missed_bounded,
    sd_missed_unbounded,
    sd_missed_coverage,
    sd_missed_Chao2009,
    sd_missed_bounded_new,
    sd_missed_unbounded_new,
    mean_extra_bounded,
    mean_extra_unbounded,
    mean_extra_coverage,
    mean_extra_Chao2009,
    mean_extra_bounded_new,
    mean_extra_unbounded_new,
    sd_extra_bounded,
    sd_extra_unbounded,
    sd_extra_coverage,
    sd_extra_Chao2009,
    sd_extra_bounded_new,
    sd_extra_unbounded_new
  )
)

# Inspect
summary_df

## ------------------------------------------------------------
## Save results
## ------------------------------------------------------------

tail_label <- paste(hyperparams$distros, collapse = "-")
q_label    <- paste(hyperparams$q_vec, collapse = "_")

eps_str   <- gsub("\\.", "p", as.character(hyperparams$eps))
alpha_str <- gsub("\\.", "p", as.character(hyperparams$alpha))
C_str     <- gsub("\\.", "p", as.character(hyperparams$C_target))
g_str     <- gsub("\\.", "p", as.character(hyperparams$g_Chao2009))

file_name <- sprintf(
  "sim_tail-%s_M-%d_eps-%s_alpha-%s_C-%s_batch-%d_nmax-%d_reps-%d_q-%s.RData",
  tail_label,
  hyperparams$M,
  eps_str,
  alpha_str,
  C_str,
  hyperparams$batch_size,
  hyperparams$n_max,
  hyperparams$n_reps,
  q_label
)

# Collect everything in a single object
sim_out <- list(
  hyperparams   = hyperparams,
  res_heavy     = res_heavy,
  res_unif_small= res_unif_small,
  res_geom_025  = res_geom_005,
  res_all       = res_all,
  summary_df    = summary_df
)

# Save to disk
save(sim_out, file = file_name)
cat(sprintf(
  "[%s] RESULTS SAVED | file=%s\n",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  file_name
))
flush.console()

summary_df

## ------------------------------------------------------------
## Tables in LaTeX: type I error and N_stop (mean ± sd)
## ------------------------------------------------------------

library(knitr)
library(kableExtra)

q_levels <- sort(unique(summary_df$q_error))
q_cols   <- paste0("q_", gsub("\\.", "p", as.character(q_levels)))
scenarios <- sort(unique(summary_df$distro))

## 1) Type I error table: rows = (scenario, method), cols = q

rows_type1 <- list()

for (tt in scenarios) {
  df_tt <- subset(summary_df, distro == tt)
  df_tt <- df_tt[match(q_levels, df_tt$q_error), ]
  
  # bounded
  row_b <- as.data.frame(t(df_tt$type1_bounded))
  names(row_b) <- q_cols
  row_b$scenario <- tt
  row_b$method   <- "bounded"
  row_b <- row_b[, c("scenario", "method", q_cols)]
  
  # unbounded
  row_u <- as.data.frame(t(df_tt$type1_unbounded))
  names(row_u) <- q_cols
  row_u$scenario <- tt
  row_u$method   <- "unbounded"
  row_u <- row_u[, c("scenario", "method", q_cols)]
  
  # coverage
  row_c <- as.data.frame(t(df_tt$type1_coverage))
  names(row_c) <- q_cols
  row_c$scenario <- tt
  row_c$method   <- "coverage"
  row_c <- row_c[, c("scenario", "method", q_cols)]
  
  # Chao 2009
  row_d <- as.data.frame(t(df_tt$type1_coverage))
  names(row_d) <- q_cols
  row_d$scenario <- tt
  row_d$method   <- "Chao2009"
  row_d <- row_d[, c("scenario", "method", q_cols)]

  # bounded new
  row_bn <- as.data.frame(t(df_tt$type1_bounded_new))
  names(row_bn) <- q_cols
  row_bn$scenario <- tt
  row_bn$method   <- "bounded new"
  row_bn <- row_bn[, c("scenario", "method", q_cols)]

  # unbounded new
  row_un <- as.data.frame(t(df_tt$type1_unbounded_new))
  names(row_un) <- q_cols
  row_un$scenario <- tt
  row_un$method   <- "unbounded new"
  row_un <- row_un[, c("scenario", "method", q_cols)]
  
  rows_type1[[length(rows_type1) + 1L]] <- row_b
  rows_type1[[length(rows_type1) + 1L]] <- row_u
  rows_type1[[length(rows_type1) + 1L]] <- row_c
  rows_type1[[length(rows_type1) + 1L]] <- row_d
  rows_type1[[length(rows_type1) + 1L]] <- row_bn
  rows_type1[[length(rows_type1) + 1L]] <- row_un
}

type1_table <- do.call(rbind, rows_type1)
type1_table[ , q_cols] <- lapply(type1_table[ , q_cols], function(x) round(x, 3))

type1_latex <- kable(
  type1_table,
  format  = "latex",
  booktabs = TRUE,
  caption = "Estimated type I error by scenario, method, and contamination level $q$",
  label   = "tab:type1"
) %>%
  kable_styling(full_width = FALSE, position = "center", latex_options = c("hold_position"))

cat(type1_latex)

## 2) Mean N_stop table (with sd and '-' if type I error > 2 alpha)

alpha <- hyperparams$alpha
rows_Nstop <- list()

for (tt in scenarios) {
  df_tt <- subset(summary_df, distro == tt)
  df_tt <- df_tt[match(q_levels, df_tt$q_error), ]
  q_names <- paste0("q_", gsub("\\.", "p", as.character(df_tt$q_error)))
  
  ## ---- bounded ----
  row_b_mean <- as.data.frame(t(df_tt$mean_Nstop_bounded))
  row_b_sd   <- as.data.frame(t(df_tt$sd_Nstop_bounded))
  row_b_chr <- as.data.frame(
    setNames(
      lapply(seq_along(q_levels), function(j) {
        ns <- row_b_mean[1, j]
        ss <- row_b_sd[1, j]
        t1 <- df_tt$type1_bounded[j]
        if (is.na(ns)) {
          NA_character_
        } else if (!is.na(t1) && t1 > 2 * alpha) {
          "-"
        } else {
          sprintf("%.1f (%.1f)", ns, ss)
        }
      }),
      q_names
    )
  )
  row_b_chr$scenario <- tt
  row_b_chr$method   <- "bounded"
  q_cols_this <- setdiff(colnames(row_b_chr), c("scenario", "method"))
  row_b_chr   <- row_b_chr[, c("scenario", "method", q_cols_this)]
  
  ## ---- unbounded ----
  row_u_mean <- as.data.frame(t(df_tt$mean_Nstop_unbounded))
  row_u_sd   <- as.data.frame(t(df_tt$sd_Nstop_unbounded))
  row_u_chr <- as.data.frame(
    setNames(
      lapply(seq_along(q_levels), function(j) {
        ns <- row_u_mean[1, j]
        ss <- row_u_sd[1, j]
        t1 <- df_tt$type1_unbounded[j]
        if (is.na(ns)) {
          NA_character_
        } else if (!is.na(t1) && t1 > 2 * alpha) {
          "-"
        } else {
          sprintf("%.1f (%.1f)", ns, ss)
        }
      }),
      q_names
    )
  )
  row_u_chr$scenario <- tt
  row_u_chr$method   <- "unbounded"
  q_cols_this <- setdiff(colnames(row_u_chr), c("scenario", "method"))
  row_u_chr   <- row_u_chr[, c("scenario", "method", q_cols_this)]
  
  ## ---- coverage ----
  row_c_mean <- as.data.frame(t(df_tt$mean_Nstop_coverage))
  row_c_sd   <- as.data.frame(t(df_tt$sd_Nstop_coverage))
  row_c_chr <- as.data.frame(
    setNames(
      lapply(seq_along(q_levels), function(j) {
        ns <- row_c_mean[1, j]
        ss <- row_c_sd[1, j]
        t1 <- df_tt$type1_coverage[j]
        if (is.na(ns)) {
          NA_character_
        } else {
          sprintf("%.1f (%.1f)", ns, ss)
        }
      }),
      q_names
    )
  )
  row_c_chr$scenario <- tt
  row_c_chr$method   <- "coverage"
  q_cols_this <- setdiff(colnames(row_c_chr), c("scenario", "method"))
  row_c_chr   <- row_c_chr[, c("scenario", "method", q_cols_this)]
  
  ## ---- Chao2009 ----
  row_d_mean <- as.data.frame(t(df_tt$mean_Nstop_Chao2009))
  row_d_sd   <- as.data.frame(t(df_tt$sd_Nstop_Chao2009))
  row_d_chr <- as.data.frame(
    setNames(
      lapply(seq_along(q_levels), function(j) {
        ns <- row_d_mean[1, j]
        ss <- row_d_sd[1, j]
        t1 <- df_tt$type1_Chao2009[j]
        if (is.na(ns)) {
          NA_character_
        } else {
          sprintf("%.1f (%.1f)", ns, ss)
        }
      }),
      q_names
    )
  )
  row_d_chr$scenario <- tt
  row_d_chr$method   <- "Chao2009"
  q_cols_this <- setdiff(colnames(row_c_chr), c("scenario", "method"))
  row_d_chr   <- row_d_chr[, c("scenario", "method", q_cols_this)]

  ## ---- Bounded new / Unbounded new ----
  make_new_Nstop_row <- function(mean_values, sd_values, type1_values, method) {
    out <- as.data.frame(setNames(
      lapply(seq_along(q_levels), function(j) {
        ns <- mean_values[j]
        ss <- sd_values[j]
        t1 <- type1_values[j]
        if (is.na(ns)) NA_character_
        else if (!is.na(t1) && t1 > 2 * alpha) "-"
        else sprintf("%.1f (%.1f)", ns, ss)
      }),
      q_names
    ))
    out$scenario <- tt
    out$method <- method
    out[, c("scenario", "method", q_names)]
  }

  row_bn_chr <- make_new_Nstop_row(
    df_tt$mean_Nstop_bounded_new,
    df_tt$sd_Nstop_bounded_new,
    df_tt$type1_bounded_new,
    "bounded new"
  )
  row_un_chr <- make_new_Nstop_row(
    df_tt$mean_Nstop_unbounded_new,
    df_tt$sd_Nstop_unbounded_new,
    df_tt$type1_unbounded_new,
    "unbounded new"
  )
  
  rows_Nstop[[length(rows_Nstop) + 1L]] <- row_b_chr
  rows_Nstop[[length(rows_Nstop) + 1L]] <- row_u_chr
  rows_Nstop[[length(rows_Nstop) + 1L]] <- row_c_chr
  rows_Nstop[[length(rows_Nstop) + 1L]] <- row_d_chr
  rows_Nstop[[length(rows_Nstop) + 1L]] <- row_bn_chr
  rows_Nstop[[length(rows_Nstop) + 1L]] <- row_un_chr
}

Nstop_table <- do.call(rbind, rows_Nstop)

Nstop_latex <- kable(
  Nstop_table,
  format  = "latex",
  booktabs = TRUE,
  caption = "Mean stopping sample size (standard deviation in parentheses) by scenario, method, and contamination level $q$",
  label   = "tab:nstop"
) %>%
  kable_styling(full_width = FALSE, position = "center", latex_options = c("hold_position"))

cat(Nstop_latex)



## ------------------------------------------------------------
## 3) Missed big species *proportion* table (mean (sd))
## ------------------------------------------------------------

# Number of "big" species per scenario: n_big(d) = sum_j 1{p_j >= eps}
distro_levels <- sort(unique(summary_df$distro))

n_big_vec <- setNames(
  sapply(distro_levels, function(d) {
    p_d <- switch(
      d,
      "zipf_heavy" = build_zipf_p(M = hyperparams$M, gamma = 1.05),
      "unif_big" = build_uniform_p(M = hyperparams$M, c = 0.05),
      "unif_small" = build_uniform_p(M = hyperparams$M, c = 0.006),
      "geom_0.05"  = build_geom_p(M = hyperparams$M, a = 0.05),
      stop(sprintf("Unknown distro '%s' when computing n_big", d))
    )
    sum(p_d >= hyperparams$eps)
  }),
  distro_levels
)

rows_missed <- list()

for (tt in scenarios) {
  df_tt <- subset(summary_df, distro == tt)
  df_tt <- df_tt[match(q_levels, df_tt$q_error), ]
  q_names <- paste0("q_", gsub("\\.", "p", as.character(df_tt$q_error)))
  
  n_big <- n_big_vec[[tt]]
  
  ## ---- bounded ----
  # convert mean / sd counts -> proportions
  row_b_mean <- as.data.frame(t(df_tt$mean_missed_bounded / n_big))
  row_b_sd   <- as.data.frame(t(df_tt$sd_missed_bounded   / n_big))
  row_b_chr <- as.data.frame(
    setNames(
      lapply(seq_along(q_levels), function(j) {
        mm <- row_b_mean[1, j]
        ss <- row_b_sd[1, j]
        if (is.na(mm)) {
          NA_character_
        } else {
          sprintf("%.3f (%.3f)", mm, ss)
        }
      }),
      q_names
    )
  )
  row_b_chr$scenario <- tt
  row_b_chr$method   <- "bounded"
  q_cols_this <- setdiff(colnames(row_b_chr), c("scenario", "method"))
  row_b_chr   <- row_b_chr[, c("scenario", "method", q_cols_this)]
  
  ## ---- unbounded ----
  row_u_mean <- as.data.frame(t(df_tt$mean_missed_unbounded / n_big))
  row_u_sd   <- as.data.frame(t(df_tt$sd_missed_unbounded   / n_big))
  row_u_chr <- as.data.frame(
    setNames(
      lapply(seq_along(q_levels), function(j) {
        mm <- row_u_mean[1, j]
        ss <- row_u_sd[1, j]
        if (is.na(mm)) {
          NA_character_
        } else {
          sprintf("%.3f (%.3f)", mm, ss)
        }
      }),
      q_names
    )
  )
  row_u_chr$scenario <- tt
  row_u_chr$method   <- "unbounded"
  q_cols_this <- setdiff(colnames(row_u_chr), c("scenario", "method"))
  row_u_chr   <- row_u_chr[, c("scenario", "method", q_cols_this)]
  
  ## ---- coverage ----
  row_c_mean <- as.data.frame(t(df_tt$mean_missed_coverage / n_big))
  row_c_sd   <- as.data.frame(t(df_tt$sd_missed_coverage   / n_big))
  row_c_chr <- as.data.frame(
    setNames(
      lapply(seq_along(q_levels), function(j) {
        mm <- row_c_mean[1, j]
        ss <- row_c_sd[1, j]
        if (is.na(mm)) {
          NA_character_
        } else {
          sprintf("%.3f (%.3f)", mm, ss)
        }
      }),
      q_names
    )
  )
  row_c_chr$scenario <- tt
  row_c_chr$method   <- "coverage"
  q_cols_this <- setdiff(colnames(row_c_chr), c("scenario", "method"))
  row_c_chr   <- row_c_chr[, c("scenario", "method", q_cols_this)]
  
  ## ---- Chao2009 ----
  row_d_mean <- as.data.frame(t(df_tt$mean_missed_Chao2009 / n_big))
  row_d_sd   <- as.data.frame(t(df_tt$sd_missed_Chao2009   / n_big))
  row_d_chr <- as.data.frame(
    setNames(
      lapply(seq_along(q_levels), function(j) {
        mm <- row_d_mean[1, j]
        ss <- row_d_sd[1, j]
        if (is.na(mm)) {
          NA_character_
        } else {
          sprintf("%.3f (%.3f)", mm, ss)
        }
      }),
      q_names
    )
  )
  row_d_chr$scenario <- tt
  row_d_chr$method   <- "Chao2009"
  q_cols_this <- setdiff(colnames(row_d_chr), c("scenario", "method"))
  row_d_chr   <- row_d_chr[, c("scenario", "method", q_cols_this)]

  ## ---- Bounded new / Unbounded new ----
  make_new_missed_row <- function(mean_values, sd_values, method) {
    out <- as.data.frame(setNames(
      lapply(seq_along(q_levels), function(j) {
        mm <- mean_values[j] / n_big
        ss <- sd_values[j] / n_big
        if (is.na(mm)) NA_character_ else sprintf("%.3f (%.3f)", mm, ss)
      }),
      q_names
    ))
    out$scenario <- tt
    out$method <- method
    out[, c("scenario", "method", q_names)]
  }

  row_bn_chr <- make_new_missed_row(
    df_tt$mean_missed_bounded_new,
    df_tt$sd_missed_bounded_new,
    "bounded new"
  )
  row_un_chr <- make_new_missed_row(
    df_tt$mean_missed_unbounded_new,
    df_tt$sd_missed_unbounded_new,
    "unbounded new"
  )
  
  rows_missed[[length(rows_missed) + 1L]] <- row_b_chr
  rows_missed[[length(rows_missed) + 1L]] <- row_u_chr
  rows_missed[[length(rows_missed) + 1L]] <- row_c_chr
  rows_missed[[length(rows_missed) + 1L]] <- row_d_chr
  rows_missed[[length(rows_missed) + 1L]] <- row_bn_chr
  rows_missed[[length(rows_missed) + 1L]] <- row_un_chr
}

missed_table <- do.call(rbind, rows_missed)

missed_latex <- kable(
  missed_table,
  format  = "latex",
  booktabs = TRUE,
  caption = "Mean proportion of missed species with $p_j \\ge \\varepsilon$ (standard deviation in parentheses) by scenario, method, and contamination level $q$.",
  label   = "tab:missed_prop"
) %>%
  kable_styling(full_width = FALSE, position = "center", latex_options = c("hold_position"))

cat(missed_latex)


## ------------------------------------------------------------
## 4) Extra species table (mean (sd))
## ------------------------------------------------------------

rows_extra <- list()

for (tt in scenarios) {
  df_tt <- subset(summary_df, distro == tt)
  df_tt <- df_tt[match(q_levels, df_tt$q_error), ]
  q_names <- paste0("q_", gsub("\\.", "p", as.character(df_tt$q_error)))
  
  ## ---- bounded ----
  row_b_mean <- as.data.frame(t(df_tt$mean_extra_bounded))
  row_b_sd   <- as.data.frame(t(df_tt$sd_extra_bounded))
  row_b_chr <- as.data.frame(
    setNames(
      lapply(seq_along(q_levels), function(j) {
        mm <- row_b_mean[1, j]
        ss <- row_b_sd[1, j]
        if (is.na(mm)) {
          NA_character_
        } else {
          sprintf("%.2f (%.2f)", mm, ss)
        }
      }),
      q_names
    )
  )
  row_b_chr$scenario <- tt
  row_b_chr$method   <- "bounded"
  q_cols_this <- setdiff(colnames(row_b_chr), c("scenario", "method"))
  row_b_chr   <- row_b_chr[, c("scenario", "method", q_cols_this)]
  
  ## ---- unbounded ----
  row_u_mean <- as.data.frame(t(df_tt$mean_extra_unbounded))
  row_u_sd   <- as.data.frame(t(df_tt$sd_extra_unbounded))
  row_u_chr <- as.data.frame(
    setNames(
      lapply(seq_along(q_levels), function(j) {
        mm <- row_u_mean[1, j]
        ss <- row_u_sd[1, j]
        if (is.na(mm)) {
          NA_character_
        } else {
          sprintf("%.2f (%.2f)", mm, ss)
        }
      }),
      q_names
    )
  )
  row_u_chr$scenario <- tt
  row_u_chr$method   <- "unbounded"
  q_cols_this <- setdiff(colnames(row_u_chr), c("scenario", "method"))
  row_u_chr   <- row_u_chr[, c("scenario", "method", q_cols_this)]
  
  ## ---- coverage ----
  row_c_mean <- as.data.frame(t(df_tt$mean_extra_coverage))
  row_c_sd   <- as.data.frame(t(df_tt$sd_extra_coverage))
  row_c_chr <- as.data.frame(
    setNames(
      lapply(seq_along(q_levels), function(j) {
        mm <- row_c_mean[1, j]
        ss <- row_c_sd[1, j]
        if (is.na(mm)) {
          NA_character_
        } else {
          sprintf("%.2f (%.2f)", mm, ss)
        }
      }),
      q_names
    )
  )
  row_c_chr$scenario <- tt
  row_c_chr$method   <- "coverage"
  q_cols_this <- setdiff(colnames(row_c_chr), c("scenario", "method"))
  row_c_chr   <- row_c_chr[, c("scenario", "method", q_cols_this)]
  
  ## ---- Chao2009 ----
  row_d_mean <- as.data.frame(t(df_tt$mean_extra_Chao2009))
  row_d_sd   <- as.data.frame(t(df_tt$sd_extra_Chao2009))
  row_d_chr <- as.data.frame(
    setNames(
      lapply(seq_along(q_levels), function(j) {
        mm <- row_d_mean[1, j]
        ss <- row_d_sd[1, j]
        if (is.na(mm)) {
          NA_character_
        } else {
          sprintf("%.2f (%.2f)", mm, ss)
        }
      }),
      q_names
    )
  )
  row_d_chr$scenario <- tt
  row_d_chr$method   <- "Chao2009"
  q_cols_this <- setdiff(colnames(row_d_chr), c("scenario", "method"))
  row_d_chr   <- row_d_chr[, c("scenario", "method", q_cols_this)]

  ## ---- Bounded new / Unbounded new ----
  make_new_extra_row <- function(mean_values, sd_values, method) {
    out <- as.data.frame(setNames(
      lapply(seq_along(q_levels), function(j) {
        mm <- mean_values[j]
        ss <- sd_values[j]
        if (is.na(mm)) NA_character_ else sprintf("%.2f (%.2f)", mm, ss)
      }),
      q_names
    ))
    out$scenario <- tt
    out$method <- method
    out[, c("scenario", "method", q_names)]
  }

  row_bn_chr <- make_new_extra_row(
    df_tt$mean_extra_bounded_new,
    df_tt$sd_extra_bounded_new,
    "bounded new"
  )
  row_un_chr <- make_new_extra_row(
    df_tt$mean_extra_unbounded_new,
    df_tt$sd_extra_unbounded_new,
    "unbounded new"
  )
  
  rows_extra[[length(rows_extra) + 1L]] <- row_b_chr
  rows_extra[[length(rows_extra) + 1L]] <- row_u_chr
  rows_extra[[length(rows_extra) + 1L]] <- row_c_chr
  rows_extra[[length(rows_extra) + 1L]] <- row_d_chr
  rows_extra[[length(rows_extra) + 1L]] <- row_bn_chr
  rows_extra[[length(rows_extra) + 1L]] <- row_un_chr
}

extra_table <- do.call(rbind, rows_extra)

extra_latex <- kable(
  extra_table,
  format  = "latex",
  booktabs = TRUE,
  caption = "Mean number of extra species discovered at stopping (standard deviation in parentheses) by scenario, method, and contamination level $q$.",
  label   = "tab:extra"
) %>%
  kable_styling(full_width = FALSE, position = "center", latex_options = c("hold_position"))

cat(extra_latex)







library(ggplot2)

## ------------------------------------------------------------
## 1) Number of "big" species per scenario (p_j >= eps)
## ------------------------------------------------------------


## ------------------------------------------------------------
## 2) Build long data frame for plotting
##    metrics: N_stop, Missed proportion, Extra species
## ------------------------------------------------------------

q_levels <- sort(unique(summary_df$q_error))

plot_rows <- list()

add_method_rows <- function(df_tt, distro_name, n_big, method_name,
                            Nstop_col, missed_col, extra_col) {
  # df_tt: subset of summary_df for a given distro, ordered by q
  stopifnot(length(Nstop_col)  == nrow(df_tt))
  stopifnot(length(missed_col) == nrow(df_tt))
  stopifnot(length(extra_col)  == nrow(df_tt))
  
  data.frame(
    distro = distro_name,
    q      = df_tt$q_error,
    method = method_name,
    metric = rep(c("N_stop", "Missed proportion", "Extra species"),
                 each = nrow(df_tt)),
    value  = c(
      Nstop_col,
      missed_col / n_big,   # proportion of big species missed
      extra_col
    ),
    stringsAsFactors = FALSE
  )
}

for (d in distro_levels) {
  df_d <- subset(summary_df, distro == d)
  # ensure q is ordered
  df_d <- df_d[order(df_d$q_error), ]
  n_big <- n_big_vec[[d]]
  
  # bounded
  plot_rows[[length(plot_rows) + 1L]] <-
    add_method_rows(df_d, d, n_big,
                    "bounded",
                    df_d$mean_Nstop_bounded,
                    df_d$mean_missed_bounded,
                    df_d$mean_extra_bounded)
  
  # unbounded
  plot_rows[[length(plot_rows) + 1L]] <-
    add_method_rows(df_d, d, n_big,
                    "unbounded",
                    df_d$mean_Nstop_unbounded,
                    df_d$mean_missed_unbounded,
                    df_d$mean_extra_unbounded)
  
  # coverage
  plot_rows[[length(plot_rows) + 1L]] <-
    add_method_rows(df_d, d, n_big,
                    "coverage",
                    df_d$mean_Nstop_coverage,
                    df_d$mean_missed_coverage,
                    df_d$mean_extra_coverage)
  
  # Chao2009
  plot_rows[[length(plot_rows) + 1L]] <-
    add_method_rows(df_d, d, n_big,
                    "Chao2009",
                    df_d$mean_Nstop_Chao2009,
                    df_d$mean_missed_Chao2009,
                    df_d$mean_extra_Chao2009)

  # bounded new
  plot_rows[[length(plot_rows) + 1L]] <-
    add_method_rows(df_d, d, n_big,
                    "bounded new",
                    df_d$mean_Nstop_bounded_new,
                    df_d$mean_missed_bounded_new,
                    df_d$mean_extra_bounded_new)

  # unbounded new
  plot_rows[[length(plot_rows) + 1L]] <-
    add_method_rows(df_d, d, n_big,
                    "unbounded new",
                    df_d$mean_Nstop_unbounded_new,
                    df_d$mean_missed_unbounded_new,
                    df_d$mean_extra_unbounded_new)
}

plot_df <- do.call(rbind, plot_rows)

plot_df$metric <- factor(
  plot_df$metric,
  levels = c("N_stop", "Missed proportion", "Extra species"),
  labels = c(
    expression(N[stop]),
    "Missed big species (proportion)",
    "Extra species"
  )
)

plot_df$method <- factor(
  plot_df$method,
  levels = c("bounded", "unbounded", "coverage", "Chao2009",
             "bounded new", "unbounded new"),
  labels = c("Bounded CI", "Unbounded CI", "Coverage", "Chao2009",
             "Bounded new", "Unbounded new")
)

## Nice facet labels
plot_df$metric <- factor(
  plot_df$metric,
  levels = c("N[stop]", "Missed big species (proportion)", "Extra species"),
  labels = c("N stop", "% missed", "num. extra")
)

plot_df$distro <- factor(
  plot_df$distro,
  levels = c("zipf_heavy", "geom_0.05", "unif_small", "unif_big"),
  labels = c("Zipf, α = 1.05", "Geometric(0.05)", "p = 0.006", "p = 0.05")
)


## ------------------------------------------------------------
## 3) Plot: rows = scenarios, columns = metrics
## ------------------------------------------------------------

sqrt_transform <- function(y) sign(y) * sqrt(abs(y))

inv_sqrt_label <- function(y) {
  sign(y) * (abs(y)^2)
}

df_main <- plot_df |> dplyr::filter(metric != "% missed")
df_missed <- plot_df |> dplyr::filter(metric == "% missed") |>
  dplyr::mutate(value_trans = sqrt_transform(value))

p1 <- ggplot(df_main, aes(x = q, y = value, colour = method)) +
  geom_line() +
  geom_point(size = 1.5) +
  facet_grid(metric ~ distro, scales = "free_y") +
  labs(
    x = expression(q),
    y = NULL,
    colour = "Stopping rule"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "none",      # legend only in p2
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p2 <- ggplot(df_missed, aes(x = q, y = value_trans, colour = method)) +
  geom_line() +
  geom_point(size = 1.5) +
  facet_grid(metric ~ distro, scales = "free_y") +
  labs(
    x = expression(q),
    y = NULL,
    colour = "Stopping rule"
  ) +
  scale_y_continuous(labels = inv_sqrt_label) +
  theme_bw() +
  theme(
    # Remove only column strips (top)
    strip.background.x = element_blank(),
    strip.text.x       = element_blank(),
    strip.background.y = element_rect(fill = "grey90"),
    strip.text.y       = element_text(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

library(patchwork)

p2_adj <- p2 + theme(plot.margin = margin(t = -30, r = 5, b = 0, l = 5))

p_final <- p1 / p2_adj + 
  plot_layout(heights = c(2, 1))

print(p_final)

ggsave("stopping_rules_figure.pdf", p_final, width = 9, height = 7)

experiment_end_time <- Sys.time()
cat(sprintf(
  "[%s] EXPERIMENT DONE | total elapsed=%s | output=stopping_rules_figure.pdf\n",
  format(experiment_end_time, "%Y-%m-%d %H:%M:%S"),
  format_elapsed_time(as.numeric(difftime(
    experiment_end_time, experiment_start_time, units = "secs"
  )))
))
flush.console()

