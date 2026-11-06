# Lower confidence bound for the largest unseen prevalence.
# Model: independent N_j ~ Bin(n, p_j), fixed p, sum(p_j) < Inf.
# counts may contain all counts or only the positive counts; M is not used.

.lb_check_level <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a number strictly between zero and one.")
  }
}

.lb_check_data <- function(n, counts) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 4 || n != floor(n)) {
    stop("n must be an integer >= 4.")
  }
  if (!is.numeric(counts) || !is.null(dim(counts)) ||
      any(!is.finite(counts)) || any(counts < 0 | counts > n) ||
      any(counts != floor(counts))) {
    stop("counts must be a numeric vector of integer counts between 0 and n.")
  }
}

# Ayed et al. (2019), Remark A.1: an upper confidence bound for E[K].
# Applies to K_n and K_{n,k}, each a sum of independent Bernoulli variables.
ayed_expected_count_upper <- function(K, alpha) {
  .lb_check_level(alpha)
  if (!is.numeric(K) || any(!is.finite(K)) || any(K < 0) ||
      any(K != floor(K))) {
    stop("K must contain nonnegative integer counts.")
  }
  x <- -log(alpha)
  (sqrt(x / 2) + sqrt(K + 7 * x / 6))^2
}

# One-sided version obtained explicitly from Proposition A.3 / (A.22)-(A.23).
# Split alpha between the count-concentration and missing-mass tail events.
# This makes the error allocation explicit rather than relying on the
# simultaneous two-sided interval in Theorem 3.1.
ayed_missing_mass_upper <- function(n, counts, alpha = 0.05) {
  .lb_check_data(n, counts)
  .lb_check_level(alpha)
  x <- log(2) - log(alpha)
  C <- ayed_expected_count_upper(sum(counts > 0), alpha / 2)
  sum(counts == 1) / n + x / (n - 1) +
    sqrt(8 * x * C / ((n - 1)^2 * (1 - 2 / n)))
}

# With probability >= 1-alpha, Mmax >= the returned endpoint.
# Five outer events, each with error alpha/5:
#   E[K_{n,2}], E[K_{n,3}], E[K_{n,4}] upper bounds;
#   the right tail of D_2 = K_{n,2}/choose(n,2) - M_2;
#   the upper bound for M_1 (internally split into two events).
# No independence between these events is assumed.
lower_bound_mmax <- function(n, counts, alpha = 0.05, details = FALSE) {
  .lb_check_data(n, counts)
  .lb_check_level(alpha)
  if (!is.logical(details) || length(details) != 1L || is.na(details)) {
    stop("details must be TRUE or FALSE.")
  }

  epsilon <- alpha / 5
  x <- -log(epsilon)
  K <- c(Kn = sum(counts > 0), K1 = sum(counts == 1),
         K2 = sum(counts == 2), K3 = sum(counts == 3),
         K4 = sum(counts == 4))
  C <- ayed_expected_count_upper(K[c("K2", "K3", "K4")], epsilon)
  # Inverse binomial coefficients evaluated on the log scale.
  a <- exp(-lchoose(n, 2))
  bias_upper <- unname(2 * C["K3"] * exp(-lchoose(n, 3)))
  variance_upper <- unname(a^2 * C["K2"] +
                            C["K4"] * exp(-lchoose(n, 4)))
  tail_radius <- sqrt(2 * variance_upper * x) + 2 * a * x / 3
  M2_hat <- unname(a * K["K2"])
  M2_lower_raw <- M2_hat - bias_upper - tail_radius
  M2_lower <- max(0, M2_lower_raw)
  M1_upper <- ayed_missing_mass_upper(n, counts, epsilon)
  lower <- if (M1_upper > 0) M2_lower / M1_upper else 0

  # Do not clip at one: the formula itself should be checked in simulations.
  if (!details) return(lower)
  list(lower = lower, M2_lower = M2_lower, M1_upper = M1_upper,
       M2_lower_raw = M2_lower_raw, M2_hat = M2_hat,
       bias_upper = bias_upper, variance_upper = variance_upper,
       tail_radius = tail_radius, expected_counts_upper = C,
       counts = K, n = n, alpha = alpha,
       error_budget = c(K2 = epsilon, K3 = epsilon, K4 = epsilon,
                        D2_tail = epsilon, M1_upper = epsilon))
}

# Simultaneous diagonal-ellipsoid theorem (thm:simultaneous) in Notes/theory.tex.
# Each threshold is evaluated by a rationalized root, rounded downward.
.lb_joint_roots <- function(T, B, V, c) {
  r <- 2:length(T)
  A <- T[r] - B[r]
  D <- T[r - 1L]
  z <- c * sqrt(V[r])
  w <- c * sqrt(V[r - 1L])
  result <- numeric(length(r))
  ok <- which(A > z)
  if (length(ok)) {
    scale <- pmax(A[ok], D[ok], z[ok], w[ok])
    aa <- A[ok] / scale; dd <- D[ok] / scale
    vv <- (z[ok] / scale)^2; ww <- (w[ok] / scale)^2
    N <- (aa - sqrt(vv)) * (aa + sqrt(vv))
    denominator <- aa * dd + sqrt(dd^2 * vv + ww * N)
    root <- N / denominator
    root[denominator == 0] <- 1
    result[ok] <- pmin(1, root) * (1 - 16 * .Machine$double.eps)
    # Return a feasible side of the root even in unusual floating-point cases.
    for (i in ok) {
      score <- function(t) A[i] - t * D[i] - sqrt(z[i]^2 + (t * w[i])^2)
      if (!is.finite(result[i]) || score(result[i]) < 0) {
        lo <- 0; hi <- 1
        for (step in 1:60) {
          mid <- (lo + hi) / 2
          if (score(mid) > 0) lo <- mid else hi <- mid
        }
        result[i] <- lo
      }
    }
  }
  names(result) <- paste0("r", r)
  result
}

lower_bound_mmax_simultaneous_spectrum <- function(n, K, alpha = 0.05, R = 6L,
                                                  details = FALSE) {
  .lb_check_data(n, integer())
  .lb_check_level(alpha)
  if (!is.numeric(R) || length(R) != 1L || !is.finite(R) || R < 2 ||
      R != floor(R) || 2 * R > n) stop("R must be an integer >= 2 with 2*R <= n.")
  if (!is.numeric(K) || length(K) != 2 * R || any(!is.finite(K)) ||
      any(K < 0 | K != floor(K))) stop("K must contain K_1,...,K_(2R).")
  if (!is.logical(details) || length(details) != 1L || is.na(details)) stop("Invalid details.")
  j <- seq_len(R)
  a <- exp(-lchoose(n, seq_len(R + 1L)))
  C <- ayed_expected_count_upper(K, alpha / (4 * R))
  T <- a[j] * K[j]
  B <- j * C[j + 1L] * a[j + 1L]
  V <- a[j]^2 * C[j] + C[2L * j] * exp(-lchoose(n, 2L * j))
  by_r <- .lb_joint_roots(T, B, V, sqrt(2 * R / alpha))
  selected <- which.max(by_r)
  lower <- unname(by_r[selected])
  if (!details) return(lower)
  list(lower = lower, by_r = by_r, selected_r = selected + 1L,
       T = T, bias_upper = B, variance_upper = V, expected_counts_upper = C,
       K = K, R = R, alpha = alpha, version = "simultaneous-ellipse-v1")
}

lower_bound_mmax_simultaneous <- function(n, counts, alpha = 0.05, R = 6L,
                                         details = FALSE) {
  .lb_check_data(n, counts)
  # Validate R before using it as a vector length.
  if (!is.numeric(R) || length(R) != 1L || !is.finite(R) || R < 2 ||
      R != floor(R) || 2 * R > n) stop("R must be an integer >= 2 with 2*R <= n.")
  K <- tabulate(counts[counts > 0 & counts <= 2 * R], nbins = 2 * R)
  lower_bound_mmax_simultaneous_spectrum(n, K, alpha, R, details)
}
