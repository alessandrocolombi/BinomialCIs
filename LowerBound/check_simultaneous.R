# Deterministic checks only: no Monte Carlo sampling.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
source(file.path(script_dir, "lower_bound.R"))
source(file.path(script_dir, "results_report.R"))

# Empty spectra, fully observed features, and a nonzero-endpoint fixture.
stopifnot(lower_bound_mmax_simultaneous(12, integer()) == 0,
          lower_bound_mmax_simultaneous(12, rep(12, 10)) == 0)
counts <- rep(1:12, each = 10000)
fit <- lower_bound_mmax_simultaneous(2000, counts, details = TRUE)
stopifnot(fit$lower > 0, fit$lower <= 1,
          fit$lower == max(fit$by_r), fit$selected_r %in% 2:6,
          identical(fit$lower, lower_bound_mmax_simultaneous(2000, c(counts, rep(0, 10)))))

# The rationalized root agrees with numerical inversion at diverse scales.
for (scale in c(1, 1e-15, 1e-100)) {
  T <- c(2, 1) * scale
  B <- c(0.1, 0.1) * scale
  V <- c(0.005, 0.003) * scale^2
  cvalue <- 3
  root <- unname(.lb_joint_roots(T, B, V, cvalue))
  score <- function(t) (T[2] - B[2] - t * T[1] -
                         cvalue * sqrt(V[2] + t^2 * V[1])) / scale
  ref <- uniroot(score, c(0, 1), tol = 1e-14)$root
  stopifnot(abs(root - ref) < 1e-12, score(root) >= -1e-14)
}

# Exact enumeration validates the common-event argument without drawing data.
n <- 12L; R <- 6L; alpha <- 0.05; p <- c(0.02, 0.12, 0.4)
states <- as.matrix(expand.grid(rep(list(0:n), length(p))))
prob <- rep(1, nrow(states))
for (j in seq_along(p)) prob <- prob * dbinom(states[, j], n, p[j])
eta <- vapply(1:12, function(j) sum(dbinom(j, n, p)), numeric(1))
a <- exp(-lchoose(n, 1:12))
D <- M <- matrix(0, nrow(states), R)
for (j in 1:R) {
  M[, j] <- as.vector((states == 0) %*% p^j)
  D[, j] <- a[j] * rowSums(states == j) - M[, j]
}
b <- colSums(D * prob)
v <- a[1:R]^2 * eta[1:R] + a[2*(1:R)] * eta[2*(1:R)]
Q <- rowSums(sweep(sweep(D, 2, b), 2, sqrt(v), "/")^2)
stopifnot(sum(prob[Q > 2*R/alpha]) <= alpha/2 + 1e-12)
Mmax <- apply(sweep(states == 0, 2, p, "*"), 1, max)
lower <- apply(states, 1, function(x) lower_bound_mmax_simultaneous(n, x, alpha, R))
stopifnot(sum(prob[lower > Mmax]) <= alpha + 1e-12)

# Saved spectra can enrich a cell without regenerating observations.
cell <- list(task = data.frame(n=2000), summary = data.frame(alpha=alpha),
             runs = matrix(0, 1, 2, dimnames=list(NULL,c("lower","Mmax"))),
             spectrum = matrix(rep(10000,12), nrow=1))
enriched <- lb_add_joint_results(cell, alpha)
stopifnot("lower_joint" %in% colnames(enriched$runs),
          enriched$summary$lb_joint_median == fit$lower)
cat("Passed: root inversion, edge cases, spectrum replay and exact enumeration.\nFINITO\n")
