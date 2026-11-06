# Mathematical and edge-case checks; base R only.
# Run: Rscript LowerBound/check_lower_bound.R
args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(args)) stop("Run this check script with Rscript.")
script_dir <- dirname(normalizePath(sub("^--file=", "", args[1])))
source(file.path(script_dir, "lower_bound.R"))

# A non-vacuous fixture; independent calculations of the proposed formula.
n <- 2000
counts <- rep(1:4, c(1300, 1400, 800, 500))
fit <- lower_bound_mmax(n, counts, details = TRUE)
eps <- 0.05 / 5
x <- -log(eps)
C <- (sqrt(x / 2) + sqrt(c(1400, 800, 500) + 7 * x / 6))^2
bb <- 2 * C[2] / choose(n, 3)
vv <- C[1] / choose(n, 2)^2 + C[3] / choose(n, 4)
L2 <- max(0, 1400 / choose(n, 2) - bb - sqrt(2 * vv * x) -
            2 * x / (3 * choose(n, 2)))
y <- log(2 / eps)
CK <- (sqrt(y / 2) + sqrt(length(counts) + 7 * y / 6))^2
U1 <- 1300 / n + y / (n - 1) +
  sqrt(8 * y * CK / ((n - 1)^2 * (1 - 2 / n)))
stopifnot(fit$lower > 0,
          isTRUE(all.equal(fit$M2_lower, L2)),
          isTRUE(all.equal(fit$M1_upper, U1)),
          isTRUE(all.equal(fit$lower, L2 / U1)),
          isTRUE(all.equal(sum(fit$error_budget), 0.05)),
          identical(fit$lower, lower_bound_mmax(n, c(counts, rep(0, 100)))),
          identical(fit$lower, lower_bound_mmax(n, rev(counts))))

# Empty observations, fully observed populations, and the smallest valid n.
stopifnot(lower_bound_mmax(4, integer()) == 0,
          lower_bound_mmax(4, rep(0, 10)) == 0,
          lower_bound_mmax(4, rep(4, 10)) == 0)
fails <- function(expr) inherits(tryCatch(expr, error = identity), "error")
stopifnot(fails(lower_bound_mmax(3, 1)),
          fails(lower_bound_mmax(4, c(0, 5))),
          fails(lower_bound_mmax(4, 1.5)),
          fails(lower_bound_mmax(4, NA_real_)),
          fails(lower_bound_mmax(4, 1, alpha = 0)),
          fails(lower_bound_mmax(4, 1, alpha = 1)))

# Exact enumeration of a heterogeneous three-feature model. Verify the bias,
# variance, deterministic moment inequality, and Bernstein tail with no MC.
n <- 20
p <- c(0.01, 0.08, 0.25)
states <- as.matrix(expand.grid(rep(list(0:n), length(p))))
prob <- rep(1, nrow(states))
for (j in seq_along(p)) prob <- prob * dbinom(states[, j], n, p[j])
unseen <- states == 0
M1 <- as.vector(unseen %*% p)
M2 <- as.vector(unseen %*% p^2)
Mmax <- apply(sweep(unseen, 2, p, "*"), 1, max)
D2 <- rowSums(states == 2) / choose(n, 2) - M2
eta <- vapply(2:4, function(k) sum(dbinom(k, n, p)), numeric(1))
bias <- sum(prob * D2)
variance <- sum(prob * (D2 - bias)^2)
b <- 2 * eta[2] / choose(n, 3)
v <- eta[1] / choose(n, 2)^2 + eta[3] / choose(n, 4)
stopifnot(abs(sum(prob) - 1) < 1e-12, bias >= 0, bias <= b,
          variance <= v, all(M2 <= Mmax * M1 + 1e-14))
for (epsilon in c(0.01, 0.05, 0.2, 0.5)) {
  x <- -log(epsilon)
  threshold <- b + sqrt(2 * v * x) + 2 * x / (3 * choose(n, 2))
  stopifnot(sum(prob[D2 > threshold]) <= epsilon + 1e-12)
}
cat("Passed: formula, nonzero fixture, zero-count invariance, input checks,\n",
    "and exact bias/variance/tail checks on 9261 possible samples.\n", sep = "")
