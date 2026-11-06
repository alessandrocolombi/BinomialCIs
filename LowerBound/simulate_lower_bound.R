# Run from any working directory:
# Rscript LowerBound/simulate_lower_bound.R [Brep=1000] [output_dir] [seed=42]
# In RStudio: source("LowerBound/simulate_lower_bound.R").
# Only base R is required; no Rcpp compilation or parallel workers.

.simulation_script_dir <- function() {
  # Prefer source()'s filename when this script is sourced from another script.
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    if (!is.null(frames[[i]]$ofile)) {
      return(dirname(normalizePath(frames[[i]]$ofile, mustWork = TRUE)))
    }
  }
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]),
                                 mustWork = TRUE)))
  }
  stop("Run with Rscript or source() so the script location can be resolved.")
}

script_dir <- .simulation_script_dir()
repo_dir <- dirname(script_dir)
source(file.path(script_dir, "lower_bound.R"))
source(file.path(repo_dir, "R", "Rfunctions.R"))

args <- if (sys.nframe() == 0L) commandArgs(trailingOnly = TRUE) else character()
if (length(args) > 3L) stop("Usage: simulate_lower_bound.R [Brep] [output_dir] [seed]")
Brep <- if (length(args) >= 1L) as.numeric(args[1]) else 1000
output_dir <- if (length(args) >= 2L) args[2] else file.path(script_dir, "results")
seed <- if (length(args) >= 3L) as.numeric(args[3]) else 42
if (length(Brep) != 1L || !is.finite(Brep) || Brep < 1 || Brep != floor(Brep)) {
  stop("Brep must be a positive integer.")
}
if (length(seed) != 1L || !is.finite(seed) || seed < 0 ||
    seed > .Machine$integer.max || seed != floor(seed)) stop("Invalid seed.")
alpha <- 0.05
set.seed(seed)

# Each pair belongs to an actual arm of Exp1_varying_M_n.R, not a cross-grid.
settings <- data.frame(M = c(2000L, 10000L, 5000L, 5000L),
                       n = c(2000L, 2000L, 1000L, 5000L))
# Parameters match Exp1_varying_M_n_parallel_main.R.
benchmarks <- data.frame(
  distribution = rep(c("Zipfs", "Geom", "Constant"), each = 3),
  param = c(0.25, 0.5, 1.02, 0.005, 0.1, 0.25, 0.001, 0.05, 0.5),
  label = c("Z 0.25", "Z 0.5", "Z 1.02", "G 0.005", "G 0.1", "G 0.25",
            "C 0.001", "C 0.05", "C 0.5"), stringsAsFactors = FALSE)

mean_or_na <- function(x) if (length(x)) mean(x) else NA_real_
summary_rows <- list()
replicate_rows <- list()
cell <- 0L
started <- proc.time()["elapsed"]

for (s in seq_len(nrow(settings))) {
  M <- settings$M[s]
  n <- settings$n[s]
  for (d in seq_len(nrow(benchmarks))) {
    cell <- cell + 1L
    name <- benchmarks$distribution[d]
    param <- benchmarks$param[d]
    # Fixed prevalences across repetitions, using the existing generators.
    # Finite support is an exact special case of the infinite-alphabet model:
    # all probabilities beyond M are zero. No truncation approximation is used.
    p <- sim_generic(name, M, param)
    cell_rows <- vector("list", Brep)
    for (b in seq_len(Brep)) {
      counts <- rbinom(M, size = n, prob = p)
      unseen <- counts == 0
      Mmax <- if (any(unseen)) max(p[unseen]) else 0
      M1 <- sum(p[unseen])
      M2 <- sum(p[unseen]^2)
      # Deliberately give the estimator only observed counts, not M or p.
      fit <- lower_bound_mmax(n, counts[counts > 0], alpha, details = TRUE)
      oracle_ratio <- if (M1 > 0) M2 / M1 else 0
      cell_rows[[b]] <- data.frame(
        setting = s, benchmark = benchmarks$label[d], distribution = name,
        param = param, M = M, n = n, repetition = b, S = sum(p),
        Mmax = Mmax, M1 = M1, M2 = M2, oracle_ratio = oracle_ratio,
        lower = fit$lower, M2_lower = fit$M2_lower, M1_upper = fit$M1_upper,
        M2_lower_raw = fit$M2_lower_raw, M2_hat = fit$M2_hat,
        bias_upper = fit$bias_upper, variance_upper = fit$variance_upper,
        tail_radius = fit$tail_radius,
        Kn = unname(fit$counts["Kn"]), K1 = unname(fit$counts["K1"]),
        K2 = unname(fit$counts["K2"]), K3 = unname(fit$counts["K3"]),
        K4 = unname(fit$counts["K4"]),
        covered = fit$lower <= Mmax,
        M2_covered = fit$M2_lower <= M2,
        M1_covered = M1 <= fit$M1_upper)
    }
    runs <- do.call(rbind, cell_rows)
    has_unseen <- runs$Mmax > 0
    positive <- runs$lower > 0
    coverage <- mean(runs$covered)
    ci <- binom.test(sum(runs$covered), Brep)$conf.int
    quantiles <- quantile(runs$lower, c(0.025, 0.5, 0.975), names = FALSE)
    summary_rows[[cell]] <- data.frame(
      setting = s, benchmark = benchmarks$label[d], distribution = name,
      param = param, M = M, n = n, Brep = Brep, alpha = alpha, S = sum(p),
      coverage = coverage, coverage_mcse = sqrt(coverage * (1 - coverage) / Brep),
      coverage_ci_low = ci[1], coverage_ci_high = ci[2],
      failure_count = sum(!runs$covered), unseen_rate = mean(has_unseen),
      coverage_given_unseen = mean_or_na(runs$covered[has_unseen]),
      positive_rate = mean(positive),
      positive_rate_given_unseen = mean_or_na(positive[has_unseen]),
      mean_lower = mean(runs$lower), q025_lower = quantiles[1],
      median_lower = quantiles[2], q975_lower = quantiles[3],
      mean_Mmax = mean(runs$Mmax), mean_oracle_ratio = mean(runs$oracle_ratio),
      mean_relative_lower = mean_or_na(runs$lower[has_unseen] / runs$Mmax[has_unseen]),
      mean_relative_oracle = mean_or_na(runs$oracle_ratio[has_unseen] / runs$Mmax[has_unseen]),
      M2_coverage = mean(runs$M2_covered), M1_coverage = mean(runs$M1_covered))
    replicate_rows[[cell]] <- runs
    cat(sprintf("[%2d/%d] M=%d n=%d %-9s p=%g: coverage=%.3f, positive=%.3f\n",
                cell, nrow(settings) * nrow(benchmarks), M, n, name, param,
                coverage, mean(positive)))
  }
}

summary_results <- do.call(rbind, summary_rows)
replicate_results <- do.call(rbind, replicate_rows)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(summary_results, file.path(output_dir, "summary.csv"), row.names = FALSE)
saveRDS(list(config = list(Brep = Brep, alpha = alpha, seed = seed,
                          settings = settings, benchmarks = benchmarks,
                          RNGkind = RNGkind(), session = sessionInfo()),
             summary = summary_results, replicates = replicate_results),
        file.path(output_dir, "simulation.rds"))

# Coverage, non-vacuity, and tightness are separate diagnostics. In particular,
# coverage=1 with Mmax=0 in every repetition is not evidence of a useful bound.
png(file.path(output_dir, "diagnostics.png"), width = 2400, height = 2200, res = 180)
par(mfrow = c(nrow(settings), 3), mar = c(5, 4, 3, 1), mgp = c(2.4, 0.6, 0))
for (s in seq_len(nrow(settings))) {
  z <- summary_results[summary_results$setting == s, ]
  xx <- seq_len(nrow(z))
  title <- sprintf("M=%d, n=%d", settings$M[s], settings$n[s])
  plot_panel <- function(y, ylab, suffix) {
    plot(xx, y, xaxt = "n", xlim = c(0.5, length(xx) + 0.5), ylim = c(0, 1.03),
         xlab = "", ylab = ylab, main = paste(title, suffix), pch = 16)
    axis(1, at = xx, labels = z$benchmark, las = 2, cex.axis = 0.75)
  }
  plot_panel(z$coverage, "Coverage (95% MC interval)", "- coverage")
  arrows(xx, z$coverage_ci_low, xx, z$coverage_ci_high,
         angle = 90, code = 3, length = 0.025, col = "grey40")
  abline(h = 1 - alpha, col = "firebrick", lty = 2)
  plot_panel(z$positive_rate, "Proportion of repetitions", "- nonzero bound")
  points(xx, z$unseen_rate, pch = 1, col = "steelblue")
  legend("bottomleft", c("Lower > 0", "Mmax > 0"), pch = c(16, 1),
         col = c("black", "steelblue"), bty = "n", cex = 0.75)
  # Initialize with a finite vector even if every conditional ratio is NA.
  plot_panel(rep(NA_real_, length(xx)), "Mean ratio to Mmax, given Mmax > 0",
             "- tightness")
  points(xx, z$mean_relative_lower, pch = 16)
  points(xx, z$mean_relative_oracle, pch = 1, col = "darkgreen")
  legend("bottomleft", c("Confidence lower", "True M2 / M1"), pch = c(16, 1),
         col = c("black", "darkgreen"), bty = "n", cex = 0.75)
}
invisible(dev.off())
cat(sprintf("\nSaved %d cells and %d repetitions to %s (%.1f seconds).\n",
            nrow(summary_results), nrow(replicate_results),
            normalizePath(output_dir), proc.time()["elapsed"] - started))
