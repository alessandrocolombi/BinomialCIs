# VM launch, from the repository root:
# ./run_job.sh LowerBound/Exp_LB_M_n.R choose_wd=2
# Optional: wd=/absolute/path/to/BinomialCIs workers=15 Brep=1000 seed=42
# Rebuild figures/tables without simulating: choose_wd=2 report_only=TRUE
# This script starts the experiment when launched; prepared for the VM.
# Only base R and parallel (bundled with R) are required.

main <- function() {
  # Working-directory profiles. The VM profile uses this checkout's location,
  # so it also works when the VM user's home is not /home/colombi.
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list()
  allowed <- c("choose_wd", "wd", "workers", "Brep", "seed", "alpha", "report_only")
  for (arg in args) {
    eq <- regexpr("=", arg, fixed = TRUE)[1]
    if (eq <= 1) stop("Arguments must have the form key=value: ", arg)
    key <- substr(arg, 1, eq - 1)
    value <- substring(arg, eq + 1)
    if (!key %in% allowed || !nzchar(value)) stop("Invalid argument: ", arg)
    if (!is.null(opts[[key]])) stop("Repeated argument: ", key)
    opts[[key]] <- value
  }
  integer_option <- function(key, default, minimum = 1) {
    value <- if (is.null(opts[[key]])) default else suppressWarnings(as.numeric(opts[[key]]))
    if (length(value) != 1L || !is.finite(value) || value != floor(value) ||
        value < minimum || value > .Machine$integer.max) {
      stop(key, " must be an integer >= ", minimum)
    }
    as.integer(value)
  }
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) != 1L) stop("Launch this script with Rscript or run_job.sh.")
  script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
  checkout <- dirname(dirname(script_path))
  wd_vec <- c(pc = "C:/Users/colom/BinomialCIs", vm = checkout)
  choose_wd <- integer_option("choose_wd", 2L)
  if (choose_wd > length(wd_vec)) stop("choose_wd must be 1 (PC) or 2 (VM).")
  wd <- if (is.null(opts$wd)) wd_vec[choose_wd] else opts$wd
  wd <- normalizePath(wd, mustWork = TRUE)
  if (!file.exists(file.path(wd, "LowerBound", "lower_bound.R")) ||
      !file.exists(file.path(wd, "R", "Rfunctions.R"))) {
    stop("wd must be the BinomialCIs repository root: ", wd)
  }
  setwd(wd)
  report_only <- if (is.null(opts$report_only)) FALSE else {
    value <- tolower(opts$report_only)
    if (!value %in% c("true", "false")) stop("report_only must be TRUE or FALSE.")
    value == "true"
  }
  if (report_only && any(c("Brep", "alpha", "seed", "workers") %in% names(opts))) {
    stop("report_only uses the saved configuration; omit Brep, alpha, seed and workers.")
  }

  # The master also needs get_first3digits() for output names. Previously
  # Rfunctions.R was sourced only on workers, which do not share environments.
  # Keep a local copy of the user's naming helper to avoid that dependency.
  get_first3digits <- function(x, d = 3) {
    s <- format(x, scientific = FALSE, trim = TRUE)
    s <- gsub("\\.", "", s)
    substr(s, 1, d)
  }

  # Defaults retained from the preceding lower-bound simulation.
  Brep <- integer_option("Brep", 1000L)
  seed <- integer_option("seed", 42L, minimum = 0)
  alpha <- if (is.null(opts$alpha)) 0.05 else suppressWarnings(as.numeric(opts$alpha))
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be strictly between zero and one.")
  }
  cores <- parallel::detectCores(logical = FALSE)
  if (is.na(cores)) cores <- 2L
  # Leave one core free, with the same 15-worker ceiling as other VM jobs.
  workers <- integer_option("workers", min(15L, max(1L, cores - 1L)))

  # Original Exp1 grids: two one-dimensional arms, not a Cartesian product.
  n_fixed <- 2000L
  Mgrid <- seq(2000L, 10000L, by = 250L)
  M_fixed <- 5000L
  Ngrid <- seq(1000L, 5000L, by = 250L)
  benchmarks <- data.frame(
    distribution = rep(c("Zipfs", "Geom", "Constant"), each = 3),
    param = c(0.25, 0.5, 1.02, 0.005, 0.1, 0.25, 0.001, 0.05, 0.5),
    stringsAsFactors = FALSE)
  # Resolve names before starting the expensive parallel experiment.
  benchmarks$tag <- vapply(benchmarks$param, get_first3digits, character(1), d = 4)
  if (anyDuplicated(paste(benchmarks$distribution, benchmarks$tag))) {
    stop("Benchmark parameters produce duplicate output names.")
  }
  grid <- rbind(data.frame(arm = "nfix", M = Mgrid, n = n_fixed),
                data.frame(arm = "Mfix", M = M_fixed, n = Ngrid))
  tasks <- merge(benchmarks, grid, by = NULL, sort = FALSE)
  tasks$task_id <- seq_len(nrow(tasks))
  workers <- min(workers, nrow(tasks))

  lower_dir <- file.path(wd, "LowerBound")
  for (name in unique(benchmarks$distribution)) {
    dir.create(file.path(lower_dir, "img", name), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(lower_dir, "save", name), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(lower_dir, "save", "Coverage", name),
               recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(lower_dir, "res", name), recursive = TRUE, showWarnings = FALSE)
  }

  # Independent reproducible L'Ecuyer streams per task. Results do not depend
  # on worker count or dynamic scheduling. Avoid O(number of workers) copies
  # of a full repetition-by-alphabet matrix: each worker samples one vector.
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  streams <- vector("list", nrow(tasks))
  streams[[1]] <- .Random.seed
  if (length(streams) > 1L) {
    for (i in 2:length(streams)) streams[[i]] <- parallel::nextRNGStream(streams[[i - 1L]])
  }
  jobs <- lapply(seq_len(nrow(tasks)), function(i) {
    list(task = tasks[i, , drop = FALSE], rng = streams[[i]])
  })

  # Exact oracle: largest deterministic lower endpoint available when p is
  # known, at the inclusive lower-coverage level 1-alpha. For finite support,
  # Pr(Mmax <= t) = product_{j:p_j>t} [1 - (1-p_j)^n].
  # Use the alpha quantile on the discrete support, not the 1-alpha quantile
  # used for an UPPER bound, and not the moment ratio M2/M1.
  oracle_lower <- function(p, n, alpha) {
    positive <- sort(p[p > 0], decreasing = TRUE)
    if (!length(positive)) return(0)
    groups <- rle(positive)
    values <- groups$values
    log_seen <- log(-expm1(n * log1p(-values)))
    group_log <- groups$lengths * log_seen
    cumulative <- cumsum(group_log)
    # Strict lower tails at each descending support point.
    log_strict <- cumulative
    # The largest support point with Pr(Mmax < t) <= alpha is optimal for
    # the inclusive guarantee Pr(Mmax >= t) >= 1-alpha. This also handles
    # exact equality at an atom (where a standard alpha quantile can be smaller).
    eligible <- log_strict <= log(alpha)
    if (any(eligible)) return(values[which(eligible)[1L]])
    0
  }

  simulate_cell <- function(job, Brep, alpha, oracle_lower) {
    assign(".Random.seed", job$rng, envir = .GlobalEnv)
    task <- job$task
    M <- task$M
    n <- task$n
    p <- sim_generic(task$distribution, M, task$param)
    if (length(p) != M || any(!is.finite(p)) || any(p < 0 | p > 1)) {
      stop("Invalid benchmark prevalences in task ", task$task_id)
    }
    oracle <- oracle_lower(p, n, alpha)
    runs <- matrix(NA_real_, Brep, 12L, dimnames = list(NULL,
      c("lower", "Mmax", "M1", "M2", "moment_ratio", "M1_upper", "M2_lower",
        "M2_hat", "M2_lower_raw", "bias_upper", "variance_upper", "tail_radius")))
    for (b in seq_len(Brep)) {
      counts <- rbinom(M, size = n, prob = p)
      unseen <- counts == 0
      Mmax <- if (any(unseen)) max(p[unseen]) else 0
      M1 <- sum(p[unseen])
      M2 <- sum(p[unseen]^2)
      fit <- lower_bound_mmax(n, counts[counts > 0], alpha, details = TRUE)
      runs[b, ] <- c(fit$lower, Mmax, M1, M2, if (M1 > 0) M2 / M1 else 0,
                     fit$M1_upper, fit$M2_lower, fit$M2_hat, fit$M2_lower_raw,
                     fit$bias_upper, fit$variance_upper, fit$tail_radius)
    }
    q <- quantile(runs[, "lower"], c(0.025, 0.5, 0.975), names = FALSE)
    has_unseen <- runs[, "Mmax"] > 0
    covered <- runs[, "lower"] <= runs[, "Mmax"]
    summary <- cbind(task, data.frame(
      Brep = Brep, alpha = alpha, S = sum(p), oracle = oracle,
      lb_q025 = q[1], lb_median = q[2], lb_q975 = q[3],
      coverage = mean(covered), failure_count = sum(!covered),
      coverage_mcse = sqrt(mean(covered) * (1 - mean(covered)) / Brep),
      oracle_coverage = mean(oracle <= runs[, "Mmax"]),
      positive_rate = mean(runs[, "lower"] > 0), unseen_rate = mean(has_unseen),
      mean_Mmax = mean(runs[, "Mmax"]), mean_moment_ratio = mean(runs[, "moment_ratio"]),
      mean_relative_lower = if (any(has_unseen)) {
        mean(runs[has_unseen, "lower"] / runs[has_unseen, "Mmax"])
      } else NA_real_,
      M1_coverage = mean(runs[, "M1"] <= runs[, "M1_upper"]),
      M2_coverage = mean(runs[, "M2_lower"] <= runs[, "M2"])))
    cat(sprintf("Task %d: %s param=%g %s M=%d n=%d coverage=%.4f positive=%.4f\n",
                task$task_id, task$distribution, task$param, task$arm,
                M, n, summary$coverage, summary$positive_rate))
    list(task = task, summary = summary, runs = runs)
  }
  # These functions use only their arguments and functions loaded on workers;
  # avoid serializing main()'s task/stream lists with every function closure.
  environment(oracle_lower) <- .GlobalEnv
  environment(simulate_cell) <- .GlobalEnv

  cat(sprintf("wd=%s\nchoose_wd=%d workers=%d Brep=%d alpha=%g seed=%d\n",
              wd, choose_wd, workers, Brep, alpha, seed))
  if (report_only) {
    cat("Modalita report_only: caricamento del checkpoint, senza simulazioni.\n")
  } else {
    cat(sprintf("%d tasks; %d repetitions; outputs: %s\n",
                nrow(tasks), nrow(tasks) * Brep, lower_dir))
  }
  flush.console()
  started <- proc.time()["elapsed"]
  config <- list(wd = wd, choose_wd = choose_wd, workers = workers,
                 Brep = Brep, alpha = alpha, seed = seed,
                 Mgrid = Mgrid, Ngrid = Ngrid, n_fixed = n_fixed, M_fixed = M_fixed,
                 benchmarks = benchmarks, RNGkind = RNGkind(), session = sessionInfo())

  if (report_only) {
    checkpoint_file <- file.path(lower_dir, "save", "Exp_LB_checkpoint.rds")
    if (!file.exists(checkpoint_file)) stop("Checkpoint not found: ", checkpoint_file)
    checkpoint <- readRDS(checkpoint_file)
    config <- checkpoint$config
    results <- checkpoint$results
    if (!length(results) || !is.data.frame(config$benchmarks)) stop("Invalid checkpoint.")
    benchmarks <- config$benchmarks
    if (is.null(benchmarks$tag)) {
      benchmarks$tag <- vapply(benchmarks$param, get_first3digits, character(1), d = 4)
    }
    tasks <- do.call(rbind, lapply(results, `[[`, "task"))
    alpha <- config$alpha
    Brep <- config$Brep
    seed <- config$seed
    workers <- config$workers
    cat("REPORT ONLY: rebuilding saved results; no simulation or workers started.\n")
  } else {
    cluster <- parallel::makePSOCKcluster(workers, outfile = "")
    on.exit(if (!is.null(cluster)) parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterCall(cluster, function(root) {
      setwd(root)
      source(file.path(root, "LowerBound", "lower_bound.R"), local = .GlobalEnv)
      source(file.path(root, "R", "Rfunctions.R"), local = .GlobalEnv)
      NULL
    }, wd)
    results <- parallel::parLapplyLB(cluster, jobs, simulate_cell,
                                    Brep = Brep, alpha = alpha, oracle_lower = oracle_lower)
    parallel::stopCluster(cluster)
    cluster <- NULL
  }
  summary_results <- do.call(rbind, lapply(results, `[[`, "summary"))
  # Preserve completed simulations even if a subsequent figure cannot be saved.
  if (!report_only) {
    saveRDS(list(config = config, results = results),
            file.path(lower_dir, "save", "Exp_LB_checkpoint.rds"))
  }

  # Add diagnostics to both new and existing checkpoints before exporting.
  for (i in seq_along(results)) {
    runs <- results[[i]]$runs
    seen <- runs[, "Mmax"] > 0
    results[[i]]$summary$coverage_given_unseen <- if (any(seen)) {
      mean(runs[seen, "lower"] <= runs[seen, "Mmax"])
    } else NA_real_
    results[[i]]$summary$positive_rate_given_unseen <- if (any(seen)) {
      mean(runs[seen, "lower"] > 0)
    } else NA_real_
    # Exact probability that all positive-prevalence features have been seen.
    # It can be >= alpha even when some samples still have unseen features.
    task <- results[[i]]$task
    p <- switch(task$distribution,
                Zipfs = (seq_len(task$M) + 1)^(-task$param),
                Geom = (1 - task$param)^seq_len(task$M),
                Constant = rep(task$param, task$M),
                stop("Unsupported benchmark: ", task$distribution))
    p <- p[p > 0]
    log_p0 <- sum(log(-expm1(task$n * log1p(-p))))
    results[[i]]$summary$prob_Mmax_zero <- exp(log_p0)
    results[[i]]$summary$oracle_is_zero <- results[[i]]$summary$oracle == 0
    results[[i]]$summary$oracle_equals_lb_median <-
      results[[i]]$summary$oracle == results[[i]]$summary$lb_median
    covered <- sum(runs[, "lower"] <= runs[, "Mmax"])
    ci <- binom.test(covered, nrow(runs))$conf.int
    results[[i]]$summary$coverage_ci_low <- ci[1]
    results[[i]]$summary$coverage_ci_high <- ci[2]
  }
  summary_results <- do.call(rbind, lapply(results, `[[`, "summary"))

  # Readable, offline HTML tables, plus full-precision CSV exports.
  html_escape <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    gsub(">", "&gt;", x, fixed = TRUE)
  }
  write_results_table <- function(z, title, stem) {
    columns <- c("distribution", "param", "arm", "M", "n", "Brep", "alpha", "oracle",
                 "lb_q025", "lb_median", "lb_q975", "coverage", "coverage_ci_low",
                 "coverage_ci_high", "coverage_mcse", "oracle_coverage", "failure_count",
                 "positive_rate", "unseen_rate", "coverage_given_unseen",
                 "positive_rate_given_unseen", "prob_Mmax_zero", "oracle_is_zero",
                 "oracle_equals_lb_median", "mean_Mmax", "mean_moment_ratio",
                 "mean_relative_lower", "M1_coverage", "M2_coverage")
    table <- z[, columns, drop = FALSE]
    write.csv(table, paste0(stem, ".csv"), row.names = FALSE, na = "NA")
    cells <- lapply(table, function(x) {
      if (is.numeric(x)) {
        vapply(x, function(v) if (is.na(v)) "NA" else format(v, digits = 7, trim = TRUE),
               character(1))
      } else as.character(x)
    })
    rows <- vapply(seq_len(nrow(table)), function(i) {
      paste0("<tr>", paste0("<td>", html_escape(vapply(cells, `[`, character(1), i)),
                            "</td>", collapse = ""), "</tr>")
    }, character(1))
    page <- c('<!DOCTYPE html><html lang="it"><head><meta charset="utf-8">',
      paste0("<title>", html_escape(title), "</title>"),
      '<meta name="viewport" content="width=device-width,initial-scale=1">',
      '<style>body{font:16px/1.5 system-ui,sans-serif;margin:30px;color:#202b36}',
      'table{border-collapse:collapse;white-space:nowrap}th,td{padding:8px 12px;border:1px solid #ddd}',
      'th{background:#eef2f6;position:sticky;top:0}tr:nth-child(even){background:#fafafa}',
      'td{text-align:right}td:first-child{text-align:left}.scroll{overflow:auto;max-height:75vh}',
      'a{color:#125d9a}code{background:#f0f3f6;padding:2px 4px}</style></head><body>',
      paste0("<h1>", html_escape(title), "</h1>"),
      paste0('<p><a href="', basename(stem), '.csv">Scarica CSV</a>. ',
             'Valori di oracle e LB nella scala originale delle probabilita; nei grafici sono moltiplicati per 1000.</p>'),
      '<p>Ogni riga e una coppia (M,n). coverage = frequenza di LB &lt;= Mmax; ',
      'oracle_coverage = frequenza di oracle &lt;= Mmax. Gli intervalli della coverage sono ',
      'intervalli binomiali Monte Carlo al 95%. positive_rate = frequenza di LB &gt; 0; ',
      'unseen_rate = frequenza di Mmax &gt; 0. Le quantita given_unseen sono condizionate ',
      'a Mmax &gt; 0 e valgono NA se non ci sono tali campioni.</p>',
      '<p>prob_Mmax_zero e la probabilita esatta di Mmax = 0; oracle_is_zero segnala ',
      'un oracle nullo, oracle_equals_lb_median segnala curve coincidenti in quel punto. ',
      'M1_coverage e M2_coverage riportano le coperture delle due componenti del bound.</p>',
      '<div class="scroll"><table><thead><tr>',
      paste0("<th>", html_escape(names(table)), "</th>", collapse = ""),
      '</tr></thead><tbody>', rows, '</tbody></table></div></body></html>')
    writeLines(enc2utf8(page), paste0(stem, ".html"), useBytes = TRUE)
  }

  # Original figure style: base-R PDF, grey grid, bty='l', lwd=3,
  # 1000 * bound on the vertical axis; only the median curve is displayed.
  # The 2.5% and 97.5% repetition quantiles are saved as in the original study.
  plot_arm <- function(z, arm, filename) {
    x <- if (arm == "nfix") z$M else z$n
    ord <- order(x)
    x <- x[ord]
    z <- z[ord, , drop = FALSE]
    ymax <- max(c(z$lb_q975, z$oracle), na.rm = TRUE) * 1.05
    if (ymax <= 0) ymax <- 1e-3
    pdf(filename)
    on.exit(dev.off(), add = TRUE)
    par(mfrow = c(1, 1), mar = c(4, 4, 3, 0.5), mgp = c(2.5, 0.5, 0), bty = "l")
    plot(0, 0, type = "n", xaxt = "s", yaxt = "n", xlim = range(x),
         ylim = c(-0.04 * ymax, ymax), xlab = if (arm == "nfix") "M" else "n",
         ylab = "1000 * bound", main = " ")
    grid(lty = 1, lwd = 1, col = "gray90")
    ticks <- pretty(c(0, ymax * 1e3))
    ticks <- ticks[ticks >= 0 & ticks <= ymax * 1e3]
    axis(2, at = ticks * 1e-3, labels = ticks, las = 1, cex.axis = 1)
    lines(x, z$lb_median, lwd = 3, lty = 1, col = "darkorange")
    lines(x, z$oracle, lwd = 3, lty = 1, col = "black")
    points(x, z$oracle, pch = 1, cex = 0.7, col = "black")
    points(x, z$lb_median, pch = 4, cex = 0.7, col = "darkorange")
    note <- if (all(z$oracle == 0) && all(z$lb_median == 0)) {
      "Oracle e mediana LB entrambi nulli"
    } else if (all(z$oracle == 0)) "Oracle nullo su tutta la griglia" else {
      if (all(z$oracle == z$lb_median)) "Oracle e mediana LB coincidenti" else ""
    }
    mtext(sprintf("Oracle: min=%.4g; max=%.4g (scala originale)",
                  min(z$oracle), max(z$oracle)), side = 3, line = 1.1, cex = 0.7)
    if (nzchar(note)) mtext(note, side = 3, line = 0.1, cex = 0.7)
    legend("topright", c("Oracle lower", "Lower bound"),
           lwd = 3, pch = c(1, 4), col = c("black", "darkorange"), bty = "n")
  }

  for (d in seq_len(nrow(benchmarks))) {
    name <- benchmarks$distribution[d]
    param <- benchmarks$param[d]
    tag <- benchmarks$tag[d]
    for (arm in c("nfix", "Mfix")) {
      ids <- which(tasks$distribution == name & tasks$param == param & tasks$arm == arm)
      z <- summary_results[ids, , drop = FALSE]
      x <- if (arm == "nfix") z$M else z$n
      ord <- order(x)
      z <- z[ord, , drop = FALSE]
      selected <- results[ids[ord]]
      exp_name <- paste0("Exp_LB_", arm, "_", name, "_", tag)
      save_res <- list(config = config, summary = z,
                       lower_mat = do.call(rbind, lapply(selected, function(r) r$runs[, "lower"])),
                       oracle = z$oracle,
                       oracle_mat = do.call(rbind, lapply(selected, function(r) r$runs[, "Mmax"])),
                       cells = selected)
      save(save_res, file = file.path(lower_dir, "save", name, paste0(exp_name, ".Rdat")))
      write.csv(z, file.path(lower_dir, "save", name, paste0(exp_name, ".csv")), row.names = FALSE)
      dir.create(file.path(lower_dir, "res", name), recursive = TRUE, showWarnings = FALSE)
      write_results_table(z, sprintf("%s, parametro %g, %s", name, param, arm),
                          file.path(lower_dir, "res", name, exp_name))
      cov_mat <- cbind(if (arm == "nfix") z$M else z$n,
                       Lower = z$coverage, Oracle = z$oracle_coverage)
      colnames(cov_mat)[1] <- if (arm == "nfix") "M" else "n"
      save(cov_mat, file = file.path(lower_dir, "save", "Coverage", name,
                                    paste0(exp_name, ".Rdat")))
      plot_arm(z, arm, file.path(lower_dir, "img", name, paste0(exp_name, ".pdf")))
    }
  }
  write.csv(summary_results, file.path(lower_dir, "save", "Exp_LB_summary.csv"), row.names = FALSE)
  write_results_table(summary_results, "Lower bound: tutti gli scenari",
                      file.path(lower_dir, "res", "Exp_LB_summary"))
  saveRDS(config, file.path(lower_dir, "save", "Exp_LB_config.rds"))
  recap <- aggregate(cbind(coverage, positive_rate, unseen_rate) ~ distribution + param,
                     data = summary_results, FUN = mean)
  cat("\nRIEPILOGO RISULTATI\n")
  cat(sprintf("Scenari: %d | Repliche per scenario: %d | Repliche totali: %d\n",
              nrow(summary_results), Brep, sum(summary_results$Brep)))
  cat(sprintf("alpha=%g | seed=%d | workers=%d | Tempo: %.1f secondi\n",
              alpha, seed, workers, proc.time()["elapsed"] - started))
  cat(sprintf("Copertura empirica: min=%.4f, max=%.4f | Violazioni totali: %d\n",
              min(summary_results$coverage), max(summary_results$coverage),
              sum(summary_results$failure_count)))
  cat("Medie sulle celle di ciascun benchmark (coverage, positive_rate, unseen_rate):\n")
  print(recap, row.names = FALSE, digits = 4)
  cat("Figure PDF salvate: ", 2L * nrow(benchmarks), "\n", sep = "")
  cat("Figure: ", file.path(lower_dir, "img"), "\n", sep = "")
  cat("Risultati: ", file.path(lower_dir, "save"), "\n", sep = "")
  cat("Tabelle HTML/CSV: ", file.path(lower_dir, "res"), "\n", sep = "")
  cat(sprintf("Oracle nullo in %d/%d celle; coincidente con mediana LB in %d/%d celle.\n",
              sum(summary_results$oracle_is_zero), nrow(summary_results),
              sum(summary_results$oracle_equals_lb_median), nrow(summary_results)))
  cat("Summary: ", file.path(lower_dir, "save", "Exp_LB_summary.csv"), "\n", sep = "")
  cat("FINITO\n")
  flush.console()
}

main()
