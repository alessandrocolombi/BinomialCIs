# Rebuild only the single HTML report, never simulate or start workers.
# Rscript LowerBound/render_results.R
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Run this script with Rscript.")
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE))
source(file.path(script_dir, "results_report.R"))
source(file.path(script_dir, "lower_bound.R"))
checkpoint_file <- file.path(script_dir, "save", "Exp_LB_checkpoint.rds")
csv <- file.path(script_dir, "res", "Exp_LB_summary.csv")
if (file.exists(checkpoint_file)) {
  checkpoint <- readRDS(checkpoint_file)
  summary_results <- lb_summary_from_runs(checkpoint$results)
} else {
  if (!file.exists(csv)) stop("No checkpoint or summary CSV found.")
  summary_results <- lb_upgrade_summary(read.csv(csv, stringsAsFactors = FALSE))
}
html <- file.path(script_dir, "res", "Exp_LB_summary.html")
recap <- write_lb_results_report(summary_results, html)
write.csv(summary_results, csv, row.names = FALSE, na = "NA")
groups <- unique(summary_results[, c("distribution", "param", "arm")])
for (i in seq_len(nrow(groups))) {
  group <- groups[i, ]
  z <- summary_results[summary_results$distribution == group$distribution &
                         summary_results$param == group$param &
                         summary_results$arm == group$arm, ]
  varying <- if (group$arm == "nfix") "M" else "n"
  z <- z[order(z[[varying]]), ]
  tag <- substr(gsub(".", "", format(group$param, scientific = FALSE, trim = TRUE), fixed = TRUE), 1, 4)
  stem <- file.path(script_dir, "res", group$distribution,
                    paste0("Exp_LB_", group$arm, "_", group$distribution, "_", tag))
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  write.csv(z, paste0(stem, ".csv"), row.names = FALSE, na = "NA")
  # Retire only the old generated per-scenario HTML: one consolidated file.
  obsolete <- paste0(stem, ".html")
  if (file.exists(obsolete) && unlink(obsolete) != 0) stop("Cannot remove old report: ", obsolete)
}
cat(sprintf("Report: %d tabelle, %d righe; coverage non disponibile in %d righe.\n",
            recap$tables, recap$rows, recap$missing))
cat("HTML: ", html, "\nFINITO\n", sep = "")
