# Reporting only: no random sampling or simulation is performed here.

# Older checkpoints and new cells may have different diagnostic columns.
lb_bind_rows <- function(rows) {
  columns <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(x) {
    for (name in setdiff(columns, names(x))) x[[name]] <- NA
    x[, columns, drop = FALSE]
  })
  do.call(rbind, rows)
}

lb_task_key <- function(tasks) {
  paste(tasks$distribution, sprintf("%.17g", tasks$param), tasks$arm,
        tasks$M, tasks$n, sep = "|")
}

lb_summary_from_runs <- function(results) {
  lb_bind_rows(lapply(results, function(cell) {
    z <- cell$summary
    runs <- cell$runs
    below_oracle <- runs[, "lower"] <= z$oracle[1]
    below_Mmax <- runs[, "lower"] <= runs[, "Mmax"]
    unseen <- runs[, "Mmax"] > 0
    z$coverage_definition <- "LB <= Mmax"
    z$coverage <- mean(below_Mmax)
    z$coverage_vs_oracle <- mean(below_oracle)
    z$coverage_Mmax <- mean(below_Mmax)
    z$failure_count <- sum(!below_Mmax)
    z$failure_count_Mmax <- sum(!below_Mmax)
    z$coverage_mcse <- sqrt(z$coverage * (1 - z$coverage) / nrow(runs))
    ci <- binom.test(sum(below_Mmax), nrow(runs))$conf.int
    z$coverage_ci_low <- ci[1]
    z$coverage_ci_high <- ci[2]
    z$coverage_given_unseen <- if (any(unseen)) mean(below_Mmax[unseen]) else NA_real_
    z$coverage_Mmax_given_unseen <- if (any(unseen)) mean(below_Mmax[unseen]) else NA_real_
    z
  }))
}

# Normalize legacy CSVs while keeping the oracle comparison separate.
lb_upgrade_summary <- function(z) {
  was_oracle <- "coverage_definition" %in% names(z) &&
    all(z$coverage_definition == "LB <= oracle")
  if (was_oracle) {
    z$coverage_vs_oracle <- z$coverage
    if (!"coverage_Mmax" %in% names(z)) stop("Missing original Mmax coverage; provide checkpoint.")
    z$coverage <- z$coverage_Mmax
    if ("coverage_Mmax_given_unseen" %in% names(z)) {
      z$coverage_given_unseen <- z$coverage_Mmax_given_unseen
    }
  }
  z$coverage_Mmax <- z$coverage
  z$coverage_definition <- "LB <= Mmax"
  known <- !is.na(z$coverage)
  z$failure_count <- round(z$Brep * (1 - z$coverage))
  z$coverage_mcse <- sqrt(z$coverage * (1 - z$coverage) / z$Brep)
  z$coverage_ci_low <- z$coverage_ci_high <- NA_real_
  for (i in which(known)) {
    ci <- binom.test(round(z$Brep[i] * z$coverage[i]), z$Brep[i])$conf.int
    z$coverage_ci_low[i] <- ci[1]
    z$coverage_ci_high[i] <- ci[2]
  }
  z
}

write_lb_results_report <- function(summary, filename) {
  required <- c("distribution", "param", "arm", "M", "n", "Brep", "alpha",
                "oracle", "lb_median", "coverage", "coverage_definition")
  if (!all(required %in% names(summary))) stop("Incomplete report summary.")
  if (any(summary$coverage_definition != "LB <= Mmax")) stop("Incorrect coverage definition.")
  esc <- function(x) {
    x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    gsub(">", "&gt;", x, fixed = TRUE)
  }
  number <- function(x) format(x, digits = 7, scientific = FALSE, trim = TRUE)
  probability <- function(x) {
    if (is.na(x)) return('<span class="missing">Non disponibile</span>')
    format(x, digits = 7, trim = TRUE)
  }
  blocks <- list()
  navigation <- character()
  groups <- unique(summary[, c("distribution", "param", "arm")])
  groups <- groups[order(groups$distribution, groups$param, groups$arm), ]
  for (g in seq_len(nrow(groups))) {
    group <- groups[g, ]
    z <- summary[summary$distribution == group$distribution &
                   summary$param == group$param & summary$arm == group$arm, ]
    varying <- if (group$arm == "nfix") "M" else "n"
    fixed <- if (varying == "M") "n" else "M"
    z <- z[order(z[[varying]]), ]
    if (length(unique(z[[fixed]])) != 1L || length(unique(z$Brep)) != 1L ||
        length(unique(z$alpha)) != 1L) stop("Inconsistent scenario metadata.")
    id <- paste0("scenario-", g)
    caption <- sprintf("distribution: %s | param: %s | arm: %s | %s = %s (fisso) | %s = %s ... %s (variabile) | Brep = %s | alpha = %s",
                       group$distribution, number(group$param), group$arm, fixed,
                       number(z[[fixed]][1]), varying, number(min(z[[varying]])),
                       number(max(z[[varying]])), number(z$Brep[1]), number(z$alpha[1]))
    navigation <- c(navigation, sprintf('<a href="#%s">%s %s &mdash; varia %s</a>',
                                       id, esc(group$distribution), number(group$param), varying))
    rows <- vapply(seq_len(nrow(z)), function(i) {
      paste0('<tbody><tr class="coordinate"><th colspan="3" scope="rowgroup">',
             varying, ' = ', number(z[[varying]][i]), '</th></tr><tr>',
             '<td>', probability(z$oracle[i]), '</td>',
             '<td>', probability(z$lb_median[i]), '</td>',
             '<td>', if (is.na(z$coverage[i])) probability(NA_real_) else sprintf("%.4f", z$coverage[i]),
             '</td></tr></tbody>')
    }, character(1))
    blocks[[g]] <- c(sprintf('<section id="%s"><table>', id),
                    paste0('<caption>', esc(caption), '</caption>'),
                    '<thead><tr><th scope="col"><strong>oracle</strong></th>',
                    '<th scope="col"><strong>lb_median</strong></th><th scope="col">coverage</th></tr></thead>',
                    rows, '</table></section>')
  }
  missing <- sum(is.na(summary$coverage))
  note <- if (missing) sprintf('<aside>Coverage non disponibile per %d righe; fornire il checkpoint con le repliche MC.</aside>', missing) else ""
  page <- c('<!DOCTYPE html><html lang="it"><head><meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width,initial-scale=1">',
    '<title>Lower bound - risultati per scenario</title><style>',
    'body{max-width:900px;margin:auto;padding:32px 24px;color:#202b36;font:16px/1.55 system-ui,sans-serif}',
    'h1{font-size:1.8rem}nav{display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:6px;background:#f4f6f8;padding:16px}',
    'a{color:#125d9a}section{margin:36px 0;scroll-margin-top:16px}table{width:100%;border-collapse:collapse;table-layout:fixed}',
    'caption{text-align:left;font-weight:600;background:#eef2f6;padding:14px;overflow-wrap:anywhere}',
    'th,td{padding:9px 14px;text-align:right;border-bottom:1px solid #dce3e9;font-variant-numeric:tabular-nums}',
    'thead th{font-weight:400;background:#fafafa}.coordinate th{text-align:left;padding:8px 14px 2px;color:#536275;font-size:.85rem;background:#fafafa}',
    '.missing{color:#995200;font-size:.85rem}aside{padding:14px;border-left:3px solid #d98916;background:#fff8ea;margin:18px 0}',
    '@media print{nav{display:none}tbody{break-inside:avoid}section{break-before:page}}</style></head><body>',
    '<h1>Lower bound: risultati per scenario</h1>',
    '<p>Ogni tabella fissa distribuzione e parametro e varia solo M oppure n. ',
    'Le tre colonne sono nella scala originale, senza il fattore 1000 dei grafici.</p>',
    '<p><strong>coverage</strong> = numero di repliche con LB<sub>b</sub> &le; M<sub>max,b</sub>, diviso Brep. ',
    'Si confrontano il LB e il vero M<sub>max</sub> della stessa replica, non le mediane. ',
    'Il confronto con l\'oracle e conservato solo nei dati diagnostici come <code>coverage_vs_oracle</code>.</p>',
    note, '<nav aria-label="Scenari">', navigation, '</nav>', unlist(blocks), '</body></html>')
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  writeLines(enc2utf8(page), filename, useBytes = TRUE)
  invisible(list(tables = nrow(groups), rows = nrow(summary), missing = missing))
}
