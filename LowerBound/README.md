# LowerBound: implementation and execution

This guide covers the R functions, simulation settings, VM commands and generated files. Mathematical definitions, derivations and interpretation are in [Notes/README.html](Notes/README.html).

## Use

```r
source("LowerBound/lower_bound.R")
L <- lower_bound_mmax(n, counts, alpha = 0.05)
fit <- lower_bound_mmax(n, counts, alpha = 0.05, details = TRUE)
```

`counts` can contain all category counts or just observed positive counts.
An empty numeric vector is allowed. The sample size must satisfy `n >= 4`.
Only the r=2 construction is implemented. Definitions, assumptions and the confidence guarantee are in the [theory notes](Notes/README.html).

## Simulation

```powershell
& 'C:/Program Files/R/R-4.4.2/bin/Rscript.exe' LowerBound/simulate_lower_bound.R
# Optional: repetitions, output directory, seed
& 'C:/Program Files/R/R-4.4.2/bin/Rscript.exe' LowerBound/simulate_lower_bound.R 1000 LowerBound/results 42
```

The script also works with `source("LowerBound/simulate_lower_bound.R")`.
It resolves source paths without changing the working directory and uses no
packages beyond base R. Relative output paths are relative to the caller's
working directory; the default output directory is next to the script.

Four `(M,n)` pairs are selected from actual arms of
`Rscripts/SimulationStudy_firstSub/Exp1_varying_M_n.R`:

| M | n | Original arm |
|---:|---:|---|
| 2000 | 2000 | n fixed |
| 10000 | 2000 | n fixed |
| 5000 | 1000 | M fixed |
| 5000 | 5000 | M fixed |

The nine benchmarks use the existing `sim_generic()` generators and the
parameters in `Exp1_varying_M_n_parallel_main.R`:

- Zipf: `p_j=(j+1)^(-s)`, `s` in `0.25, 0.5, 1.02`.
- Geometric: `p_j=(1-a)^j`, `a` in `0.005, 0.1, 0.25`.
- Constant: `p_j=c`, `c` in `0.001, 0.05, 0.5`.

In each cell prevalences are fixed across repetitions. Sampling is directly
binomial, without constructing an `n x M` incidence matrix. The estimator
receives only observed positive counts. The model interpretation of these
finite-support benchmarks is described in the theory notes.

Defaults: 1000 repetitions per cell, 36 cells, `alpha=0.05`, seed 42.
Outputs:

- `summary.csv`: empirical coverage, exact 95% binomial Monte Carlo intervals,
  failure counts, positive-bound rate, frequency of `Mmax>0`, conditional
  coverage and positive-bound rate, bound quantiles, tightness, and component
  coverage for `M_2` and `M_1`.
- `simulation.rds`: configuration, RNG/session information, summary, and
  per-repetition results including counts and bias/variance/tail diagnostics.
- `diagnostics.png`: coverage, nonzero-bound frequency, and conditional
  tightness, with the true `M_2/M_1` as a deterministic comparison.

Conditional tightness is averaged only over samples with `Mmax>0` and is `NA`
when none occur. High empirical coverage alone does not establish a useful
bound: some benchmarks have no unseen categories, and others yield zero
endpoints despite positive `Mmax`. The formal coverage guarantee is
unconditional; conditional coverage is reported only as a diagnostic.

## Verification

```powershell
& 'C:/Program Files/R/R-4.4.2/bin/Rscript.exe' LowerBound/check_lower_bound.R
```

Checks include a positive-endpoint formula fixture, invariance to unobserved
zero entries, empty observations, invalid inputs, and exact enumeration of
9261 samples from a heterogeneous model to check the bias, variance,
deterministic moment inequality, and Bernstein tail.

## Parallel VM experiment: curves over M and n

Prepared in `Exp_LB_M_n.R`; this experiment has not been executed locally.
From the repository root on the VM:

```bash
./run_job.sh LowerBound/Exp_LB_M_n.R choose_wd=2
```

`run_job.sh` now forwards all arguments after the R script name, preserving
argument boundaries. The VM profile (`choose_wd=2`) sets the working directory
to the root of the checkout containing the launched script. It does not assume
a particular VM username or home directory. `choose_wd=1` selects the local
`C:/Users/colom/BinomialCIs` profile. An explicit repository root is also allowed:

```bash
./run_job.sh LowerBound/Exp_LB_M_n.R choose_wd=2 wd=/home/colombi/BinomialCIs workers=15 Brep=1000
```

For paths containing spaces, quote the entire `wd=...` argument. Other optional
arguments are `seed=42` and `alpha=0.05`. The launcher retains detached execution,
logging in `logs/`, PID recording, and numeric-library thread limits. The script
uses a PSOCK cluster from base `parallel`, defaulting to one fewer physical core
with at most 15 workers. Independent L'Ecuyer streams are assigned per task, so
results are reproducible across worker counts and scheduling orders.

`workers=15` means up to 15 R processes running concurrently. Each worker
handles one simulation cell (one distribution, parameter, and pair \((M,n)\))
at a time, running that cell's `Brep` repetitions sequentially. When a worker
finishes, it receives the next cell. The simulation cells are distributed across these
workers; this does not launch 15 copies of the entire experiment.

The experiment restores the grids from the original `Exp1_varying_M_n.R`:

- Fixed `n=2000`: `M=2000,2250,...,10000`.
- Fixed `M=5000`: `n=1000,1250,...,5000`.

It uses the initial nine benchmarks plus Zipf s=1.2 and s=2, `alpha=0.05`, seed 42, and the preceding
lower-bound simulation's `Brep=1000`. The full current experiment has 550 cells / 550,000 repetitions,
not a Cartesian product of the two grids. Parallelism is across cells; each
worker generates one count vector per repetition, without nested clusters or
large incidence matrices. Results are written after the parallel tasks finish.
Concurrent launches targeting the same checkout overwrite the same result
filenames; use separate checkouts for separate simultaneous experiments.

For each distribution and parameter, two figures are saved, 22 in total:

```text
LowerBound/img/Zipfs/Exp_LB_nfix_Zipfs_*.pdf
LowerBound/img/Zipfs/Exp_LB_Mfix_Zipfs_*.pdf
LowerBound/img/Geom/...
LowerBound/img/Constant/...
```

They use the original base-R figure style: grey grid, `bty="l"`, line width 3,
and vertical axis `1000 * bound`. The darkorange curve is the median confidence
lower endpoint over repetitions; the 2.5% and 97.5% quantiles are saved in the
data, as in the original experiment.

The black curve is the exact deterministic oracle lower endpoint. It is calculated without Monte Carlo oracle repetitions; no Bor parameter is needed. Definitions and interpretation are in the [theory notes](Notes/README.html#deterministic-oracle-lower-endpoint).

Results are organized analogously to the original study:

- `LowerBound/save/<distribution>/Exp_LB_<arm>_<distribution>_<tag>.Rdat`:
  `save_res` containing the summary, lower-bound matrix, exact oracle vector,
  realized-`Mmax` matrix (`oracle_mat`), and all per-cell diagnostics.
- Matching CSV files in the same directory.
- `LowerBound/save/Coverage/<distribution>/...Rdat`: `cov_mat` with empirical
  lower-bound coverage and oracle coverage.
- `LowerBound/save/Exp_LB_summary.csv` and `Exp_LB_config.rds`: aggregate
  summary and run configuration/session information.

### Output names and completion log

The parameter-name helper `get_first3digits()` is defined locally in
`Exp_LB_M_n.R`, with the same formatting as the user's original function.
The master process prepares and checks all parameter tags before starting
the simulations; it does not depend on functions loaded only inside workers.

After the parallel computations finish, workers are stopped and all raw cell
results are saved in `LowerBound/save/Exp_LB_checkpoint.rds`, before generating
figures. This checkpoint preserves the simulations if a later output step fails.
It contains `config` and `results`; use `resume=TRUE` to add missing cells, or `report_only=TRUE` to regenerate outputs from existing cells.

On successful completion the log prints a brief recap: scenario and repetition
counts, settings, elapsed time, minimum and maximum empirical coverage, total
coverage violations, mean coverage/nonzero-bound/unseen rates by benchmark,
and output locations. Its final line is exactly:

```text
FINITO
```

`FINITO` is printed only after all outputs are saved and worker cleanup succeeds;
it is not printed on errors. The current fixes have not been run as a simulation
locally.

The earlier four-pair script `simulate_lower_bound.R` is retained unchanged.

### Oracle visibility and readable result tables

The vertical axis already includes the maximum oracle value, so the oracle
is not cut off above the plot. In the earlier drawing order, an orange zero
LB curve could cover the black zero oracle curve. Plots now distinguish them
with black circles for the oracle and darkorange crosses for the LB. A little
space is left below zero, and each plot states the minimum and maximum oracle
in the original probability scale. A note explicitly flags an entirely zero
oracle or two coincident curves.

The single readable report is `LowerBound/res/Exp_LB_summary.html`. It contains
22 separate tables after the full incremental run, one per distribution, parameter and arm, with only three
base columns **oracle**, **lb_median**, and coverage, plus the new
`lb_joint_median (R=6)` and `coverage_joint`. All metadata (distribution,
parameter, arm, fixed M or n, range of the varying quantity, Brep and alpha)
appear in the caption. Each row group is labelled with its varying M or n value.
The HTML is self-contained and works offline.

The coverage column is computed as `mean(lower <= Mmax)` over individual repetitions, using the true Mmax from the same repetition. `coverage_Mmax` is an alias in the CSV; `coverage_vs_oracle` is the separate diagnostic `mean(lower <= oracle)`. See the theory notes for their interpretation.

Full numerical diagnostics remain in `res/Exp_LB_summary.csv` and the
distribution-specific CSVs: quantiles, Monte Carlo confidence intervals,
violation counts, positive-bound and unseen rates, component coverage, and
tightness. Values are in the original probability scale; plotted bounds are
multiplied by 1000. Old per-scenario HTML files are retired on report rebuild.

To rebuild just this HTML and its CSVs from the saved MC repetitions, without
generating figures, starting workers or drawing any new samples, run:

```bash
Rscript LowerBound/render_results.R
```

The reporting-only script uses the checkpoint when available, otherwise the
existing CSV summary. It prints a short recap followed by `FINITO`.

To regenerate figures and tables from an existing checkpoint without running
the experiment again, use:

```bash
./run_job.sh LowerBound/Exp_LB_M_n.R choose_wd=2 report_only=TRUE
```

This reads `LowerBound/save/Exp_LB_checkpoint.rds`; it starts no workers and
generates no new samples. The original saved configuration determines alpha,
seed, benchmarks and repetition count. Do not pass overrides for those settings
in report-only mode. If the checkpoint is absent, the script stops with a clear
error rather than starting a simulation. The output directory is the selected
checkout. A recap and `FINITO` are printed after all outputs are written.

### Add Zipf exponents 1.2 and 2 without repeating completed cells

The default benchmark list now also includes Zipf \(s=1.2\) and \(s=2\).
The original nine benchmarks and the new two use exactly the same M/n grids,
sample counts, alpha and seed. The complete experiment has 550 cells and
22 figures/tables (11 distribution/parameter combinations, two arms each).

After updating the code on the VM, run from the repository root:

```bash
git pull
./run_job.sh LowerBound/Exp_LB_M_n.R choose_wd=2 resume=TRUE workers=15
```

`resume=TRUE` requires `LowerBound/save/Exp_LB_checkpoint.rds`. It matches saved
cells by distribution, exact parameter, arm, M and n, preserving their samples.
It checks that Brep, alpha, seed and both grids agree with the saved run. With
the original 450-cell checkpoint it simulates only the 100 missing cells
(100,000 new repetitions with Brep=1000). A repeated incremental launch after
those cells are saved simulates no cells.

Before replacing the checkpoint, the script copies the previous file to
`save/Exp_LB_checkpoint_before_resume_<timestamp>_<pid>.rds`. New cells are
appended after the original task list so the original random-stream assignments
remain stable. Parallel scheduling and worker count do not change the streams.
Checkpoints from earlier script versions may have fewer diagnostic columns;
these are normalized when old and new results are combined.

The same command automatically rebuilds all PDFs in `img/`, CSV results and
coverage diagnostics in `save/` and `res/`, and the single readable
`res/Exp_LB_summary.html`. The recap reports reused and simulated cell counts;
wait for the final `FINITO` line in the log indicated by `run_job.sh`.

Defaults remain Brep=1000, alpha=0.05, seed=42. If the previous run used
different values, provide the same saved values on the incremental command.
`resume=TRUE` and `report_only=TRUE` cannot be combined. Omitting both flags
runs the full current 550-cell experiment; using `report_only=TRUE` refreshes
only cells already present in the checkpoint and does not add new benchmarks.
The incremental simulation has been prepared, not executed locally.

## Updating the documentation

From the repository root in PowerShell:

```powershell
./LowerBound/render_readme.ps1
```

This rebuilds both `LowerBound/README.html` (implementation) and
`LowerBound/Notes/README.html` (theory), with embedded styles and offline
formula rendering. Markdown files remain the editable sources. No simulation
or reporting job is started by this command.

## Simultaneous lower bound with R=6

`Exp_LB_M_n.R` now computes both the original `lower_bound_mmax()` and
`lower_bound_mmax_simultaneous(..., R=6)`, at the same per-method error level
`alpha`. The new method uses `K1,...,K12`, chooses among orders 2 through 6
inside a single calibrated event, and saves its per-order thresholds and
selected order. It is not combined by a further maximum with the old method.

The single results HTML now has the original oracle, `lb_median`, `coverage`,
plus `lb_joint_median (R=6)` and `coverage_joint`. Each coverage compares the
method's endpoint against the true Mmax in the same repetition. CSV exports
also contain new-bound quantiles, positive rates, failure counts and tightness.
The figures add the new median as a dashed darkgreen curve, keeping the old
darkorange bound and the black oracle.

After updating code on the VM, launch:

```bash
./run_job.sh LowerBound/Exp_LB_M_n.R choose_wd=2 resume=TRUE workers=15
```

Older checkpoints saved neither the new endpoint nor the first 12 spectrum
counts. Such cells cannot be enriched from their aggregate summaries. The
updated resume mode replays their original task-specific RNG streams, computing
both methods on the same MC samples. It verifies that the old lower endpoint,
Mmax, M1 and M2 match the saved repetitions exactly before replacing any
checkpoint. A mismatch stops the run. The previous checkpoint is backed up.
New benchmark cells, if any, are simulated as usual. Future checkpoints store
the spectrum and can be enriched without regenerating samples.

`report_only=TRUE` and `render_results.R` never draw samples. On old cells
without spectra, they show the new result as unavailable rather than copying
or guessing it. Once the VM replay finishes, figures, CSVs and HTML are updated
automatically and the log ends with a recap and `FINITO`.

Deterministic verification (including exact enumeration, not MC sampling):

```powershell
& 'C:/Program Files/R/R-4.4.2/bin/Rscript.exe' LowerBound/check_simultaneous.R
```

The formal TeX source is `Notes/theory.tex`; the documentation builder also
creates its offline HTML reading version `Notes/theory.html`.
