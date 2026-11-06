# Lower bound for the largest unseen prevalence

Base-R implementation for independent \(N_j\sim\operatorname{Bin}(n,p_j)\) with fixed
prevalences and \(S(p)=\sum_jp_j<\infty\). No knowledge of the alphabet size,
the prevalences, or `S(p)` is supplied to the estimator.

## Use

```r
source("LowerBound/lower_bound.R")
L <- lower_bound_mmax(n, counts, alpha = 0.05)
fit <- lower_bound_mmax(n, counts, alpha = 0.05, details = TRUE)
```

`counts` can contain all category counts or just observed positive counts.
An empty numeric vector is allowed. The sample size must satisfy `n >= 4`.
The guarantee is pointwise in the fixed sampling time:

\[
\Pr_p\{M_{\max}\ge L\}\ge1-\alpha.
\]

This is not a simultaneous guarantee over sample sizes or an optional-stopping
guarantee. Only the proposed `r = 2` bound is implemented; there is no selection
over several moment orders. The moment-matrix extension is not implemented.

## Construction and proof

Write
\[
K_{n,k}=\sum_j\mathbf 1\{N_j=k\},\qquad
K_n=\sum_j\mathbf 1\{N_j>0\},\qquad
M_r=\sum_{j:N_j=0}p_j^r,\qquad M_1=\mathcal M_0.
\]
The deterministic inequality is \(M_2\le M_{\max}M_1\).
The output is \(L=L_2/U_1\), with \(L_2\) a lower bound for \(M_2\)
and \(U_1\) an upper bound for \(M_1\); use zero when the denominator is zero.

Let \(a=\binom n2^{-1}\), \(\eta_k=\mathbb E_p[K_{n,k}]\), and
\(D_2=aK_{n,2}-M_2\). Independence across categories gives

\[
0\le ED_2\le b=\frac{2\eta_3}{\binom n3},\qquad
\operatorname{Var}(D_2)\le v=a^2\eta_2+\frac{\eta_4}{\binom n4}.
\]

Indeed, setting \(q_j=1-p_j\), the coordinate mean is
\(p_j^2q_j^{n-2}(1-q_j^2)\), bounded by \(2p_j^3q_j^{n-3}\).
The coordinate second moment is
`a^2 Pr(N_j=2) + p_j^4 Pr(N_j=0)`; summing and replacing
`q_j^n` by `q_j^(n-4)` gives the variance bound.
Each centered coordinate is at most `a`, since its mean is nonnegative.
A one-sided Bernstein bound therefore implies, for `x=log(1/epsilon)`,

\[
\Pr\{D_2>b+\sqrt{2vx}+2ax/3\}\le\epsilon.
\]

These bounds apply to the countable sum: `S(p)<Inf` ensures finite count
expectations, finite variance, and convergence of the sums.

Use the expected-count bound from Ayed et al. (2019), Remark A.1:

\[
C_k(\epsilon)=
\left(\sqrt{\tfrac12\log(1/\epsilon)}+
\sqrt{K_{n,k}+\tfrac76\log(1/\epsilon)}\right)^2,
\quad \Pr\{\eta_k\le C_k\}\ge1-\epsilon.
\]

Replace `eta_2, eta_3, eta_4` by `C_2, C_3, C_4` in `b,v`, obtaining
`bhat,vhat`. Set

\[
L_2=\left[aK_{n,2}-\widehat b-\sqrt{2\widehat v x}-2ax/3\right]_+.
\]

For `U_1`, the code explicitly uses the one-sided construction in Ayed et al.,
Proposition A.3 and equations (A.22)-(A.23). If its error budget is `gamma`,
let `y=log(2/gamma)` and apply the same expected-count bound to `K_n`,
with error `gamma/2`, obtaining `C_K`. Then

\[
U_1=\frac{K_{n,1}}n+\frac y{n-1}
+\sqrt{\frac{8yC_K}{(n-1)^2(1-2/n)}}.
\]

The two internal events each have error `gamma/2`. This explicit one-sided
allocation differs from directly substituting the simultaneous interval in
Theorem 3.1 and avoids ambiguity about its internal error allocation.

Finally take `epsilon = gamma = alpha/5`: the three count events, the
Bernstein event, and the `M_1` event have total error at most `alpha`.
No independence between these events is needed. The implementation returns
the formula without clipping the endpoint at one.

Reference: Ayed, Battiston, Camerlenghi and Favaro (2019), *A Good-Turing
estimator for feature allocation models*, Electronic Journal of Statistics
13, 3775-3804, DOI: 10.1214/19-EJS1614.
The `M_2` construction above is the proposed extension, not a theorem quoted
from that article.

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
binomial, without constructing an `n x M` incidence matrix. Finite-support
benchmarks are exact members of the infinite-alphabet model after appending
zero prevalences; this is not a simulation of infinitely many positive
prevalences. The estimator receives only observed positive counts.

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
finishes, it receives the next cell. The 450 cells are distributed across these
workers; this does not launch 15 copies of the entire experiment.

The experiment restores the grids from the original `Exp1_varying_M_n.R`:

- Fixed `n=2000`: `M=2000,2250,...,10000`.
- Fixed `M=5000`: `n=1000,1250,...,5000`.

It retains the nine benchmarks above, `alpha=0.05`, seed 42, and the preceding
lower-bound simulation's `Brep=1000`. This is 450 cells / 450,000 repetitions,
not a Cartesian product of the two grids. Parallelism is across cells; each
worker generates one count vector per repetition, without nested clusters or
large incidence matrices. Results are written after the parallel tasks finish.
Concurrent launches targeting the same checkout overwrite the same result
filenames; use separate checkouts for separate simultaneous experiments.

For each distribution and parameter, two figures are saved, 18 in total:

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

The black curve is an **oracle deterministic lower endpoint**, analogous to
the deterministic upper oracle in the original study, with the direction of
the tail reversed. Given the entire prevalence vector, independence gives

\[
\Pr_p(M_{\max}<t)=\prod_{j:p_j\ge t}\{1-(1-p_j)^n\},\qquad t>0.
\]

The oracle is the largest point in `{0} union {p_j:p_j>0}` such that this strict
lower tail is at most `alpha`. It is thus the largest deterministic endpoint
with `Pr_p(Mmax >= oracle) >= 1-alpha` (the discrete lower-alpha quantile,
with exact atom-boundary handling). Products are computed on the log scale.
No Monte Carlo oracle sample or `Bor` parameter is needed.
This oracle is not the realized `Mmax` for each sample, which would be exactly
known if both `p` and that sample were supplied, nor is it the moment ratio
`M_2/M_1`. The per-repetition true `Mmax` and moment ratio are saved separately.
An adaptive lower-bound median can exceed this deterministic oracle without
contradiction; coverage is a statement about its joint behavior with `Mmax`.

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
It contains `config` and `results`; it is not an automatic resume mechanism.

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

The oracle lower endpoint can genuinely be zero. In particular, if
\(\Pr_p(M_{\max}=0)>\alpha\), every strictly positive deterministic endpoint
fails the required lower-coverage guarantee. This does not require
\(M_{\max}=0\) in every repetition.

The experiment also creates `LowerBound/res/`, subdivided into `Zipfs`, `Geom`
and `Constant`. For each figure there is a matching HTML table and CSV file,
with one row per \((M,n)\). Thus the 18 figures have 18 corresponding tables;
the aggregate `res/Exp_LB_summary.html` and `.csv` contain all 450 rows.
HTML files are readable offline, and CSV files retain numerical precision.

Tables contain the oracle, LB median and 2.5%/97.5% repetition quantiles,
empirical coverage, exact binomial 95% Monte Carlo intervals, coverage standard
error, violation count, oracle coverage, positive-bound and unseen rates,
conditional coverage and positive-bound rate given \(M_{\max}>0\), component
coverage for \(M_1,M_2\), and tightness diagnostics. Values are in the original
probability scale; plotted bounds are multiplied by 1000.
Additional columns identify zero oracles, oracle/LB median coincidences, and
the exact probability of \(M_{\max}=0\), so a hidden line is unambiguous.

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
