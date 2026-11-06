# Lower bound for the largest unseen prevalence: theory

The complete formal derivation is now in [theory.tex](theory.tex), with
propositions, theorems and proofs for the oracle hierarchy, the fixed-order
ratio, the joint pairwise threshold, and simultaneous control over orders.
An [HTML reading version](theory.html) is generated from that TeX source.
The earlier derivations below are retained as background notes.

Definitions, concentration bounds and interpretation of the lower endpoint.
For execution and file formats, see the [implementation guide](../README.html).

## Model and target

We observe independent counts
\[
N_j\sim\operatorname{Bin}(n,p_j),\qquad j\ge1,
\]
with fixed prevalences satisfying
\[
S(p)=\sum_{j\ge1}p_j<\infty.
\]
Define
\[
M_{\max}=\max_{j:N_j=0}p_j,\qquad
M_r=\sum_{j:N_j=0}p_j^r,\qquad
\mathcal M_0=M_1.
\]
Use zero for the maximum of an empty set. For r at least one, summability
ensures that the moments are finite. Finite-support prevalence sequences,
extended by zero prevalences, are included in this model.

For n at least four, the construction below yields a nonnegative lower
endpoint satisfying
\[
\inf_{p:S(p)<\infty}\Pr_p\{L_n\le M_{\max}\}\ge1-\alpha.
\]
This is a fixed-sample-size confidence guarantee, not a simultaneous-in-time
or optional-stopping guarantee.

## Why a missing-mass lower bound alone does not suffice

The deterministic inequality is \(M_{\max}\le M_1\). A fixed unseen mass s
can be spread over k unseen categories of prevalence s/k, making the maximum
arbitrarily small as k grows. Thus a lower bound for M1 alone cannot imply
a strictly positive numerical lower bound for Mmax under summability alone.

For \(r>s\ge1\), however,
\[
M_r\le M_{\max}^{r-s}M_s,
\qquad
M_{\max}\ge\left(\frac{M_r}{M_s}\right)^{1/(r-s)}
\]
when Ms is positive. In particular, \(M_{\max}\ge M_2/M_1\).
If M1 is zero, Mmax is zero as well, and the lower endpoint is set to zero.

## Construction and proof

Write
\[
K_{n,k}=\sum_j\mathbf 1\{N_j=k\},\qquad
K_n=\sum_j\mathbf 1\{N_j>0\},\qquad
M_r=\sum_{j:N_j=0}p_j^r,\qquad M_1=\mathcal M_0.
\]
The deterministic inequality is \(M_2\le M_{\max}M_1\).
The lower endpoint is \(L_n=L_2/U_1\), with \(L_2\) a lower bound for \(M_2\)
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
\(a^2\Pr_p(N_j=2)+p_j^4\Pr_p(N_j=0)\); summing and replacing
\(q_j^n\) by \(q_j^{n-4}\) gives the variance bound.
Each centered coordinate is at most \(a\), since its mean is nonnegative.
A one-sided Bernstein bound therefore implies, for \(x=\log(1/\epsilon)\),

\[
\Pr\{D_2>b+\sqrt{2vx}+2ax/3\}\le\epsilon.
\]

These bounds apply to the countable sum: \(S(p)<\infty\) ensures finite count
expectations, finite variance, and convergence of the sums.

Use the expected-count bound from Ayed et al. (2019), Remark A.1:

\[
C_k(\epsilon)=
\left(\sqrt{\tfrac12\log(1/\epsilon)}+
\sqrt{K_{n,k}+\tfrac76\log(1/\epsilon)}\right)^2,
\quad \Pr\{\eta_k\le C_k\}\ge1-\epsilon.
\]

Use upper bounds for the three expected counts to define
\[
\widehat b=\frac{2C_3(\epsilon)}{\binom n3},\qquad
\widehat v=a^2C_2(\epsilon)+\frac{C_4(\epsilon)}{\binom n4}.
\]
Set

\[
L_2=\left[aK_{n,2}-\widehat b-\sqrt{2\widehat v x}-2ax/3\right]_+.
\]

For \(U_1\), use the one-sided construction in Ayed et al.,
Proposition A.3 and equations (A.22)-(A.23). If its error budget is \(\gamma\),
let \(y=\log(2/\gamma)\) and apply the same expected-count bound to \(K_n\),
with error \(\gamma/2\), obtaining
\[
C_K=\left(\sqrt{y/2}+\sqrt{K_n+7y/6}\right)^2.
\]
Then

\[
U_1=\frac{K_{n,1}}n+\frac y{n-1}
+\sqrt{\frac{8yC_K}{(n-1)^2(1-2/n)}}.
\]

The two internal events each have error \(\gamma/2\). This explicit one-sided
allocation differs from directly substituting the simultaneous interval in
Theorem 3.1 and avoids ambiguity about its internal error allocation.

Finally take \(\epsilon=\gamma=\alpha/5\): the three count events, the
Bernstein event, and the \(M_1\) event have total error at most \(\alpha\).
No independence between these events is needed. On their intersection,
\[
M_{\max}\ge\frac{M_2}{M_1}\ge\frac{L_2}{U_1}=L_n
\]
when \(M_1>0\). If \(M_1=0\), then \(M_2=M_{\max}=0\) and the valid
lower bound also vanishes. This proves the stated coverage guarantee.

Reference: Ayed, Battiston, Camerlenghi and Favaro (2019), *A Good-Turing
estimator for feature allocation models*, Electronic Journal of Statistics
13, 3775-3804, DOI: 10.1214/19-EJS1614.
The \(M_2\) construction above is the proposed extension, not a theorem quoted
from that article.


## Deterministic oracle lower endpoint

An **oracle deterministic lower endpoint** is the largest nonrandom lower
endpoint satisfying the required coverage when the prevalence vector is known.
For a finite-support prevalence vector, independence gives

\[
\Pr_p(M_{\max}<t)=\prod_{j:p_j\ge t}\{1-(1-p_j)^n\},\qquad t>0.
\]

The oracle is the largest point in \(\{0\}\cup\{p_j:p_j>0\}\) such that this
strict lower tail is at most \(\alpha\). It is thus the largest deterministic
endpoint with \(\Pr_p(M_{\max}\ge\mathrm{oracle})\ge1-\alpha\), with exact
handling of equality at probability atoms.
This oracle is not the realized \(M_{\max}\) for each sample, which would be
exactly known if both p and that sample were supplied, nor is it the moment
ratio \(M_2/M_1\).
An adaptive lower-bound median can exceed this deterministic oracle without
contradiction; coverage is a statement about its joint behavior with `Mmax`.

The oracle lower endpoint can genuinely be zero. In particular, if
\(\Pr_p(M_{\max}=0)>\alpha\), every strictly positive deterministic endpoint
fails the required lower-coverage guarantee. This does not require
\(M_{\max}=0\) in every repetition.


## Interpretation and loss from the moment ratio

The ratio \(R_n=M_2/M_1\) is a prevalence-weighted mean of the unseen
prevalences, rather than their maximum. It equals Mmax when all positive
unseen prevalences are equal, but can be substantially smaller otherwise.
When M2 is positive,
\[
\frac{L_n}{M_{\max}}=
\frac{L_2}{M_2}\frac{M_1}{U_1}\frac{M_2}{M_1M_{\max}}.
\]
This separates the statistical losses in the two confidence bounds from
the intrinsic loss of using the second-to-first moment ratio.

When \(M_1>0\), the deterministic lower bounds
\[
B_r=\left(\frac{M_r}{M_1}\right)^{1/(r-1)},\qquad r>1,
\]
are nondecreasing in r and converge to Mmax. Indeed, with
\(w_j=p_j/M_1\) on the unseen categories, they are weighted power means:
\[
B_r=\left(\sum_{j:N_j=0}w_jp_j^{r-1}\right)^{1/(r-1)},
\qquad\sum_jw_j=1.
\]
Optimizing estimated lower bounds over r requires simultaneous confidence
control. Monotonicity of the true moments does not ensure that estimated
endpoints improve with r.

## Relation to the lower bounds in the submitted manuscript

Theorems 3 and 6 in *Binomial_CIs_jasa_sub* lower-bound the attainable size
of any uniformly valid upper confidence bound Tn along least-favourable
sequences of prevalence vectors. They are minimax limitations on procedures,
not lower confidence bounds for Mmax for every prevalence vector.

The new construction instead bounds Mmax from below, for every fixed
summable prevalence sequence and at finite n. It complements an upper
confidence bound, rather than strengthening the minimax theorem.
If an upper endpoint Tn has error delta/2 and the lower endpoint has error
delta/2, the union bound gives
\[
\Pr_p\{L_n(\delta/2)\le M_{\max}\le T_n(\delta/2)\}\ge1-\delta.
\]
A lower bound on the upper endpoint Tn is not automatically a lower bound
on the length Tn minus Ln of this two-sided interval. No optimality claim
for the new lower endpoint or interval length has been established.

## Two-sided confidence bounds for an unseen moment of fixed order

Fix an integer \(r\ge2\), \(n\ge2r\), and \(\delta\in(0,1)\).
The following construction is valid for every fixed summable prevalence
sequence; no regular variation or knowledge of S(p) is needed.
It is an explicit conservative interval, not an optimality claim.

Set
\[
a_r=\binom nr^{-1},\qquad
T_r=a_rK_{n,r},\qquad D_r=T_r-M_r,\qquad
\eta_k=\mathbb E_pK_{n,k}.
\]
Independence across categories implies
\[
0\le\mathbb E_pD_r\le b_r:=\frac{r\eta_{r+1}}{\binom n{r+1}},\qquad
\operatorname{Var}_p(D_r)\le v_r:=a_r^2\eta_r+
\frac{\eta_{2r}}{\binom n{2r}}.
\]

Indeed, writing \(q_j=1-p_j\), the coordinate mean is
\[
\mathbb E_pd_{r,j}=p_j^rq_j^{n-r}(1-q_j^r)
\le r p_j^{r+1}q_j^{n-r-1},\qquad
d_{r,j}=a_r\mathbf1\{N_j=r\}-p_j^r\mathbf1\{N_j=0\}.
\]
For the variance, use disjointness of the two indicator events:
\[
\operatorname{Var}(D_r)\le\sum_j\mathbb E[d_{r,j}^2]
=a_r^2\eta_r+\sum_jp_j^{2r}q_j^n
\le a_r^2\eta_r+\frac{\eta_{2r}}{\binom n{2r}}.
\]
Since \(d_{r,j}\le a_r\) and its mean is nonnegative, each centered
coordinate is at most \(a_r\). One-sided Bernstein yields
\[
\Pr_p\{D_r>b_r+\sqrt{2v_rx}+2a_rx/3\}\le e^{-x}.
\]
For the other tail, nonnegativity of the mean and Cantelli's inequality give
\[
\Pr_p\left\{D_r<-\sqrt{\frac{1-\epsilon}{\epsilon}v_r}\right\}
\le\epsilon.
\]
If the variance is zero, this conclusion follows directly since Dr is a
nonnegative constant. Summability ensures convergence of the countable sums
and finite variance, so the inequalities apply to the infinite alphabet.

To eliminate the unknown expected counts, let
\[
\epsilon=\delta/5,\qquad x=\log(5/\delta),\qquad
C_k=\left(\sqrt{x/2}+\sqrt{K_{n,k}+7x/6}\right)^2,
\quad k\in\{r,r+1,2r\}.
\]
Each count bound \(\eta_k\le C_k\) has error at most epsilon, as in
Ayed et al. (2019), Remark A.1. Define
\[
\widehat b_r=\frac{rC_{r+1}}{\binom n{r+1}},\qquad
\widehat v_r=a_r^2C_r+\frac{C_{2r}}{\binom n{2r}}.
\]
Fully observable endpoints are
\[
\boxed{L_r=\left[T_r-\widehat b_r-\sqrt{2\widehat v_rx}
-\frac{2a_rx}{3}\right]_+},\qquad
\boxed{U_r=T_r+\sqrt{\frac{1-\epsilon}{\epsilon}\widehat v_r}}.
\]
The three count events, the Bernstein event and the Cantelli event have total
failure probability at most \(5\epsilon=\delta\). On their intersection,
\[
\boxed{L_r\le M_r\le U_r}.
\]
Hence \(\inf_{p:S(p)<\infty}\Pr_p(L_r\le M_r\le U_r)\ge1-\delta\).
The argument does not require independence among those five events.
The upper radius has a square-root inverse-error dependence because it uses
Cantelli; it is not asserted to be the sharpest upper bound available.

If n is smaller than 2r, a conservative fallback valid for any positive n is
\([0,\overline S]\), since \(M_r\le S(p)\) for r at least one. With
\(\widehat S=n^{-1}\sum_jN_j\) and \(z=\log(1/\delta)\), the left-tail
concentration of this sum gives
\[
\overline S=\left(\sqrt{\widehat S+z/(2n)}+\sqrt{z/(2n)}\right)^2,
\qquad\Pr_p\{S(p)\le\overline S\}\ge1-\delta.
\]

For noninteger \(r\ge2\), let \(k=\lfloor r\rfloor\). Since
\(M_{k+1}\le M_r\le M_k\), construct the above intervals for k and k+1,
each with error delta/2. Then \([L_{k+1},U_k]\) covers Mr with probability
at least one minus delta, provided \(n\ge2(k+1)\).

Non-strict inequalities are necessary for a nonnegative lower endpoint when
Mr can be zero. Strict inequalities can be obtained by replacing Lr by
Lr minus any positive constant and Ur by Ur plus a positive constant, at
the cost of allowing a negative lower endpoint.

## Monte Carlo coverage and oracle comparison

Empirical coverage of the random lower endpoint is
\[
\frac1{B_{\mathrm{rep}}}\sum_{b=1}^{B_{\mathrm{rep}}}
\mathbf 1\{L_{n,b}\le M_{\max,b}\}.
\]
Both quantities use the same Monte Carlo sample. This differs from the
fraction of lower endpoints below a deterministic oracle, and neither
comparison can be replaced by comparing their medians. High coverage alone
does not imply a tight lower endpoint: an identically zero endpoint has
coverage one.
