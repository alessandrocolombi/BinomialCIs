# Proofreading report — JASA draft

Source reviewed: `pasted-text.txt` (2,536 lines of active LaTeX plus comments).

The line numbers below refer to the pasted source. I have separated mathematical issues from language/style corrections. Items marked **critical** can affect the validity or meaning of a stated result; they are not merely editorial.

## 1. Mathematical and logical issues

### Definition and global notation

1. **Undefined maximum over the empty set (lines 324–327).** If every category has been observed, the set \(\{j:N_j=0\}\) is empty. Explicitly define
   \[
   M_{\max}(N)=\max_{j:N_j=0}p_j,\qquad \max\varnothing:=0.
   \]
   This convention is used implicitly later, especially in the lower-bound proofs.

2. **Inconsistent notation for the target.** The draft alternates among `M_{\max}`, `M_{max}`, and `M_{\text{max}}`, and sometimes writes \(M_{\max}(N,p)\). Standardize throughout, preferably with the existing `\Mmax` macro and an explicit definition of any dependence on \(p\).

3. **“Confidence interval” versus “upper confidence bound.”** Most displayed results produce \([0,T_n(N)]\), so “upper confidence bound” is generally more precise. Reserve “confidence interval length” for the endpoint only after stating that the lower endpoint is zero.

### Bounded-alphabet results

4. **Theorem 1 needs \(r^\star\ge1\) (lines 395–405).** The proof requires \(r\ge1\), but \(r^\star=\log(M/\alpha)\) can be below 1 (e.g. \(M=1\) and \(\alpha>e^{-1}\)). Either restrict the parameters or set \(r^\star=\max\{1,\log(M/\alpha)\}\).

5. **The interpretation after Theorem 1 is stronger than the calculation (lines 412–424).** \(T_{\mathrm{hom}}(r)\) is the exact homogeneous/worst-coordinate expression, whereas \(T_{\mathrm{wc}}\) uses the exponential relaxation
   \((n/(n+r))^{n/r}\le \exp\{-n/(n+r)\}\). Therefore it does not “yield precisely” the same bound, and the worst-case bound does not exactly “correspond to constant \(p_i\).” Say that the two are asymptotically equivalent under \(r=o(n)\).

6. **The remainder \(O(n^{-2})\) is not justified when \(M=M_n\) grows (lines 427–435 and 1207–1211).** The remainder depends on \(r^\star=\log(M/\alpha)\). A fixed \(O(n^{-2})\) statement is appropriate only when \(M,\alpha\) are fixed. With growing \(M\), give the remainder explicitly in terms of \(\log(M/\alpha)\), or state only \(T_{\mathrm{wc}}\sim\log(M/\alpha)/n\).

7. **Typo inside the maximization (line 444).** `\max_{p\in[0,1]} p_j^r(1-p_j)^n` should use the maximization variable: `\max_{x\in[0,1]}x^r(1-x)^n`.

8. **Theorem 2 is presented as nonasymptotic but its endpoint contains an unspecified \(O(n^{-2})\) term (lines 461–484).** An unspecified remainder cannot define an implementable finite-sample confidence bound and does not indicate its sign. State the exact endpoint obtained from line 1219 (with the data-dependent upper bound substituted), then give the logarithmic expression as a corollary/asymptotic expansion.

9. **Theorem 2 requires additional domain conditions.** At minimum, \(n>b_n\). Also the optimizing value \(r^\star=\log\{m_b(p)/\alpha\}\) used in the proof need not be at least 1. If the log argument is at most \(e\), the claimed optimizer is outside the allowed domain. Use the constrained optimizer or state sufficient conditions. A negative logarithmic endpoint is impossible for a nonnegative target unless it is truncated/handled separately.

10. **Inconsistent \(\epsilon\) in the special case (lines 468–483).** The general definition contains \(\log(1/\delta)\), whereas the \(b_n=\log n\) version suddenly contains \(\log(2/\delta)\). One of these is wrong; no factor 2 is introduced in the displayed proof.

11. **Corollary 1 loses \(\delta\) without explanation (lines 499–512, 1241–1250).** Theorem 2 gives coverage \(1-\alpha-\delta\), but the corollary claims \(1-\alpha\). The proof does not reallocate the error budget or take a valid limit. Specify a sequence \(\delta_n\), invoke Theorem 2 with appropriately split levels, and retain the resulting finite-sample coverage statement.

12. **Corollary 1 is missing growth conditions involving \(M\) (lines 500–511).** The condition \(b_n=o(\sqrt n)\) alone does not make
    \(b_n\sqrt{M/n\,\log(1/\delta)}\) negligible if \(M=M_n\) grows. Either declare \(M\) fixed or state the needed joint condition.

13. **Incorrect parentheses in the proof of Corollary 1 (line 1248).**
    `\log(m_b(\hat p))/\alpha` must be `\log(m_b(\hat p)/\alpha)`.

14. **Finite lower-bound parameter notation (line 644).** Write \(p_{n,j}^\dagger=q_n^\dagger\), \(j=1,\dots,M_n\), not \(p_j^\dagger\), because the vector is an \(n\)-indexed sequence.

15. **Set inclusion is reversed in the proof of Theorem 3 (line 1441).** The correct implication is
    \[
    \{T_n\ge q_n^\dagger\}\cap\{q_n^\dagger\ge L_n\}\subseteq\{T_n\ge L_n\}.
    \]
    The displayed probability argument below uses this corrected direction.

16. **The threshold is multiplied twice by \(1-\eta\) (lines 1439–1458).** Since \(L_n=(1-\eta)\log((K_n\vee1)/\alpha)/n\), line 1444 should be `T_n(N)\ge L_n`, not `(1-\eta)L_n`. The later passage to \(\eta_n\to0\) should redefine \(L_n\) using \(\eta_n\); with fixed \(\eta\), \(L_n/q_n^\dagger\to1-\eta\), not 1.

### Unbounded-alphabet results

17. **Theorem 4 is false as stated (critical; lines 698–706).** The deterministic bound \(U_n(\alpha)=1\) is always valid, so it is not true that every candidate fails. The result must concern a *nontrivial* bound, e.g. \(0\le U_n(\alpha)<1\), and failure of coverage is the event \(M_{\max}>U_n\), not \(M_{\max}\ge U_n\). Equality does not violate \(M_{\max}\le U_n\).

18. **The proof of Theorem 4 contains a reversed equivalence (lines 1618–1620).** For \(x=1-(1-c\alpha)^{1/n}\), the condition \(c<1/\alpha\) ensures \(x<1\), not \(x>1\). In addition, choose \(x>U_n(\alpha)\), not merely \(x\ge U_n(\alpha)\), if the failure event is correctly written with `>`.

19. **Malformed theorem statement (line 699).** `U_n(\alpha\in(0,1)` should be `U_n(\alpha)` (with the closing parenthesis), followed by a statement of its range.

20. **Lemma 1 needs special treatment at \(r=1\) (lines 714–720).** The formula contains \(0^0\). Either state \(r>1\), handle \(r=1\) separately, or explicitly adopt the continuous convention.

21. **Lemma 2’s asymptotic expansion is not valid under its stated assumptions (critical; lines 726–740 and 1678–1715).** In line 1693, replacing \((r-1)/r\) by 1 discards
    \[-\frac1r\log\delta=\frac1r\log\frac{n}{r-1},\]
    which is not \(o(1)\) for the proposed choice \(r\asymp\log n\); it converges to approximately 1. A suitable first-order expression is
    \[
    U_n(r,S;\alpha)
    =\frac{r}{en}\left(\frac{Sn}{\alpha r}\right)^{1/r}\{1+o(1)\},
    \]
    under conditions that control the remainder. This missing factor explains why line 738 has \((\log n-\log\log n)/n\), rather than that quantity divided by \(e\). The lemma, its proof, and the optimization discussion must be revised together.

22. **Confidence-level mismatch in Theorem 5’s proof (lines 744–756, 1725–1759).** Equation (for `r_star`) defines \(r^\star\) with \(\alpha\), but the proof needs
    \[
    r_\beta^\star=\log\{S/(\alpha-\beta)\}+\log n-\log\log n
    \]
    to compare it with \(R^\star\). Define this oracle value explicitly in the theorem and impose the monotonicity condition on that same value.

23. **Incorrect monotonicity sentence (line 1757).** For a base greater than 1, \((S/\gamma)^{1/r}\) decreases as \(r\) increases. The intended proof can instead compare the coupled expressions directly: because the prescribed relation gives \(U=eF(r)\), show \(F(R^\star)\ge F(r_\beta^\star)\). Do not claim that both individual factors are nondecreasing.

24. **Notation for \(U\) is reversed in the proof (lines 1727–1777).** The lemma defines `U_n(r,S;alpha)`, but the proof repeatedly writes `U_r(n,S;alpha)`. Standardize to `U_n`.

25. **LaTeX/math typo in the infinite lower bound (line 1868).** `= \le` should be just `\le`.

### Stopping rules

26. **Undefined and inconsistent quantities in the theorem statement (critical; lines 895–948).** In the bounded rule, `delta` is not introduced; it should apparently be `beta`. The sum must run from \(j=1\) to \(M\), not \(j\ge1\) (which diverges because nonexistent coordinates have \(\hat p_j=0\)). Define \(b=b_n\), \(\hat p_{n,j}\), and \(\hat S_n\) explicitly. Also state the convention `inf emptyset = infinity`.

27. **The final assertion uses an undefined generic \(\tau\) (line 947).** State the inequality separately for \(\tau_{\mathrm{bd}}\) and \(\tau_{\mathrm{ub}}\).

28. **The bounded stopping-rule proof drops the unseen indicators (critical; lines 2115–2125).** The event should be bounded through
    \(M_r(N)=\sum_jp_j^r\mathbf1\{N_j=0\}\). As written, \(\sum_jp_j^r\) is deterministic, and the Markov step is invalid. In the last factor, the numerator should be \(r\), not \(n\):
    \[
    \left(\frac{r}{n+r-b}\right)^r
    \left(\frac{n-b}{n+r-b}\right)^{n-b}.
    \]

29. **The oracle quantity varies with \(n\) (lines 2133–2150).** With \(b_n=\log n\), write \(m_{b_n}(p)\), not a fixed \(m(p)\), throughout the oracle rule and event \(\Sigma\). Also, \(b_1=0\) violates \(b\ge1\); start after an admissible minimum sample size or use a modified sequence.

30. **Undefined \(n^\star\) (lines 2207–2214).** These occurrences should be \(\tau^\star\), and the exponent/bound must be evaluated at the same stopping index.

31. **Multiple copy/paste errors in the unbounded stopping proof (lines 2264–2289).** The union-bound sum should be over `n=1,...,n_max`, not `m=1,...,m_max`; the per-time confidence level is \(\beta/n_{\max}\), not \(n_{\max}/\beta\); choose either \(\Sigma\) or \(\Sigma_S\); and every `tau_bd` in this part should be `tau_ub`.

32. **Theorem number hard-coded (line 2217).** Replace “Theorem 5” with `Theorem~\ref{thm:rnorm_unbounded}` so numbering remains correct.

### Numerical sections

33. **Incorrect inversion using Lambert \(W\) (critical; lines 1053–1056).** From
    \[
    S<c\frac Mn\log(M/\alpha),\qquad c=\frac{-\log(1-\alpha)}{\alpha},
    \]
    set \(x=M/\alpha\). The boundary solves
    \(x\log x=Sn/[-\log(1-\alpha)]\), hence
    \[
    M_{\rm crit}
    =\alpha\exp\!\left[W\!\left(\frac{Sn}{-\log(1-\alpha)}\right)\right]
    =\frac{\alpha Sn/[-\log(1-\alpha)]}{W(Sn/[-\log(1-\alpha)])}.
    \]
    The unbounded rule is favored for sufficiently **large** \(M\), not for \(M<5000\). Lines 1055 and 1056 currently contradict one another.

34. **Numerical threshold mismatch (line 976).** With \(n=1000\), \(M=1500\), and \(\alpha=0.05\), the displayed formula gives approximately **15.8634**, not 15.46. Check whether the reported value or one of the inputs differs in the code.

35. **Inconsistent geometric model (lines 1012 and 2296).** The main text uses \(p_j=(1-a)^{j-1}\) (currently with index `i`), while the supplement uses \(p_j=a^j\). These are radically different, especially for \(a=0.005\). Confirm the implementation and use one definition everywhere.

36. **Potentially reversed tail terminology (lines 995 and 2326).** For a standard Zipf/power-law sequence \((j+1)^{-\gamma}\), smaller \(\gamma\) means slower decay/heavier tail. The draft calls smaller \(\gamma\) “light-tailed” and larger \(\gamma\) “heavy-tailed.” If you mean “more extremely rare categories” rather than a heavier mathematical tail, say that explicitly.

37. **BBS interpretation overstates the guarantee (line 1089).** Controlling \(M_{\max}\le\varepsilon\) guarantees that every unseen species has prevalence at most \(\varepsilon\); it does **not** guarantee that the aggregate probability of encountering *any* unseen species is below \(\varepsilon\). Rewrite this sentence in terms of the maximum unseen prevalence.

38. **Table header has two errors (line 2511).** The contamination grid in line 1014 is \(0,10^{-4},5\times10^{-4},10^{-3},2.5\times10^{-3},5\times10^{-3}\). Thus the third header should be `q=0.0005`, not `q=0.005`, and the last should be `q=0.005`, not `q0.005`.

## 2. Language, grammar, and typographical corrections

### Front matter and introduction

- Line 263: `Depertment` → `Department`; `Isreal` → `Israel`.
- Lines 253–281: the comment describing `\blind` is reversed. `\blind=1` currently prints authors; `\blind=0` produces the anonymous title page.
- Title: “with application to” → preferably “with **an** application to” (one application) or “with applications to” (several).
- Line 293: there are 8 keyword phrases. The current JASA instruction is 3–5; reduce the list. The abstract is 161 words and is within the 200-word limit.
- Line 335: “the smallest sample size that guarantees that the probability of the maximum unseen probability ... is smaller” is not meaningful. Prefer the older commented formulation: “the first sample size \(n\) at which the upper confidence bound satisfies \(U_n(\alpha)\le\varepsilon\).”
- Line 372: `approximaetly` → `approximately`.
- Line 382: “It might be then tempting” → “It might then be tempting.”
- Line 385: “data independent” → “data-independent.”
- Line 424: `p_i` → `p_j`, to match the rest of the section.

### Unbounded regime and stopping section

- Line 708: `Propostion` → `Proposition`; “a shift from worst-case, data-independent approach” → “a shift from **a** worst-case, data-independent approach.”
- Line 723: “turn ... into” is acceptable but “convert ... into” is more formal.
- Line 725: lowercase “lemma” in running prose.
- Line 893: `lifter` → `lifted`; “namely there is” → “namely, there is.”
- Line 947: “both ... are stopping times satisfying” → “both ... are stopping times and each satisfies ...”.

### Numerical illustrations and applications

- Line 955: replace “where we compare” with “in which we compare”; replace “either the bounds” with “the bounds.”
- Line 956: “shorter or equal than” → “shorter than or equal to”; “can either yield ... but” → “can yield ...; however,”; “as we experienced that” → “because it performs.”
- Line 964: use either American `favor` consistently or British `favour`; ASA requests American spelling.
- Line 1020: `propotion` → `proportion`.
- Line 1024: hyphenate `M_max-based`; remove the comma after “rules”; “differ for the type of employed estimator, either using using...” → “differ in the estimator used: one uses the estimator of ... and the other that of ...”.
- Line 1025: Suggested rewrite: “Unlike our proposal, the competing stopping rules are heuristic: they do not account for the randomness of the stopping time and do not guarantee that the estimated coverage exceeds 0.99 with arbitrarily high probability.” Check whether “arbitrarily high probability” is actually the intended claim.
- Line 1031: `Tables ... reports` → `Table ... reports`.
- Line 1033: use plural “coverage-based rules” if both competitors are being discussed.
- Line 1045: “analyze each cancer type separated” → “analyze each cancer type separately”; “the different symbols ... corresponds” → “the symbols ... correspond.”
- Line 1046: remove the space before the period; “Even if” → “Although”; remove the double space after the comma.
- Line 1050: consider deleting the question mark and making this an indirect question: “we ask how likely it is ...”.
- Line 1062: “datasets, that are” → “datasets—Breast ..., Lung ..., Skin ..., and Esophageal ...—”; `smaller that` → `smaller than`.
- Line 1063: `neglibile` → `negligible`; the opening “While looking at” clause is awkward. Prefer “Although the accumulation curves ... suggest that new variants will continue to appear ..., the bounds ... indicate that each unseen variant has relatively low prevalence ...”.
- Line 1075: remove the double space before `TCGA`.
- Line 1083: “Application to BBS dataset” → “Application to **the** BBS dataset.”
- Line 1089: “Theorem ..., Section ...” → “Theorem ... **in** Section ...”.
- Line 1090: “which in turn substantially outperforms” → “which, in turn, substantially outperform”; hyphenate `M_max-based`.
- Line 1097: “discovered birds species” → “discovered bird species”; “in BBS dataset, data refer to 2019 only” → “in the BBS dataset; data are from 2019 only.”
- Line 1119: `It includes` → `They include` (or list the supplement explicitly).
- Line 1127: `Fundings` → `Funding`.

### Supplementary proofs

- Line 1169: “The \(r\)-norm ... follows” → “The corresponding \(r\)-norm is”. Also the displayed object should be checked: an \(r\)-norm of the vector \(p_i\mathbf1\{N_i=0\}\) is \((\sum_i p_i^r\mathbf1\{N_i=0\})^{1/r}\), not the norm of the already powered entries as currently written.
- Line 1193: capitalize “Where” only if starting a new sentence; better attach it to the preceding display: “where”.
- Line 1195: `N_u(X^n)` is an unexplained leftover and should be `N_i`.
- Line 1206: use \(M\), not `m`, for the alphabet size.
- Line 1211: `last in equality` → `last inequality`.
- Line 1215: finite-alphabet sum should be over `i=1,...,M`, not `i in N`.
- Line 1243: the domain is \(u\ge0\), not merely \(u\ge-1\), for the way it is used here (the standard inequality is valid for \(u>-1\), with care at \(-1\)).
- Line 1444 onward: alternate use of `T(N)` and `T_n(N)`; standardize to `T_n(N)`.
- Line 1491: add punctuation after the display.
- Line 1553: delete duplicated `any any`.
- Lines 1555–1610: the index `m` and the first `m-1` coordinates equal to one are unnecessary for the construction; a one-active-coordinate example is clearer.
- Line 1580: “\(U(\alpha,n)\) implies that” → “There are two cases.”
- Lines 1585–1591: `T` is undefined; use `U_n(\alpha)` consistently.
- Line 1719: `\hat S` is unbiased, but cite the concentration inequality—not unbiasedness—as the reason for the displayed bound.
- Line 2111: `Firstly` → `First`; “probability of events” → “probabilities of the events.”
- Line 2138: `can not` → `cannot`; “requires the oracle knowledge” → “requires knowledge.”
- Line 2158: “arbitrary close” → “arbitrarily close.”
- Line 2194: “yields to” → “yields.”
- Lines 2228–2229: delete duplicated “yielding / This yields.”
- Line 2262: “arbitrary large” → “arbitrarily large.”

### Supplementary simulations

- Line 2326: add commas around “i.e.”; remove the extra space before “the”; use closed-up em dashes in LaTeX (`---`) consistently.
- Line 2343: lowercase “geometric” unless it is a proper distribution name.
- Line 2371: `upper-bouding` → `upper-bounding`; use sentence case: “Bounded versus unbounded intervals when upper-bounding \(M\).”
- Line 2391: “\(x\)-axis is in \(\log_{10}\) scale” → “The \(x\)-axis is on a \(\log_{10}\) scale.”

## 3. LaTeX and submission cleanup

1. `amsmath` and `graphicx` are each loaded twice. Remove duplicates.
2. The source still contains visible author macros (`\ACtext`, `\ACnote`) and TODOs in the main text. Remove colored revision markup before submission.
3. Labels such as `r-norm lim`, `CI r-norm`, and `E(p) explicit` contain spaces/punctuation. They compile in many setups but are fragile; use labels such as `eq:rnorm-limit`.
4. Prefer `\mathbb P` and `\mathbb E` consistently instead of alternating with bare `P` and `E`.
5. Use `\operatorname{Bin}` or `\mathrm{Bin}` consistently; the draft alternates with `\operatorname{Bin}` in proofs.
6. Replace hard-coded Appendix references with `Supplementary Material` if the supplement will be uploaded as a separate file.
7. The source currently places the supplement after the bibliography in the same document. JASA asks that appendices, proofs, and supporting material be uploaded as a separate anonymized supplement; split the source before submission.
8. The manuscript uses `\spacingset{1.8}` after initially setting 1. JASA asks for fully double-spaced text. Confirm the output has 26 lines of text per page; `1.8` is not automatically equivalent to the requested layout.
9. The custom margins are asymmetric (`oddsidemargin` and `evensidemargin` are changed by different amounts). For a one-sided submission file, use a standard, explicit geometry rather than incremental adjustments.
10. JASA requests 3–5 keywords, but the draft has 8.

## 4. Suggested order of revision

1. Repair Lemma 2 and the resulting unbounded-bound derivation.
2. Restate and reprove the deterministic impossibility theorem with the correct nontrivial range and strict failure event.
3. Repair Theorem 2/Corollary 1 so the advertised endpoint is an exact, computable finite-sample bound.
4. Rewrite the stopping-rule theorem and proof with indexed quantities \(b_n,m_{b_n},\hat p_n,\bar S_n\), then recheck the simulation code against those formulas.
5. Correct the Lambert-\(W\) inversion and numeric threshold.
6. Resolve the geometric-model definition and table headers against the code used to generate the figures.
7. Apply the language edits and clean the LaTeX for submission.
