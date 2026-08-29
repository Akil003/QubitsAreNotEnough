# Correctness Review — *Qubits Are Not Enough* (IEEE TQE submission)

**Reviewer:** independent re-derivation and re-execution of the submission package
**Date:** 2026-08-09
**Package:** `IEEE_TQE/` @ manuscript `Qubits_Are_Not_Enough_TQE.tex` (1218 lines), scripts
`reproduce_study.py`, `revision_study.py`, `minor_revision_study.py`,
`validate_outputs.py`, `dimension_generalization/mesh_dimension_study.py`

**Environment used for re-checks:** Python 3.12.3 (`.venv`), NumPy 2.3.5, SciPy 1.17.0,
pandas 2.2.3. The archived results were produced under Python 3.13.5; all deterministic
quantities re-derived below matched to the digits reported. The package has since been
regenerated end to end under a single environment (see `environment_versions.txt`).

**Baseline integrity at the start of Round 1:** `validate_outputs.py` → *"All citation,
numerical, and complexity integrity checks passed."* `shasum -a 256 -c MANIFEST.sha256` →
79/79 OK. That validator has since been shown vacuous and rebuilt (R2-3).

**Current state:** `validate_outputs.py --full` passes · manifest 81/81 · PDF 19 pages,
0 errors, 0 undefined references · abstract 246/250 · citations 33/33 · every manifest CSV
has a generator. All Round 1 and Round 2 findings are resolved; §6 lists what remains open.

**Round 1 status log** (external-review round is tabulated separately below)

| Finding | State |
|---|---|
| **F1** damage element-stiffness convention | **RESOLVED** — code, data, figure, and all seven manuscript sites updated; see the banner in §F1 |
| **F2** Davis–Kahan factor and precondition | **RESOLVED** — factor 2 restored, grid refined, §VII-A rewritten; see the banner in §F2 |
| **F6.6** SNR quoted outside linearization validity | **RETIRED** — dissolved by the F1 fix (SNR is now 1.44, not 0.13) |
| **F3** two CSVs with no generating code | **RESOLVED** — `restart_study()` and `theory_beta_study()` shipped; every CSV in the manifest now has a generator |
| **F4** Algorithm 1 never executed | **RESOLVED** — executed under the Gershgorin policy; α certification guard added |
| **F5** two exponents for one curve | **RESOLVED** — fit ranges disclosed at both sites |
| **F6.1–F6.5, F6.7, F6.8** | **RESOLVED** — see the banner in §F6 |
| **F7** Monte Carlo depends on Pauli term ordering | **RESOLVED** — `kind="stable"` on both Pauli routines; numbers re-based |
| **S1** lumped-mass Pauli closed form | **ADDED** — Proposition 4, with proof |
| **S2** normalized gradient statistic | **ADDED** — §VII-E |
| Test E scope | **ADDED** — Limitations now states the damaged states are exact eigenvectors of a known model |

---

## Round 2 — independent external review (2026-08-09)

Three independent reviewers audited the *corrected* package under separate lenses
(mathematics; code and reproducibility; structural-dynamics and quantum-algorithm domain).
Every CRITICAL and MAJOR finding below was re-verified locally before being accepted; that
process rejected two reviewer claims outright. Nine findings trace to Round 1 work.

### Accepted and fixed

| ID | Finding | Resolution |
|---|---|---|
| R2-1 | **The deflation threshold is prior art.** Higgott et al.\ (arXiv:1805.08138) Eq. (2) is the deflated cost and §3 states "it suffices to choose $\beta_i > E_k - E_i$" — the manuscript's Eq. (22), while §I claimed it as *the* contribution | Attribution added at four sites; novelty repositioned onto the necessity converse, the padded-block extension, the certified rule, and the Davis–Kahan bound |
| R2-2 | **The "smallest certified penalty" rule is backwards.** $\beta$ enters the *denominator* via $g_r$; the bound is minimized at $\beta^\star=(\varepsilon_r-\varepsilon_j)+G_r$ and diverges at the threshold. Verified: actual error falls monotonically in $\beta$, and $\beta=3$ beats the recommended value 10× | §V-C rewritten with the correct optimum (new Eq. betastar); §V-B warns the Ritz policy is *worse* for robustness; §IX-B and Appendix B corrected |
| R2-3 | **`validate_outputs.py` could not detect a code regression.** 7 injected bugs → "All checks passed". Independently reproduced with 2 | New `validate_regeneration()` re-executes generators into a temp dir and diffs at `rtol=0`. Now catches 7/7 across two tiers (`--full`), no false positives, 7.5 s |
| R2-4 | **SHM section has no deployment path** — preparing $\psi^d_i$ needs the damaged model, or reduces to an $O(1)$ classical dot product | Caveat paragraph added at Eq. (qmse); §VII-G retitled to "estimable, but only given the damaged model" |
| R2-5 | **$3.1\times10^8$ over-counts ~4.3×.** All element supports ⊆ the Hamiltonian's 19 strings; $\sum_e\Hq_e=\widetilde\A$ to 1.1e-16; distinct settings 152 | §IX-D now brackets $7.3\times10^7$–$3.1\times10^8$; abstract and conclusion use ${\approx}10^8$ |
| R2-6 | **"2950 Pauli terms" is set by an undocumented absolute cut.** 2316/2950/3640/4251 at $\tau=10^{-8..-14}$; exponent moves 1.36→1.65 | Threshold and its sensitivity now stated in §VI-A |
| R2-7 | **Shot cost scales as $O(N^4)$**, not $N$-independent $\varepsilon$ (relative frequency accuracy needs $\varepsilon=O(N^{-2})$). Verified: mode ratio 31.9 vs predicted 31.6 | Added to §VI-D as $O(kN^6N_{\rm eval})$ |
| R2-8 | **Table 1 omits state preparation** despite Remark 3 proving $L=\Omega(N/\log N)$ | Row plus footnote added; columns declared non-commensurable |
| R2-9 | **Theorem 1 fails for degenerate spectra** — the case §VIII-C showcases | New Remark 1 on degenerate targets and block deflation |
| R2-10 | **Contributions (iv), (v) never implemented** | Demoted to "proposed and analyzed here, but not evaluated experimentally" |
| R2-11 | **Lumped-vs-consistent discretization error never measured.** Computed: 0.06% on modes 1–4 at $N=128$, vs 2.45% shot noise — the expensive mass model buys unusable accuracy | Added to §VIII-C as an argument *for* the lumped operator |
| R2-12 | **Damage signal vs measurement noise never compared.** Modes 1–2 have SNR 0.49 / 0.43 at the reported budget | Added to §VII-G |
| R2-13 | **Eq. (ansatz) does not describe the implemented circuit** (as typeset it has $n_q$ redundant parameters and ends with an entangler) | Rewritten to the implemented $R\cdot(\mathrm{Ent}\cdot R)^L$ form |
| R2-14 | **Variance computed by catastrophic cancellation**, returning negative values then clamped (−5.8e-18 on an exact eigenvector) | Replaced by $\lVert(\Hq-E)\psi\rVert^2$; clamps removed |
| R2-15 | **Two different six-DOF benchmarks**, neither specified in the paper | Both now defined in §IX as Benchmark S and Benchmark C |
| R2-16 | **Stubbs/Cornwell convention imprecise** — those formulations hold the undamaged stiffness in both numerator and denominator and compare fractions as a *ratio*; \cref{eq:mse} keeps the current $\K$ in the denominator and \cref{eq:damageindex} takes an absolute difference | §VI-D now states the difference explicitly, including that an absolute difference is biased toward high-strain-energy elements |
| R2-17 | **Mode tracking presented as novel**; $\MAC^{(M)}$ needs $\M$ at every DOF, which output-only identification cannot supply; the $\chi^2$ threshold is not calibrated for $1-\MAC$; "unbiased" overstated | Attributed to the automated-OMA literature; $\chi^2$ replaced by an empirical quantile; sensor-space MAC restriction and the mode-veering / missing-label failure modes now stated |
| R2-18 | **Damage location implied swept but is fixed** at element 3 | Reworded, with the exact margin range across all six locations (0.515–0.702) showing element 3 at 0.654 is mid-range, not favourably selected |
| R2-19 | **Table 4 juxtaposes Monte Carlo and analytic rows**, so the quoted 2.5–2.8 improvement factors do not match a direct row division | Footnote added distinguishing the two, plus the 400-trial relative standard error (~3.5%, i.e. mode 1 good to ±0.09) |

### Rejected after verification

| Claim | Why rejected |
|---|---|
| Element observables cost $\Theta(N\cdot N_P)$ extra settings (Reviewer C) | **False.** Their supports are subsets of the Hamiltonian's: 0 extra strings at $N=8,16$; 54 of 511 at $N=32$. Reviewer A's opposite reading is correct, and is what R2-5 implements |
| A uniform lumped chain has $2N-1$ Pauli terms (Reviewer A) | **Wrong number.** Measured $3\cdot2^{n_q-1}-1$ (5, 11, 23, 47, 95). The *substance* — genericity fails for uniform meshes — is correct and is now stated in §VI-A |
| "Frequency-based detection is not possible at that budget" (Reviewer C) | **Overstated.** Mode 3 has SNR 13.1. Corrected in the §VII-G text |

### Errors introduced in Round 1

R2-2, R2-5 (partly), R2-9, R2-11 and R2-14 predate this review, but the following were mine:
the Appendix B optimality claim (R2-2); the "$<2\times10^{-13}$%" bound, actual max $9.2\times10^{-13}$;
"indistinguishable from guessing at $10^3$" (binomial $p=0.013$); "local exponent falls
monotonically" (it rises at the first step); an incomplete F1 fix leaving `reproduce_study.py`
on the oracle convention; Proposition 4's unstated genericity gap; a "single pass" provenance
claim contradicted by timestamps; leaving the validator vacuous; and not checking Theorem 1's
attribution. All are fixed. The common failure mode was checking internal consistency
(number matches CSV) without checking external validity (is the claim about the number true).

### Also corrected in Round 2

Smaller manuscript items, each verified before change:

| Item | Fix |
|---|---|
| $\beta$ denotes both the deflation penalty and the type-II error rate | Table 8 caption disambiguates, and records that detection is two-sided while localization is one-sided with Bonferroni |
| Table 2 header reads $s_G/\lambda_{\max}$ but the text defines $U_G$ | Header corrected to $U_G/\lambda_{\max}$ |
| $\varepsilon_{N+1}:=\alpha$ undefined when $D=N$ (every power-of-two chain) | Theorem 1 now sets $\varepsilon_{N+1}:=+\infty$ in that case |
| "2D separator fill scales as $O(\sqrt N)$" confuses separator *size* with *fill* | Corrected to size $O(\sqrt N)$/$O(N^{2/3})$, fill $O(N\log N)$/$O(N^{4/3})$ |
| \cref{eq:shotsdetect} tests a null under which its own linearization fails ($D_e$ is a sum of absolute values, folded and biased upward by $\sigma\sqrt{2/\pi}$ at $D_e=0$) | Stated, with a recommendation to calibrate the null by permutation |
| Appendix C's $\sigma_H^2\le10^{-10}$ is unattainable at any shot budget, and $\langle\Hq_p^2\rangle$ is not itemized in the resource vector | Both now stated |
| The degeneracy demonstration is analytically tautological ($\varepsilon_{\rm sub}=0$ follows in closed form) | Reframed as verifying the implementation; the untested near-degenerate case is named |

Code hygiene, each confirmed numerically neutral before acceptance:

| Item | Fix | Evidence |
|---|---|---|
| `sample_operator` was dead code duplicating `prepare_operator`+`sample_prepared` | Removed | 0 references outside its own definition |
| The gradient sweep built each state twice per shift, doubling the cost of the most expensive script | Bound to locals | max abs difference **0.000e+00** at $n_q=3,5$ |
| `pauli_coefficients` recomputed the parity vector for all $4^n$ strings though it depends only on the $2^n$ distinct z-masks | Cached | validator passes at `rtol=0`; $n_q{=}8$ decomposition 6 s → **0.6 s** |
| The shared module RNG makes results order-of-call dependent | Documented at the definition, with the instruction that new studies take their own generator | — |

Documentation and provenance:

| Item | Fix |
|---|---|
| `mesh_dimension_study.py` absent from the README block, `make study`, and the provenance list, though it produces Table 5 and two figures | Added to all three, with its ~35 min runtime and its dependency on `revision_study.py` running first |
| "Regenerated in a single pass" contradicted by file timestamps | Replaced with an accurate statement, including why `dimension_generalization/data/` is older but current |
| `make clean` deleted `Qubits_Are_Not_Enough_TQE.bbl`, a manifested artefact | Removed from the clean target, with a comment |
| `PDF_PREFLIGHT.txt` and `VERIFICATION_REPORT.md` page counts stale | Synced to 19 |

### Final verified state

`validate_outputs.py --full` passes · `MANIFEST.sha256` 81/81 · PDF 19 pages, 0 errors,
0 undefined references · abstract 246/250 words · citations 33/33, none missing, none
uncited · every manifest CSV has a generator · no stale numbers in the source.

The validator was re-tested by injection after rebuilding: **7 of 7 severe bugs caught**
across its two tiers (sign flip in `pauli_expectation`, corrupted Cholesky whitening,
`qwc_groups` hard-return, off-by-one in `sample_prepared`, Gershgorin dropping its
off-diagonal sum, `tau -> 0` in the certified schedule, and a flipped chain stiffness
gradient), with no false positive on the clean package.

---

## 0. Verdict

The theory is sound. Propositions 1–3 and Theorem 1 are correct as stated, and Theorem 1 is
numerically sharp to machine precision. The central thesis — that logarithmic qubit count is
not sufficient, and that measurement is the binding constraint — is well supported and
survives everything below.

One result did not survive: **the finite-shot damage-localization outcome (§VII-G, Table 8,
abstract, conclusions) was an artifact of an element-stiffness convention that contradicts the
modal-strain-energy literature the paper cites and cannot be implemented without knowing the
damage.** Correcting it moved 10%-damage localization from 26% to 83% and cut the shot
requirement by ~118×. **This has since been fixed — see the banner in §F1.** A second issue,
the Davis–Kahan numbers, was a reporting error rather than a modelling one, **also now fixed
— see the banner in §F2**. Everything else is
minor or editorial.

---

## 1. What was verified as correct

### 1.1 Analytical results

| Claim | Status | Note |
|---|---|---|
| Prop. 1 (mass–Euclidean isometry) | correct | $\mathbf A^T=\mathbf L^{-1}\mathbf K^T\mathbf L^{-T}=\mathbf A$; $\mathbf y_i^T\mathbf y_j=\boldsymbol\phi_i^T\mathbf L\mathbf L^T\boldsymbol\phi_j=\boldsymbol\phi_i^T\mathbf M\boldsymbol\phi_j$ |
| Prop. 2 (spectral-safe padding) | correct | $\lambda_{\max}(\mathbf A/s_t)\le U_G/s_t=U<\alpha$ for any $s_t>0$ |
| Prop. 3 (scalar-scaling law) | correct | eigenvectors, $\kappa$ invariant; $\nabla E$, $\mathrm{Var}$ scale as $1/c$, $1/c^2$ |
| Theorem 1 (deflation threshold) | correct, both directions | see §1.2 |
| Eq. (25) Hoeffding radius | correct | $t=\sqrt{2\log(2/\alpha)\sum h_\ell^2/n_\ell}$ follows from Hoeffding with per-sample coefficient $h_\ell/n_\ell$, range $2h_\ell/n_\ell$ |
| Eq. (54)–(56) shot variance / optimal allocation | correct | $n_\ell\propto\lvert h_\ell\rvert\sqrt{1-\mu_\ell^2}$ is the Lagrange stationary point |
| Eq. (57) delta method | correct | $\partial f/\partial\varepsilon = s_t/(8\pi^2 f)$ |
| Eq. (52) quantum MSE identity | correct | $\mathbf y^T\mathbf L^{-1}\mathbf K_e\mathbf L^{-T}\mathbf y=\boldsymbol\phi^T\mathbf K_e\boldsymbol\phi$ |
| Eq. (63) $\mathrm{MAC}=\cos^2\theta$ | correct | confirmed numerically to machine precision |
| Eq. (65) depth bound rearrangement | correct | but see finding **F6.3** on how it is used |
| Graph-diameter estimates (4160 / 81,000) | correct | $2(\sqrt N-1)=127$, $3(N^{1/3}-1)=127$ |

### 1.2 Theorem 1 is numerically sharp

From `revision_data/penalty_threshold.csv`:

| target mode | max overlap for $\beta<\beta_{\min}$ | overlap at $\beta=\beta_{\min}$ | min overlap for $\beta>\beta_{\min}$ |
|---|---|---|---|
| 2 | $2.71\times10^{-26}$ | 0.368 (degenerate) | 1.000000 |
| 3 | $2.54\times10^{-26}$ | 0.465 (degenerate) | 1.000000 |
| 4 | $9.79\times10^{-28}$ | 0.001 (degenerate) | 1.000000 |

Exactly the necessary-and-sufficient behaviour claimed, with the expected indeterminacy at
equality. This is the strongest result in the paper and it is airtight.

### 1.3 The depth / expressibility result is robust

The manuscript's Table 3 was generated with `starts=2, maxiter=160`
(`revision_study.py:551-554`), which invites the objection that shallow-circuit failure is
optimizer failure. It is not. Re-run with 40 restarts, `maxiter=1500`, `ftol=1e-15`:

| $N$ | $n_q$ | $L$ | params | need $2^{n_q}-1$ | 40-restart best freq. error |
|---|---|---|---|---|---|
| 8 | 3 | 1 | 6 | 7 | $2.7333\times10^{1}$ % |
| 8 | 3 | 2 | 9 | 7 | $1.68\times10^{-14}$ % |
| 16 | 4 | 1 | 8 | 15 | $8.1751\times10^{1}$ % |
| 16 | 4 | 2 | 12 | 15 | $5.1585\times10^{0}$ % |
| 16 | 4 | 3 | 16 | 15 | $4.22\times10^{-12}$ % |

The failing rows reproduce to five significant figures. Only the $L=3$ row improves
(manuscript reports $4.09\times10^{-9}$%, true floor $\approx4\times10^{-12}$%), i.e. that
single entry is optimizer-limited, not expressibility-limited.

### 1.4 Implementation spot-checks

- `pauli_coefficients` computes $\mathrm{tr}(PH)/2^{n_q}$ correctly via the signed-permutation
  identity $P=i^{n_Y}X^{x}Z^{z}$; the index reshuffle $i\mapsto i\oplus x$ makes the code's sum
  identical to the trace. Odd-$Y$ strings vanish for real symmetric $H$, so the `.real` cast is safe.
- `qwc_groups` conflict mask `non_a & non_b & ((x_a^x_b)|(z_a^z_b))` is the correct qubit-wise
  commutation test.
- Parameter-shift gradient `0.5*(E(θ+π/2) − E(θ−π/2))` is correct for $R_y=e^{-i\theta Y/2}$.
- `qwc_optimal_variance_coefficient` correctly computes
  $\min\sum_g v_g/n_g$ s.t. $\sum n_g=B$ $\Rightarrow$ $(\sum_g\sqrt{v_g})^2/B$, with exact
  within-group covariances via `_qwc_product`.
- Bit ordering is MSB-first and consistent across `cnot_matrix`, `kron_all`, and `label_masks`.
- Ratio-variance and margin-variance terms neglected in `minor_revision_study.py` are genuinely
  the beneficial ones ($\sum_e a_e = b$ exactly, so $\mathrm{cov}(a_e,b)>0$), so "conservative"
  is justified.

### 1.5 Reported values that match their data

Tables 4 (finite-shot RMSE + QWC floor), 5 (2D/3D scaling), 6 (Gershgorin range 1.04–1.14),
and 8 (delta-method shot requirements) match their CSVs exactly. Wilson intervals for 52/200
and 72/200 are 20.4–32.5% and 29.7–42.9% as stated. Degeneracy demo: $\mathrm{MAC}=\cos^2\theta$
to machine precision, $\varepsilon_{\rm sub}\in[2.58,3.65]\times10^{-8}$, consistent with the
quoted $\sim3\times10^{-8}$ numerical floor.

---

## 2. Findings

Ordered by impact on the manuscript's claims.

---

### F1 — Damage-localization result is an artifact of the element-stiffness convention

> ### ✅ RESOLVED (2026-08-09)
>
> Applied as described under *Recommended action* below, with the sensitivity case retained.
>
> - `revision_study.py:643` and `minor_revision_study.py:240` now build the damaged-state
>   numerator operators from `Ke0`, each with an explanatory comment; the
>   `damage_shot_requirement_study` docstring gained a matching assumption bullet.
> - Eq. (49) and (50) now carry the $\K_e^0$ superscript, with three sentences explaining why
>   the baseline operator is reused and why the denominator is unaffected. This also removes
>   the Eq. (48)/(50) inconsistency.
> - §VII-G retitled *"Finite-shot damage localization is feasible but shot-expensive"* and
>   rewritten; the oracle-convention numbers are preserved as **Remark 3**.
> - Abstract, Test E, Table 8, §IX-D, and the conclusion updated. Abstract is 249/250 words.
> - Data regenerated on Python 3.12.3 (see `environment_versions.txt`, Note 2). Pre-fix
>   outputs archived as `revision_data/damage_shot_localization_oracle_convention.csv` and
>   `minor_revision_data/damage_shot_requirements_oracle_convention.csv`, both added to the
>   manifest so Remark 3 is backed by committed data.
> - `validate_outputs.py` constants updated; PDF rebuilt clean; `MANIFEST.sha256` regenerated
>   (81/81 OK).
>
> **Final numbers** (after the single-provenance regeneration). 10% loss at $10^5$ shots/term:
> 25.0% → **84.5%** (Wilson 78.8–88.9%), median margin −0.413 → **+0.283**. 20% loss:
> 29% → **100%**. Localization requirement at 10%: $5.71\times10^7$ → **$4.83\times10^5$**
> shots/term, i.e. $3.7\times10^{10}$ → **$3.1\times10^8$** circuit executions. SNR at $10^5$:
> 0.13 → **1.44**. Exact margin ratio **11.18×**, shot-requirement ratio **118×** — both
> deterministic and unaffected by the F7 re-basing.
>
> The oracle-convention comparison is now produced by shipped code
> (`damage_shot_study("oracle")`, `damage_shot_requirement_study("oracle")`) rather than by
> hand-copied files, so Remark 3 is reproducible.

**Severity:** major — changes abstract, §VII-G, Fig. 11, Table 8, §IX-D, conclusions.
**Location:** `revision_study.py:641-644`, `minor_revision_study.py:237-241`; Eq. (50).

#### What the code does

```python
damaged = exact_modal_data(damage_story=2, damage_fraction=severity)
_, _, Ked, _, Ld, lamd, Yd, sd, Hd = damaged
elemd = transformed_element_operators(Ked, Ld, sd, Hd.shape[0])   # <-- Ke of the DAMAGED model
```

`Ked` comes from `chain_model(..., damage_story, damage_fraction)`, whose element matrices
already carry the $(1-d_e)$ factor. So the implemented statistic is

$$\eta_e^{d}=\frac{\boldsymbol\phi_d^T(1-d_e)\mathbf K_e^{0}\boldsymbol\phi_d}{\boldsymbol\phi_d^T\mathbf K^{d}\boldsymbol\phi_d}.$$

#### Why this is wrong

1. **It is not the cited method.** Stubbs (1995) and Cornwell (1999), cited at
   `Qubits_Are_Not_Enough_TQE.tex:614`, both form the strain-energy fraction with the
   *undamaged* $\mathbf K_e^{0}$ for the baseline **and** the damaged state, precisely because
   the damaged element stiffness is the unknown.
2. **It is not implementable.** Measuring $\langle\Hq_e\rangle$ with
   $\Hq_e=\mathbf L^{-1}\mathbf K_e^{d}\mathbf L^{-T}/s_t$ requires already knowing $d_e$.
3. **It is internally inconsistent.** Eq. (48), two lines above Eq. (50), correctly writes
   $\partial\lambda_i/\partial d_e=-\boldsymbol\phi_i^T\mathbf K_e^{0}\boldsymbol\phi_i$ with the
   $\mathbf K_e^{0}$ superscript. Eq. (50) drops it, and the code follows Eq. (50) literally.
4. **It cancels most of the signal.** The $(1-d_e)$ numerator shrinkage opposes the denominator
   drop and the mode-shape change.

Note the denominator is *not* affected: $\boldsymbol\phi^T\mathbf K^{d}\boldsymbol\phi=\lambda_d$
for mass-normalised modes, i.e. it is the VQE energy itself and is legitimately available.
Only the numerator operator needs to change.

#### Quantified effect

Exact (noise-free) four-mode index $D_e$, 10% loss at element 3:

```
as coded : [0.00824 0.01070 0.01251 0.01045 0.00400 0.00693]  margin +0.00181
K_e^0    : [0.00824 0.01070 0.03094 0.01045 0.00400 0.00693]  margin +0.02025
```

Monte Carlo localization accuracy (200 trials, identical seeds, only the numerator operator
changed; random baseline 16.7%):

| severity | shots/term | as coded | with $\mathbf K_e^{0}$ |
|---|---|---|---|
| 10% | $10^4$ | 24.5% | **45.5%** |
| 10% | $10^5$ | 25.5% (paper: 26%) | **85.5%** |
| 20% | $10^4$ | 24.5% | **65.5%** |
| 20% | $10^5$ | 34.0% (paper: 36%) | **99.0%** |

Delta-method requirements (Table 8 recomputed):

| severity | detect, coded | detect, $\mathbf K_e^{0}$ | localize, coded | localize, $\mathbf K_e^{0}$ |
|---|---|---|---|---|
| 2% | $1.32\times10^{7}$ | $2.55\times10^{6}$ | $1.27\times10^{9}$ | $1.52\times10^{7}$ |
| 5% | $2.04\times10^{6}$ | $3.85\times10^{5}$ | $2.11\times10^{8}$ | $2.23\times10^{6}$ |
| 10% | $4.82\times10^{5}$ | $8.75\times10^{4}$ | $5.71\times10^{7}$ | $4.83\times10^{5}$ |
| 20% | $1.10\times10^{5}$ | $1.82\times10^{4}$ | $1.83\times10^{7}$ | $9.01\times10^{4}$ |

Localization cost falls by 84–203×. The "$\approx3.7\times10^{10}$ circuit shots" figure in
§IX-D becomes $\approx3.1\times10^{8}$.

#### Recommended action

1. Fix Eq. (50) to $\eta_{e,i}=\boldsymbol\phi_i^T\mathbf K_e^{0}\boldsymbol\phi_i/\boldsymbol\phi_i^T\mathbf K\boldsymbol\phi_i$ and Eq. (49) to $\Hq_e=\mathbf L^{-1}\mathbf K_e^{0}\mathbf L^{-T}/s_t$, stating explicitly that the *baseline* element operator is reused on the damaged structure because the damaged model is unknown.
2. Change `Ked` → `Ke0` in `revision_study.py:643` and `minor_revision_study.py:237`.
3. Regenerate `damage_shot_localization.csv`, `damage_shot_requirements.csv`, Fig. 11, Table 8.
4. Update the abstract sentence, §VII-G, and the conclusion.
5. Optionally keep the current numbers as a labelled sensitivity case ("if the damaged element operator were known exactly, cancellation reduces the margin by ~11×") — that is a legitimate and interesting secondary observation.

#### Does this break the paper?

No. $3\times10^{8}$ shots for one localization decision on a six-DOF chain is still
prohibitive, the measurement bottleneck argument in §VI-D and §IX-C is untouched, and the
$1/f_i^2$ first-mode sensitivity result is independent of this. What must go is the specific
"barely above random guessing" framing. The honest replacement is roughly: *reliable
localization of a 10% loss needs $\sim5\times10^{5}$ shots per Pauli term, i.e. $\sim3\times10^{8}$
circuit executions per decision epoch.*

---

### F2 — Reported Davis–Kahan bounds are not Eq. (33), and two of three violate its precondition

> ### ✅ RESOLVED (2026-08-09)
>
> - `revision_study.py:329` now computes `2.0 * perturb / separation`, matching the boxed
>   Eq. (33); a `weyl_condition_met` column records whether $\norm{\Delta_r}<g_r/2$ holds.
> - The angle grid was refined below $g_r/2$ (`linspace(0, 0.035, 8)` plus the original coarse
>   tail), turning **one** valid data point into **seven**. The bound holds at every one and is
>   uniformly conservative by a factor of about two — a cleaner result than the original text
>   claimed.
> - §VII-A rewritten accordingly, with a second paragraph disclosing that the perturbation is
>   adversarial by construction (the first projector is rotated toward the target mode itself,
>   `q = Yp[:, (j+r) % 6]` with `j=0, r=3`), which is what makes a 0.02 rad projector error
>   produce a 0.249 target rotation.
> - `deflation_perturbation.csv` regenerated; PDF rebuilt clean; manifest reverified 81/81.
>
> | angle | $\norm{\Delta_r}$ | $\sin\theta$ | Eq. (33) | Weyl met |
> |---|---|---|---|---|
> | 0.005 | 0.00341 | 0.0677 | 0.1363 | ✓ |
> | 0.010 | 0.00682 | 0.1330 | 0.2727 | ✓ |
> | 0.015 | 0.01023 | 0.1940 | 0.4090 | ✓ |
> | 0.020 | 0.01363 | 0.2492 | 0.5453 | ✓ |
> | 0.025 | 0.01704 | 0.2985 | 0.6816 | ✓ |
> | 0.030 | 0.02045 | 0.3418 | 0.8179 | ✓ |
> | 0.035 | 0.02386 | 0.3796 | 0.9542 | ✓ |
> | 0.040–0.20 | 0.0273–0.135 | 0.413–0.712 | vacuous | ✗ |
>
> **Note.** Rerunning `theorem_penalty_study()` also rewrites `penalty_threshold.csv`. That
> file was compared and restored: the only non-float-noise differences were the three rows at
> exactly $\beta/\beta_{\min}=1$, where the ground state is degenerate and the overlap is
> therefore arbitrary. Theorem 1's sharpness result is unchanged.

**Severity:** significant — §VII-A, one paragraph.
**Location:** `revision_study.py:329`; Eq. (32)–(33).

The manuscript boxes

$$\sin\angle(\widehat u_r,u_r)\le\frac{2\lVert\Delta_r\rVert_2}{g_r}$$

conditioned on $\lVert\Delta_r\rVert_2<g_r/2$. The code computes

```python
bound = min(1.0, perturb / max(separation, 1e-15))
```

i.e. $\lVert\Delta_r\rVert/g_r$ — **the factor 2 is missing**. The quoted "0.273, 0.545, 0.818"
are these halved values.

| angle (rad) | $\lVert\Delta_r\rVert$ | $\sin\theta$ observed | CSV "bound" $=\lVert\Delta\rVert/g$ | Eq. (33) $=2\lVert\Delta\rVert/g$ | $\lVert\Delta\rVert<g/2=0.025$? |
|---|---|---|---|---|---|
| 0.02 | 0.013633 | 0.24924 | 0.27267 | 0.54534 | yes |
| 0.04 | 0.027261 | 0.41250 | 0.54523 | 1.09045 | **no** |
| 0.06 | 0.040878 | 0.50810 | 0.81757 | 1.63514 | **no** |

Two of the three reported rows lie outside the regime in which the paper's own theorem
applies. A reader recomputing Eq. (33) will get 0.545 / 1.09 / 1.64 and conclude the numbers
were mis-transcribed.

Note the underlying physics is right: for $\theta=0.02$, the contaminating coupling
$\beta\sin\theta\cos\theta\approx0.0136$ against the gap $g_r=0.05$ predicts a mixing angle
$\tfrac12\arctan(2\times0.0136/0.05)=0.246$ rad, matching the observed 0.249. The amplification
narrative is correct; only the bound arithmetic is not.

Also worth stating explicitly in the text: the perturbation deliberately mixes projector $j=0$
with the *target* mode ($q = Y_p[:,(j+r)\%6]$ with $j=0,r=3$ gives $q=u_3=u_r$). That is a
worst-case stress test and explains the 12× amplification; readers will otherwise find the
$0.02\text{ rad}\to0.249$ jump implausible.

**Recommended action:** report $2\lVert\Delta_r\rVert/g_r$ and restrict the reported rows to
$\lVert\Delta_r\rVert<g_r/2$ (angle 0.02, and add finer angles below it), **or** state that a
different Davis–Kahan variant is used and derive it — noting that the $\lVert\Delta\rVert/\delta$
form requires $\delta$ measured against the *perturbed* spectrum, which `revision_study.py:309-312`
does not do (it uses the exact deflated spectrum). Add one sentence describing the
target-direction contamination.

---

### F3 — Two reported results have no shipped generating code

>### ✅ RESOLVED (2026-08-09)
>
> `restart_study()` and `theory_beta_study()` added to `revision_study.py` and wired into
> `main()`; **every CSV in the manifest now has a generator.** The original seeding of
> `optimizer_restarts.csv` proved unrecoverable (two plausible schemes scanned, no match) and
> the archive would not bit-reproduce here regardless, so the study was re-specified: 20
> independent initializations per depth, initialization *k* drawn from `default_rng(k)`.
> It reproduces §VII-D's claim exactly — min **127.986**, median **130.200**, max **130.200**
> (9 of 20 land on the lower attractor, so the median sits in the upper cluster). Only the
> depth-2 median moved, $4.3\times10^{-12}$ → $4.1\times10^{-12}$%.
>
> Two further CSVs were in this category and are also fixed: the `*_oracle_convention.csv`
> files, originally hand-copied during F1, are now emitted by the `convention` parameter of
> the two damage studies. `validate_outputs.py` gained assertions for both new files.


**Severity:** reproducibility / editorial policy.

`revision_data/optimizer_restarts.csv` (41 lines) and `revision_data/six_dof_theory_beta.csv`
are present in the archive but **no shipped script writes them**. The §VII-D claim —

> "Depth 1 has a minimum first-mode error of 127.99%, a median of 130.20%, and a maximum of
> 130.20%; depth 2 has a median error $4.3\times10^{-12}$%."

— comes entirely from `optimizer_restarts.csv`. Its contents do support the claim
(depth 1: n=20, min 127.986, median 130.200, max 130.200; depth 2: min $3.31\times10^{-13}$,
median $4.320\times10^{-12}$, max $6.35\times10^{-11}$), but a referee cannot regenerate it.

This contradicts the Data-and-Code-Availability statement, which claims the archive contains
the scripts required to reproduce the exact-state studies.

**Recommended action:** add a `restart_study()` function to `revision_study.py` that writes both
files, call it from `main()`, and add a `validate_outputs.py` assertion on the depth-1 median.

---

### F4 — Certified padding and certified deflation are never actually exercised

>### ✅ RESOLVED (2026-08-09)
>
> **F4.1.** `exact_modal_data` now sets `alpha = max(1.10, gersh(A)/s + 1e-3)`. The six-DOF
> chain has $U=1.0920$, so α stays at 1.10 and **no reported number changes**; the N=4 chain
> ($U=1.1102$) would previously have been silently uncertified and now is not.
>
> **F4.2.** Algorithm 1 is now executed, by `theory_beta_study()`, under the quantum-only
> Gershgorin policy — no auxiliary classical eigensolve. `solve_vqd` gained a
> backward-compatible `beta_schedule(r, j, energies, variances)` hook (verified bit-identical
> on four configurations including the `gamma≠0` path), so $\delta_j$ comes from the run's own
> residual. Certified penalties land in **0.68–0.98**, three to four times *smaller* than the
> ad-hoc $\beta=3$, which is exactly what \eqref{eq:deflperturb} predicts. At depth 2 all four
> modes are recovered to $<2\times10^{-13}$%, MAC = 1 to twelve digits, overlap
> $<4\times10^{-15}$, leakage $<10^{-15}$ — passing the Appendix C acceptance set with orders
> of margin. Reported in Appendix B.


**Severity:** claim-vs-code gap. Does not affect any reported number, but is a latent bug and
undercuts the paper's central methodological argument.

**F4.1 — Hardcoded $\alpha$.** `revision_study.py:583`:

```python
H = pad(A / s, 2 ** math.ceil(math.log2(n)), 1.10)      # s = lam[-1], the EXACT lambda_max
```

Two problems. First, $s_t$ is the exact $\lambda_{\max}$ from a full `eigh`, not the
Lanczos/power estimate §IV-A describes. Second, $U=U_G/s_t$ is never computed, and $\alpha$ is
hardcoded rather than set to $U+\delta_p$. For the 6-DOF chain $U_G/\lambda_{\max}=1.091988$,
so $\alpha=1.10$ happens to be certified with $\delta_p=0.008$ — a 0.8% margin, by luck. The
$N=4$ chain has ratio $1.110213$; the same hardcoded $\alpha$ would be **uncertified** there.

`reproduce_study.py:277` does it correctly (`scale = gershgorin_upper(A)`, `alpha=1.25`, so
$U=1$, $\delta_p=0.25$). The two scripts disagree on the paper's own procedure.

**F4.2 — Hardcoded common $\beta$.** `reproduce_study.py:159` defaults to `beta: float = 3.0`
and `run_case` never overrides it. This is precisely the "arbitrary common penalty" that §V-B
and §IX-B argue against. It is safe here (eigenvalues $\le1$ after normalisation), but
Algorithm 1's certified rule $\beta_j=U_r-\widehat\varepsilon_j+\delta_j+\tau$ is **never
executed anywhere in the package**.

**Recommended action:** replace the literals with `alpha = gersh(A)/s + delta_p` and
`beta_j = U_r - eps_hat_j + delta_j + tau`; add an assertion `alpha > gersh(A)/s`. Even one
worked execution of Algorithm 1 under the Gershgorin policy would close the gap between the
framework and the demonstration.

---

### F5 — Two different scaling exponents reported for the same curve

>### ✅ RESOLVED (2026-08-09) — §VI-D now says the 1.4–1.5 figure is the $N\ge16$ fit and adds
> the falling local exponent (1.25 at the last doubling); §VIII-B states that 1.66/1.69/1.79
> are full-range fits **for the nonzero Pauli count**, notes the 3D fit rests on three sizes,
> and cross-references the $N\ge16$ value of 1.48. Verified: 1D 1.660 (full) / 1.475 ($N\ge16$),
> 2D 1.692, 3D 1.794.


**Severity:** internal inconsistency.

- §VI-D: *"the consistent-mass density and Pauli count grow sub-quadratically (empirical exponent ${\approx}1.4$–$1.5$)"*
- §VIII-B: *"a single power law to the full tested range gives ... (1D $\approx1.66$ ...)"*

Both refer to the 1D consistent-mass series. The difference is the fit range; neither section
says so. Recomputed:

| fit range | Pauli exponent | nnz exponent |
|---|---|---|
| $N\ge4$ | 1.660 | 1.681 |
| $N\ge16$ | 1.475 | 1.422 |
| $N\ge32$ | 1.342 | 1.232 |
| local $64\to128$ | 1.249 | 1.139 |

The exponent is falling monotonically — the same "density decline sets in" behaviour §VIII
attributes to spatial dimension is already visible in 1D and is what produces the discrepancy.

**Recommended action:** state the fit range at both sites, and add the local-exponent column.
The declining trend supports §VIII's argument rather than undermining it, so this is worth
surfacing rather than hiding.

---

### F6 — Minor factual, wording, and hygiene items

>### ✅ RESOLVED (2026-08-09)
>
> **F6.1** Table 2 → 1.058 and 1.025. **F6.2** "2–7 qubits" → "2–12". **F6.3** Remark 2 now
> opens with the genericity caveat and says "consistent with" rather than "predicts"; Table 3's
> caption separates expressibility-limited from optimizer-limited rows. **F6.4** $\tau=10^{-12}$
> stated at \eqref{eq:density}. **F6.5** Table 1 footnote distinguishes the band-aware/FWHT
> asymptote from the dense reference implementation. **F6.6** retired (see below). **F6.7** COI
> pluralised; README author block, VERIFICATION_REPORT (title, two authors, 249 words, 17 pages,
> email) corrected, with a scope banner; `PDF_PREFLIGHT.txt` page count fixed. **F6.8** all nine
> uncited entries placed at natural points — **33/33 cited, 0 missing**; README documents the
> six intentionally unused diagnostic figures.


**F6.1 — Table 2 transcription errors.** Recomputed Gershgorin ratios $s_G/\lambda_{\max}$:

| mass | $N=4$ | 8 | 16 | 32 | 64 | 128 |
|---|---|---|---|---|---|---|
| lumped | 1.110213 | 1.081462 | **1.053416** | 1.034542 | **1.022358** | **1.014447** |
| consistent | 1.105318 | 1.071128 | **1.058069** | 1.038417 | **1.025231** | **1.016850** |

Table 2 prints **1.056** for consistent/16 (should be 1.058) and **1.026** for consistent/64
(should be 1.025). The four lumped entries are correct.

**F6.2 — Stale qubit range.** §VI-C (Test C) says gradient variance is sampled at
"2--7 qubits"; the study and §VII-E use 2–12. `environment_versions.txt` records the extension
to 12 qubits, so §VI-C was simply not updated.

**F6.3 — Remark 2 over-claims.** *"This necessary condition predicts every row of Table 3."*
A parameter count establishes that the circuit map cannot be *surjective*; it cannot predict
that a *specific* target state is unreachable — the ground state could lie in the reachable
set. The empirical agreement is clean (and §1.3 above confirms it is not optimizer noise), but
the logic is a genericity argument. Suggest *"is consistent with every row of Table 3"* and add
one clause noting that a dimension count bounds the reachable manifold generically.

Also: the $L=3$, $N=16$ entry ($4.09\times10^{-9}$%) is optimizer-limited, not
expressibility-limited (true floor $\approx4\times10^{-12}$%), so Table 3 mixes two regimes.
Either rerun with more restarts or footnote it.

**F6.4 — Density threshold $\tau$ is never stated.** Eq. (9) defines $\delta_A(\tau)$
parametrically; every reported density uses $\tau=10^{-12}$ hardcoded
(`revision_study.py:357`, `mesh_dimension_study.py:136`). Tables 2 and 5 are not reproducible
without it. Add "$\tau=10^{-12}$" to both captions.

**F6.5 — Table 1 describes a different algorithm than the one shipped.** The $O(Nb^2)$ Cholesky,
$O(N^2)$ whitening, and $O(N^2\log N)$ Pauli rows assume band-aware triangular solves and a fast
Walsh–Hadamard transform. The code uses dense `scipy.linalg.cholesky` / `solve_triangular`
($O(N^3)$) and a per-string Python parity loop in `pauli_coefficients`. Since
`preprocessing_scaling.pdf` plots *those* wall times, the figure and the table measure different
algorithms. The text already hedges ("wall-clock times are machine-specific"), but one sentence
saying the table is the band-aware/FWHT cost while the profiling uses dense LAPACK would remove
the ambiguity.

**F6.6 — Linearization quoted outside its validity. — RETIRED.** As originally written: §IX-D
quoted *"a localization signal-to-noise ratio of only 0.13"*, but the delta-method variance behind
it linearizes $|\cdot|$ at its exact sign, which requires $|\text{mean}|\gg\text{sd}$ — i.e. exactly
not at SNR 0.13. The **F1** fix dissolves this: the SNR is now 1.44, where the linearization is
defensible. No action needed.

**F6.7 — Stale package metadata.** `README.md` and `VERIFICATION_REPORT.md` both state:
sole author (Bikalpa Gautam), no corresponding email, 217-word abstract, title
"Exact Deflation Conditions and Measurement Limits". The current `.tex` and compiled PDF have
two authors (Gautam, Subedi), `bgautam@ethz.ch`, a **247-word** abstract (249 after the F1 fix;
still within the 150–250 limit), and the title "Deflation and Measurement Limits". Also
`\section*{Conflict of Interest}` still reads *"The author declares"* (singular).

Additionally, **`PDF_PREFLIGHT.txt` reports "Pages: 13" and `VERIFICATION_REPORT.md` claims
"13 TQE-formatted pages", but the committed PDF is 16 pages.** Verified independently: the
committed PDF is 16 pages, and the unmodified `.tex` at `HEAD` compiles to 16 pages. This was
already stale before any edits in this review.

**F6.8 — Unused artefacts.** 9 of 33 bib entries are uncited (`Cerezo2021`, `Endo2023`,
`FarrarWorden2007`, `Kandala2017`, `Maeck2001`, `PeetersZ24_2001`, `Preskill2018`, `Tilly2022`,
`Zhang2022`) — harmless, but several are load-bearing for the framing and arguably *should* be
cited. Six generated figures are never included: `frequency_error`, `mac_heatmap`,
`damage_index`, `sensitivity_heatmap`, `pauli_spectrum`, `preprocessing_scaling`. Most of
`reproduce_study.py`'s shear-building output is now unused except via the Appendix C tolerance
table; consider saying so in the README so a referee does not go looking for the corresponding
manuscript section.

---

### F7 — Monte Carlo results depend on Pauli term ordering, not on the interpreter

>### ✅ RESOLVED (2026-08-09)
>
> Both `pauli_coefficients` and `pauli_decomposition` now sort with `kind="stable"` over the
> deterministic `itertools.product` generation order, fixing the term list and removing the
> variation entirely. This **re-based every Monte Carlo number**, which was the accepted cost:
> Table 4, the abstract's first-mode RMSE (2.23% → 2.45%), §VII-F, §VII-G's accuracies, and
> Remark 3 were all updated, and the whole package was regenerated in one pass so provenance
> is now single-environment.
>
> Decisive evidence for the mechanism came from `data/pauli_decomposition.csv`: across
> environments its **coefficients are bit-identical (max difference exactly 0.0)** and only the
> row order moves, diverging at row 6 among two groups of four tied magnitudes. No numerical
> difference exists to explain it away.
>
> The heavier substream-per-term redesign proved unnecessary. Two sources of variation remain
> and are documented rather than fixed: L-BFGS-B path dependence on different BLAS builds,
> which moves only the optimizer-limited near-zero depth-study entries, and machine-specific
> timings.

**Severity:** reproducibility. Partially mitigated; the root cause remains.
**Found:** while regenerating the F1 data. **Location:** `revision_study.py:166` (the sort),
`revision_study.py:240-250` (`sample_prepared`).

The archived Monte Carlo outputs did not reproduce on a clean run of the *same code* with the
*same seed* on Python 3.12: first-mode RMSE came back 2.3583 instead of 2.233719, with up to
9.2% deviation across the table, and `stylized_noise_budget` up to 34%.

The README attributed such drift to `scipy.optimize.minimize` shifting RNG consumption across
interpreters. **That explanation is wrong for these studies.** Nothing upstream of
`finite_shot_frequency` touches the module-level generator — `solve_vqd` creates its own local
`default_rng`, and `gradient_trainability` uses dedicated streams. Nor is the seed a variable:
numpy guarantees `Generator` reproduces the same stream for a given seed across versions and
platforms.

**The actual mechanism is term ordering.** `sample_prepared` walks the Pauli DataFrame in row
order, drawing one binomial per term from one shared generator, so *which term consumes which
draw* is determined by the row order. `pauli_coefficients` sorts by `|coefficient|` using
pandas' default **non-stable** quicksort, and the six-DOF Hamiltonian contains two groups of
**four exactly tied** magnitudes (0.06146276448005274 and 0.06702688171411184). Permuting only
those tied rows — bit-identical inputs, identical seed — moves the results:

| term order | mode1@$10^5$ | mode1@$10^4$ | mode4@$10^5$ |
|---|---|---|---|
| archive | 2.233719 | 7.869265 | 0.079209 |
| quicksort (current) | 2.358298 | 7.362970 | 0.076154 |
| `kind="stable"` | 2.358298 | 7.362970 | 0.076154 |
| total order (abs, then string) | 2.447855 | 7.878742 | 0.076814 |
| generation order | 2.361483 | 7.853067 | 0.080470 |

Tie-permutation alone shifts values by 0.8–11.8%, which spans the entire observed discrepancy.

**Honest limit of this finding:** none of the four orderings reproduces the archive's 2.233719.
So tie-ordering is demonstrably *sufficient* to produce discrepancies of this magnitude, but it
has not been shown to be the *sole* cause. Something else about the original environment also
differs. Note that the deterministic quantities — eigenvalues (0 ulps), Pauli coefficients,
densities, Gershgorin ratios, and the analytic shot requirements — are bit-identical here to the
archive, and are unaffected by any of this.

**What was done.** The corrected explanation is now in `README.md` and
`environment_versions.txt` (Note 3), and `validate_outputs.py` no longer asserts four-digit
agreement on stochastic quantities: `rmse_pct` uses `rtol=0.15` and the damage accuracy uses a
`[0.70, 0.95]` band. The deterministic assertions remain tight.

**What remains.** The estimator is still order-dependent. The proper fix is a separate RNG
substream per term (`default_rng(seed).spawn()`), making results invariant to term order *and*
to execution order — roughly five lines, but it re-bases every Monte Carlo number in the paper,
so it is deliberately left as its own decision. A cheaper partial hardening is to give the sort
a total order (`kind="stable"`, or sort by `|coefficient|` then Pauli string), which removes the
tie ambiguity without touching the shared-generator design.

---

## 3. Two results the paper undersells

### S1 — The lumped-mass Pauli count has an exact closed form

> **✅ ADDED as Proposition 4**, with a full proof. The hypotheses flagged as an overclaim risk
> are stated explicitly: power-of-two axis lengths, row-major power-of-two strides, diagonal
> mass, nearest-neighbour coupling. Without them the $\log_2 L_a$ XOR-class count fails.

Verified over $N=4$–128 in 1D, and against **every** lumped row of the 2D/3D Table 5:

$$N_P^{\text{lumped}} = (n_q+2)\,2^{\,n_q-1} = \Theta(N\log N)$$

| $n_q$ | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 10 |
|---|---|---|---|---|---|---|---|---|
| measured | 8 | 20 | 48 | 112 | 256 | 576 | 1280 | 6144 |
| $(n_q+2)2^{n_q-1}$ | 8 | 20 | 48 | 112 | 256 | 576 | 1280 | 6144 |

**Why it is dimension-independent.** For a Pauli string with X-mask $x$, the nonzero
coefficients live on the XOR-diagonal $\{(i,i\oplus x)\}$. A power-of-two grid with axis lengths
$L_d$ contributes $\log_2 L_d$ distinct XOR classes per axis, and $\sum_d\log_2 L_d=n_q$ always.
So there are exactly $n_q$ nonzero-$x$ classes regardless of $d$, each contributing $2^{n_q-1}$
surviving Z-patterns (symmetry kills half), plus $2^{n_q}$ from $x=0$:
$2^{n_q}+n_q2^{n_q-1}=(n_q+2)2^{n_q-1}$.

This explains the otherwise-suspicious fact that Table 5's lumped columns (48, 256, 576, 1280,
6144) are *identical* to the 1D values at matched $N$ — a referee will flag that as a copy-paste
error unless it is explained. It is also a much stronger statement than "lumped is sparser": it
proves the $O(N^2)$ Pauli bound is attained **only** through the mass model, which sharpens the
causal chain in Eq. (61) from an empirical observation to a proposition.

### S2 — The gradient decline is a normalization effect, not a trainability effect

> **✅ ADDED to §VII-E**, as a paragraph reporting the normalized statistic alongside the
> existing power-law-versus-exponential comparison.

Gradient variance falls 27× from 2 to 12 qubits — but the energy-landscape variance falls 21×
over the same range. The ratio is flat:

| $n_q$ | 2 | 4 | 6 | 8 | 10 | 12 |
|---|---|---|---|---|---|---|
| $\mathrm{Var}[g]$ | 2.62e-2 | 4.51e-3 | 2.22e-3 | 2.04e-3 | 1.05e-3 | 9.82e-4 |
| $\mathrm{Var}[E]$ | 4.94e-2 | 1.49e-2 | 5.20e-3 | 4.75e-3 | 3.23e-3 | 2.32e-3 |
| ratio | 0.530 | 0.302 | 0.426 | 0.430 | 0.325 | 0.423 |

Regression of $\log(\mathrm{Var}[g]/\mathrm{Var}[E])$ on $n_q$: slope $-0.036$/qubit,
95% CI $[-0.087,+0.016]$, $R^2=0.17$. **Relative to the objective's own scale, the gradients do
not shrink at all** over 2–12 qubits.

This is a considerably cleaner "no barren plateau over the tested range" statement than the
current power-law-vs-exponential AIC comparison, whose tail fit ($n_q\ge7$: slope $-0.118$,
95% CI $[-0.226,-0.010]$, $R^2=0.53$) is driven by a single 8-qubit outlier in a visibly flat
series and only marginally excludes zero. I reproduced every fit statistic the manuscript
quotes (full-range exponential $R^2=0.833$, power law $R^2=0.957$, tail slope $-0.118$), so the
existing text is accurate — it is just weaker than the data allow.

---

## 4. Reproducing the evidence

All checks below were run from the package root with the bundled `.venv`.

```bash
# Baseline integrity
.venv/bin/python validate_outputs.py
shasum -a 256 -c MANIFEST.sha256
```

**F1 (damage convention).** Duplicate `revision_study.damage_shot_study`, replacing the damaged
numerator operators with baseline ones:

```python
Ke_use = Ke0                      # instead of Ked
ed = rs.transformed_element_operators(Ke_use, Ld, sd, Hd.shape[0])
```

and rerun the 200-trial loop and the `minor_revision_study` delta-method block with the same
substitution.

**F2 (Davis–Kahan).** `revision_data/deflation_perturbation.csv`; compare
`davis_kahan_bound` against `2*operator_perturbation/spectral_separation` and test
`operator_perturbation < spectral_separation/2`.

**F5 / F6.1 (exponents, Gershgorin).**

```python
from revision_study import chain_model, mass_whiten, gersh, pauli_coefficients
# loop N in (4,8,16,32,64,128) x mass in (lumped, consistent);
# ratio = gersh(A)/lambda_max ; np.polyfit(log N, log N_P, 1)[0] over sub-ranges
```

**S1 (closed form).** Compare `len(pauli_coefficients(A/lam_max, tol=1e-10))` against
`(nq+2)*2**(nq-1)` for lumped mass in 1D, and against the lumped rows of
`dimension_generalization/data/dimension_scaling.csv`.

**S2 (gradient ratio).** `revision_data/gradient_trainability.csv`, uniform rows:
`np.polyfit(qubits, np.log(gradient_variance/energy_variance), 1)`.

---

## 5. Suggested priority order

All Round 1 items below are complete; the Round 2 items are tabulated above and are also
complete. What remains open is listed under "Still open" at the end.

| # | Item | Effort | Blocks submission? | State |
|---|---|---|---|---|
| 1 | **F1** — fix Eq. (49)–(50), rerun damage studies, update abstract/§VII-G/Table 8/conclusion | high | yes | ✅ done |
| 2 | **F2** — correct the Davis–Kahan numbers and precondition | low | yes | ✅ done |
| 3 | **F3** — ship the restart-study code | low | yes (policy) | ✅ done |
| 4 | **F6.1, F6.2, F6.7** — Table 2 ratios, qubit range, author/COI metadata, PDF page count | low | yes | ✅ done |
| 5 | **F5, F6.4** — state fit ranges and $\tau$ | low | no, but referees will ask | ✅ done |
| 6 | **F7** — deterministic term order (`kind="stable"`) | low | no | ✅ done |
| 7 | **F4** — make Algorithm 1 executable in the code | medium | no | ✅ done |
| 8 | **S1** — add the closed form as a proposition | low | no — free strengthening | ✅ done |
| 9 | **S2** — add the normalized-gradient statistic | low | no — free strengthening | ✅ done |
| 10 | **F6.3, F6.5, F6.8** — wording and hygiene | low | no | ✅ done |
| — | **F6.6** | — | — | ✅ retired by F1 |


---

## 6. Still open

Nothing in this list blocks submission, but each is a real limit on what the package
currently demonstrates.

1. **Proposition 4's genericity condition is stated, not proved.** A counterexample class
   (homogeneous meshes, $1.5N-1$ terms) is now given, but no sufficient condition for
   equality is established.
2. **The near-degenerate subspace diagnostic is untested.** The exactly-degenerate case is
   analytically trivial; the informative case is a near-degenerate pair under a physical
   perturbation.
3. **The trainability sweep runs at depth 2**, which \cref{rem:depth} proves cannot represent
   the target beyond small $n_q$, and differentiates only $\theta_0$ — the least
   plateau-prone parameter. The barren-plateau mechanism is structurally undetectable by
   this experiment.
4. **Contributions (iv) and (v) remain unevaluated.** They are now labelled as such rather
   than fixed.
5. **Test E still supplies exact damaged eigenvectors.** Disclosed in the Limitations and at
   Eq. (qmse), not resolved.
6. **Cross-platform determinism is untested.** All reproducibility claims here rest on one
   machine, one BLAS, one interpreter. The package's determinism also relies on pandas'
   descending-stable tie handling, an undocumented internal.
