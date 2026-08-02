# Possible Improvements

Items to review before final submission. Each entry references the specific location in the paper and proposes a concrete change.

> **Status (all four implemented).**
> - **#1 Complexity comparison** — added subsection "End-to-end asymptotic cost" (`tab:complexity`, time + space) in Sec. VI and an illustrative 128-DOF table (`tab:concrete`) + hedged quantum-win remark in Sec. VIII. Numbers validated by `validate_outputs.py::validate_complexity`.
> - **#2 Causal chain** — `eq:trilemma` reframed from a trilemma to a two-term causal chain (Sec. VIII).
> - **#3 Gradient sweep** — extended 7 → 12 qubits with an exponential-vs-polynomial decay fit; power law (R²=0.96) beats exponential (R²=0.83), so no barren-plateau signature over the range. Regenerated `gradient_trainability.csv` + figure.
> - **#4 Depth law** — added a Remark after `tab:vqe_depth` deriving `n_q(L+1) ≥ 2^{n_q}−1 ⇒ L = Ω(N/log N)`, verified against every table row.

---

## 1. Formal end-to-end complexity comparison is missing

**Location:** Absent from the paper. Natural home: Section VI or VIII.

**The issue:** The paper quantifies individual costs (Pauli terms, shots, Cholesky, optimizer evals) but never multiplies them together into a total cost expression compared against classical. The conclusion "classical is overwhelmingly preferable" is asserted but not derived.

**What would strengthen the paper:**

The tested structures are **1D fixed-base shear chains** (Tests A–E, N ∈ {4,…,128}), so `K` and `M` are tridiagonal with half-bandwidth `b = 1` (constant). That controls the cost model:

```
For k modes of an N-DOF banded 1D chain (half-bandwidth b = O(1)):

  Classical (Lanczos/ARPACK):
    O(k · N · T_it)                    banded matvec O(N·b) = O(N); T_it iterations
      T_it ≈ 10²–10³ in practice (ARPACK shift-invert, ~const in N);
      O(√N)–O(N) for unpreconditioned Lanczos, gap-dependent (see caveat).

  Quantum (VQD) preprocessing:
    O(N·b²) = O(N)                     banded Cholesky of M   (NOT O(N³))
    + O(N²)                            whitening densifies A + Pauli decomposition
                                       (consistent mass: density → 0.31, N_P = O(N²))

  Quantum (VQD) measurement — the dominant term:
    O(k · N_P · ε⁻² · evals) = O(k · N² · ε⁻² · evals)
```

**Corrected key observation (the earlier O(N³) claim was wrong).** For these 1D chains the banded Cholesky is `O(N)`, not `O(N³)` — the "preprocessing loses before any circuit runs because O(N³) > O(kN^{3/2})" framing does **not** hold and must not go into the paper. The real, defensible asymmetry is **densification and measurement**, not factorization:

1. Classical Lanczos operates on the *sparse* `K` directly: cost scales with `nnz(K) = O(N)`.
2. Whitening turns `A = L⁻¹KL⁻ᵀ` into a *dense-ish* operator (density → 0.31 for consistent mass), so the Pauli representation has `N_P = O(N²)` terms.
3. **Merely enumerating and storing those `O(N²)` Pauli terms already costs more than the entire classical eigensolve** (`O(k·N·T_it)`) once `N` is large relative to `k` — and that is *before a single shot*. The `ε⁻²·evals` measurement factor then makes the gap enormous.

So the pipeline still "loses before any circuit runs," but because of `O(N²)` densification/Pauli construction, not `O(N³)` Cholesky.

**Caveat to resolve before quoting a classical exponent.** `T_it` depends on the relative spectral gap. For a uniform chain `λ_j = ω_j² ∝ j²/N²`, so the relative gap is `γ ~ 1/N²`; Kaniel–Paige–Saad gives `T_it ~ γ^{-1/2} ~ N`, i.e. `O(kN²)` total for unpreconditioned Lanczos — **not** the `O(kN^{3/2})` originally written (that assumed a `~1/N` gap, which is the ω-spacing, not the λ gap). ARPACK shift-invert side-steps this (≈ const iterations) at the price of one `O(N)` banded factorization. **Do not commit to a single crisp classical exponent in the paper**; state the solver (shift-invert vs plain Lanczos) and gap scaling explicitly, or report the empirical iteration count.

See [`complexity_section_addition.md`](complexity_section_addition.md) for draft tables (theoretical + practical).

---

## 2. "Three-way design problem" is really a two-way trade-off

**Location:** Section VIII (Discussion), subsection "The actual scaling bottleneck", Eq. `\label{eq:trilemma}` (line ~872 of `.tex`)

**Problem:** Sparsity of A is not an independent design variable — it is a deterministic consequence of the choice of M:

```
Choice of M  →  determines L  →  determines A = L⁻¹KL⁻ᵀ  →  determines Pauli count  →  determines measurement cost
```

There are only two things in tension:

1. **Physical accuracy of M** (consistent mass is more faithful than lumped)
2. **Measurement cost** (which increases because consistent M densifies A, increasing Pauli terms)

Classical Lanczos operates on K directly and never whitens, so the choice of M has no effect on classical solver cost.

**Suggested fix:** Reframe as a causal chain, not a trilemma:
```latex
The causal chain is:
\begin{equation}
\boxed{
\text{physical fidelity of }\M
\xrightarrow{\text{Cholesky}} \text{density of }\A
\xrightarrow{\text{Pauli decomp.}} \text{measurement cost}.
}
\end{equation}
Lumped mass may be cheaper quantum mechanically but can perturb high modes;
consistent mass may be more faithful but measurement-heavy.
```

---

## 3. Gradient/trainability study: 6 data points can't distinguish exponential from polynomial decay

**Location:** Section V (Eq. `\label{eq:gradvar}`); Section VII subsection "Gradient statistics raise, but do not settle, the barren-plateau question"

**What the paper already says (correctly):**
- "The present 2–7-qubit data cannot establish a structural barren plateau"
- "A finite $n_q \le 7$ study cannot prove or disprove an asymptotic barren plateau"

**The data:** ~16× drop in gradient variance over 5 qubit increments. Exponential decay would give ~32×. Six points on a log scale cannot distinguish exp(−αn) from 1/n² from a transient.

**Possible improvement:** Fit both models (log(V_g) vs n_q for exponential; log(V_g) vs log(n_q) for polynomial) and report R² values. This quantifies *how* inconclusive the data is rather than just asserting it. Alternatively, extend to 10–12 qubits (feasible with state-vector simulation at 2^12 = 4096 amplitudes) to get enough data points for a meaningful fit.

---

## 4. Parameter counting explains the depth study but is never stated

**Location:** Section VII, Table III (`tab:vqe_depth`)

**The issue:** Depth 2 works for 8-DOF but fails for 16-DOF. The paper shows this empirically but never gives the explanation:

- A real unit vector in D = 2^{n_q} dimensions has D−1 free parameters
- A depth-L real-amplitude circuit has n_q(L+1) parameters
- Necessary condition for expressibility: n_q(L+1) ≥ 2^{n_q} − 1

Checking against Table III:
- n_q=3, L=2: 9 params ≥ 7 free amplitudes → works ✓
- n_q=4, L=2: 12 params < 15 free amplitudes → fails ✗
- n_q=4, L=3: 16 params ≥ 15 free amplitudes → works ✓

This predicts every row. More importantly, it gives a **scaling law**: required depth L ≥ (2^{n_q} − 1)/n_q − 1, which grows exponentially in n_q. That's a one-line proof that circuit depth must grow exponentially with system size — destroying the logarithmic qubit advantage.

**Caveat:** Parameter count ≥ D−1 is necessary but not sufficient. The circuit topology (CNOT ring) may further restrict reachable states. But as a *lower bound* on required depth, it's rigorous and already explanatory for the observed data.

**Suggested addition:** A remark after Table III or in Section V.

---

# To-do (future revisions — not yet implemented)

Captured from a critical review; none block initial submission. `T1` is the highest-leverage item for first-round acceptance at IEEE TQE.

## T1. Add a 2D finite-element example (highest priority)

**Why:** every current experiment is a 1D shear chain — the *best* case for Cholesky fill. The regime that matters (2D/3D FE) is only extrapolated verbally, and the most probable first-round reviewer demand is "show at least one 2D example." A 2D case **strengthens** the negative conclusions (2D fill is worse); it does not change them.

**Blast radius (contained, additive):** the only new code is a 2D FE assembler (Q4 membrane grid or 2D spring lattice); the downstream pipeline (`mass_whiten`, `pauli_coefficients`, `qwc_groups`, density/Gershgorin) is dimension-agnostic and reused unchanged. Sweep e.g. `N ∈ {16, 64, 256}` and report density, Pauli terms, and QWC groups against the 1D chains. `tab:complexity` stays valid (already scoped to 1D); present 2D as the worse comparison and rewrite the Limitations "1D-only" paragraph from an apology into a result. Conclusions and thesis unchanged. **Effort ~half a day.** A symmetric mesh would additionally exercise repeated/degenerate modes.

## T2. Close the end-to-end noisy loop (and scale the simulator)
Run the full pipeline at 3–4 qubits optimizing **under shot noise** with deflation, so feasibility is demonstrated rather than inferred from the separated exact-state studies. Enabler: rewrite `RealAmplitudeCircuit` as a proper statevector simulator (tensor gate application, `O(2^n)` per gate instead of dense `2^n × 2^n` matmuls) so the loop — and the trainability/scaling sweeps — can reach 16–20 qubits.

## T3. Reduce cognitive load (presentation)
The paper is dense (74 numbered equations, plus propositions/definitions). Add a notation/symbol table, condense the repetitive parts of the Discussion, and consider a "validated vs proposed" contributions table (temporal tracking, the lexicographic refinement, and sensor-space recovery are proposals without experiments). *(The end-to-end workflow figure has now been added — `fig:workflow`.)*

---

*Deprioritized (revisit only if reviewers ask): a formal lower-bound proposition for the no-speedup claim; a fill-reducing DOF ordering to minimize the Pauli count; a classical-shadows measurement estimate; sharper deflation-theorem novelty positioning against classical (Wielandt) deflation.*
