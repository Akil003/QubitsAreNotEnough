# Possible Improvements

Items to review before final submission. Each entry references the specific location in the paper and proposes a concrete change.

---

## 1. Formal end-to-end complexity comparison is missing

**Location:** Absent from the paper. Natural home: Section VI or VIII.

**The issue:** The paper quantifies individual costs (Pauli terms, shots, Cholesky, optimizer evals) but never multiplies them together into a total cost expression compared against classical. The conclusion "classical is overwhelmingly preferable" is asserted but not derived.

**What would strengthen the paper:**

```
For k modes of an N-DOF banded chain:

  Classical:  O(k · N · √N)  =  O(kN^{3/2})    [Lanczos with gap ratio ~ π/(2N)]
  Quantum:    O(N³)                               [preprocessing alone]
            + O(k · N_P · ε⁻² · evals)           [measurement phase, N_P = O(N²) for consistent mass]
```

The key observation: **quantum preprocessing (O(N³)) already exceeds the entire classical solve (O(kN^{3/2})) for sparse 1D systems.** The quantum approach loses before any circuit runs. This is a logical impossibility — the pipeline requires classical work harder than the problem it's trying to solve quantumly.

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
