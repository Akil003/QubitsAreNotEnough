# Complexity Section Addition (Draft)

## Placement

- **Theoretical table:** New subsection in Section VI (Measurement and Resource Accounting)
- **Practical table:** Section VIII (Discussion), under "The actual scaling bottleneck"

The paper already defines the resource vector $\mathcal{R}$ (Eq. 67) and error vector $\mathcal{E}$ (Eq. 66) but never fills them in as a direct classical-vs-quantum comparison.

---

## 1. Theoretical Complexity

| Step | Classical (Lanczos/ARPACK) | Quantum (VQD) |
|------|---------------------------|---------------|
| FE assembly | O(n_e · dof_e²) | Same (classical preprocessing) |
| Cholesky factorization | — | **O(N·b²) = O(N)** for 1D chains (half-bandwidth b = 1); O(N^{3/2}) for 2D nested dissection; O(N³) dense worst case |
| Whitening (L⁻¹KL⁻ᵀ) | — | **O(N²)** for 1D chains — L stays banded but A densifies (density → 0.31 consistent mass); O(N³) dense worst case |
| Pauli decomposition | — | N_P = O(4^{n_q}) = O(N²) terms; ~O(N² log N) to build via fast transform |
| Per eigenvalue: core operation | O(nnz(K) · iterations) = O(N · T_it) | O(N_P · shots · optimizer_evals) |
| Lanczos iterations T_it | ≈ const with ARPACK shift-invert; O(√N)–O(N) unpreconditioned (gap-dependent, see below) | N/A |
| Optimizer iterations | N/A | O(hundreds); no convergence guarantee |
| Shots per optimizer eval | N/A | O(N_P / ε²) for ε energy accuracy |
| **Total for k modes** | **O(k · nnz(K) · iterations)** | **O(k · N_P · shots · optimizer_evals)** |
| Space | O(N·k) | O(2^{n_q}) statevector + O(N²) classical |

### Assumptions to state explicitly

- Classical: sparse K with bandwidth b, k ≪ N, well-separated low modes
- Quantum: exact state-vector simulation (no hardware noise), independent Pauli measurement
- "iterations" in Lanczos ≈ number of matrix-vector products until Ritz values converge to tolerance

### What "iterations" means in Lanczos

Each iteration = one matrix-vector product Kv, extending the Krylov subspace by one dimension. Convergence rate for eigenvalue j depends on the spectral gap ratio:

```
r_j = (λ_{j+1} - λ_j) / (λ_max - λ_j)
```

- Large gap → fast convergence (~10–30 iterations)
- Clustered eigenvalues → slow convergence (hundreds+)
- ARPACK uses implicit restarts: typically converges in ~100–300 total matvecs regardless of N

For uniform n-DOF shear chains, the **frequency (ω) spacing** of low modes ≈ π/(2n). The table below is indexed by that ω-spacing:

| n (DOF) | ω-spacing (low modes) | Typical iterations (first k modes) |
|-------------|----------------------|-------------------------------------|
| 6           | ~0.26                | ~10–20                              |
| 20          | ~0.08                | ~30–60                              |
| 100         | ~0.016               | ~150–300                            |
| 1000        | ~0.0016              | ~1500+                              |

> **Correction / caveat.** Lanczos runs on `A` (eigenvalues `λ = ω²`), so convergence is governed by the **λ relative gap**, not the ω spacing. For a uniform chain `λ_j ∝ j²/n²`, so the relative gap is `γ ~ 1/n²` — an order tighter than the `~1/n` ω-spacing tabulated above. Kaniel–Paige–Saad then gives `T_it ~ γ^{-1/2} ~ n` for unpreconditioned Lanczos, i.e. `O(k·N²)` total — **not** `O(k·N^{3/2})`. Shift-invert ARPACK removes this dependence (≈ const iterations) at the price of one banded factorization, `O(N)`. **The exact classical exponent is model-dependent (const → √N → N iterations); do not quote a single crisp value without stating the solver and gap assumption.**

Even under the most pessimistic estimate above, classical total cost at n=1000 is ≈ O(k · N · T_it) ~ 10·1000·1500 ≈ 1.5×10⁷ flops (sub-second). The quantum pipeline at the same N would need ~10⁵ Pauli terms × 10⁵ shots × hundreds of optimizer evals ≈ 10¹² circuit executions. The comparison is illustrative — the point is the ~5-order-of-magnitude gap, not the precise classical exponent.

---

## 2. Practical Complexity

| Question | Classical (Lanczos) | Quantum (VQD) |
|----------|---------------------|---------------|
| State preparation cost | N/A (no state to prepare) | Ansatz-dependent; shallow circuits fail for N>8 (Table III in paper) |
| Dominant runtime today | Sparse matvec in Lanczos (~μs each) | Shot acquisition (~ms per circuit × 10⁵ shots × N_P terms × evals) |
| What scales poorly? | Nothing problematic for k≤10, N<10⁶ | Pauli count: 576→2950 going lumped→consistent at N=128; shot cost for low-frequency modes |
| Measurements required | 0 (deterministic) | 1.9×10⁶ shots for 6-DOF, 4 modes at 10⁵/term — still 2.23% RMSE on mode 1 |
| Damage detection feasibility | Trivial once modes are known | 26% correct at 10⁵ shots (vs 16.7% random guess) — barely above chance |
| Preprocessing overhead | Negligible (assembly only) | Banded Cholesky O(N) + whitening/densification O(N²) + Pauli build O(N² log N) — dominated by O(N²), not O(N³), for 1D chains |
| Would you choose this over Lanczos? | — | **No.** For any N where Cholesky is feasible, classical wins by orders of magnitude |

---

## 3. Concrete Example: 128-DOF Chain, k=4 Modes

| Resource | Classical | Quantum (consistent mass) |
|----------|-----------|---------------------------|
| Qubits | — | 7 |
| Classical preprocessing | Assembly only: O(N) | Banded Cholesky O(N) + whitening + Pauli decomposition, O(N²)-dominated (not O(N³)) |
| Nonzero Pauli terms | — | 2950 |
| Commuting groups | — | 221 (greedy QWC) |
| Shots per mode (for ~1% freq accuracy) | — | ~10⁵ per term × 2950 terms = 2.95×10⁸ per optimizer eval |
| Optimizer evals per mode | — | ~200–500 (L-BFGS-B) |
| Total circuit executions (4 modes) | — | ~10¹¹ |
| Lanczos matvecs (4 modes) | ~1400 (measured, unpreconditioned eigsh; ≈tens with shift-invert) | — |
| Lanczos wall time (sparse, N=128) | <1 second (~30 ms measured) | — |
| Quantum wall time (1μs/circuit) | — | ~10⁵ seconds ≈ 28 hours |

---

## 4. The Fundamental Asymmetry

Classical Lanczos exploits sparsity: cost is O(nnz(K)) per iteration, and nnz(K) = O(N) for 1D chains, O(N^{4/3}) for 3D FE meshes.

Quantum VQD destroys sparsity: Cholesky whitening densifies the operator, Pauli decomposition grows as O(density × 4^{n_q}), and every Pauli term must be measured independently.

The quantum approach pays O(N²) classical preprocessing (the banded Cholesky itself is only O(N); the O(N²) comes from whitening-induced densification and Pauli enumeration) to create a problem that is then *harder* to solve quantumly than the original sparse, O(N)-per-iteration problem was classically. Simply writing down the O(N²) Pauli operator already exceeds the classical eigensolve for large N.

---

## 5. When Could Quantum Possibly Win?

For completeness, identify the regime where quantum might become competitive:

1. **N so large that even sparse Lanczos is slow** — but the quantum route needs the O(N²) dense-ish whitened operator and its Pauli list, which becomes infeasible to store/measure well before sparse Lanczos does (the banded Cholesky itself stays O(N))
2. **Hardware with negligible shot cost** — requires fault-tolerant quantum phase estimation, not VQE
3. **Problems where sparsity doesn't help classically** — dense matrices from integral equations, not typical FE
4. **Quantum phase estimation (QPE) instead of VQE** — deterministic, no optimizer, but requires fault-tolerant hardware with millions of physical qubits

The paper's conclusion holds: for structural FE eigenproblems, no VQE-based approach is competitive with classical methods, and the bottlenecks are fundamental (preprocessing, measurement, expressibility), not engineering limitations to be optimized away.

---

## Notes for LaTeX Formatting

- Theoretical table → `\begin{table}` in Section VI
- Practical table → `\begin{table*}` (full-width) in Section VIII
- Consider a consolidated "end-to-end" comparison figure showing the pipeline steps with costs annotated
- Reference existing equations: resource vector $\mathcal{R}$ (Eq. 67), error vector $\mathcal{E}$ (Eq. 66), Pauli variance (Eq. 63), frequency variance (Eq. 65)
