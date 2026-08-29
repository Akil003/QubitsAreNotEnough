# Paper Explained: A Complete Guide for CS Readers

## Who This Is For

You have a CS background (algorithms, data structures, complexity theory, maybe some ML). You may not know anything about quantum computing, structural engineering, or physics. This document explains every concept from first principles, with worked examples and proofs where relevant.

---

> **Note — kept in sync with the current paper.** This guide was written for an earlier draft. Two corrections are reflected below: (a) the **preprocessing cost** — the banded Cholesky of a 1D chain is **O(N)**, *not* O(N³) (that dense figure applies only to non-banded problems); the real bottleneck is whitening-induced **densification → O(N²) Pauli terms → measurement**. (b) The **trainability study** now spans **2–12 qubits** and finds *polynomial*, not exponential, gradient decay (no barren-plateau signature over the tested range). The end-to-end complexity is now tabulated in the paper (Section VI, "Table 1").

## Part 0: What Is This Paper About? (The 30-Second Version)

**The question:** Can a quantum computer solve a specific math problem (finding eigenvalues of large matrices) faster than a classical computer?

**The answer:** No, not currently. And this paper tells you exactly why — with theorems, experiments, and precise resource accounting.

**The math problem:** Given two large matrices K and M (both N×N, sparse, symmetric), find the smallest few values λ and corresponding vectors φ such that Kφ = λMφ. This is called a "generalized eigenvalue problem."

**Why it matters in the real world:** Those eigenvalues tell you the natural vibration frequencies of a building, bridge, or aircraft wing. If a frequency changes, something might be damaged. Engineers solve this problem constantly.

**Why people thought quantum might help:** A quantum computer stores an N-dimensional vector using only log₂(N) "qubits" (quantum bits). So a 1,000,000-dimensional vector needs only 20 qubits. That sounds like an exponential advantage.

**Why it doesn't actually help (the paper's conclusion):** Storing the vector cheaply is just one step. You also need to:
1. Preprocess the input classically (Cholesky whitening + Pauli decomposition) — whitening destroys the sparsity that makes classical solvers fast, so just *building* the operator you will measure can cost more than solving the whole problem classically
2. Have a deep enough quantum circuit (grows with N, not log N)
3. Measure the answer (costs O(N²) measurements in the worst case)
4. Repeat measurements billions of times to get useful accuracy

The net effect: the quantum approach is far worse than classical algorithms for this problem today.

---

## Part 1: Background Concepts (No Prior Knowledge Assumed)

### 1.1 What Is an Eigenvalue Problem?

**Definition:** Given a square matrix A (size N×N), an eigenvalue λ and eigenvector v satisfy:

```
Av = λv
```

This means: when you multiply A by v, you get back the same vector v, just scaled by the number λ.

**Concrete example:** Let A = [[2, 1], [1, 2]].

Try v = [1, 1]:
- Av = [2·1 + 1·1, 1·1 + 2·1] = [3, 3] = 3·[1, 1]
- So λ = 3, v = [1, 1] is an eigenpair.

Try v = [1, -1]:
- Av = [2·1 + 1·(-1), 1·1 + 2·(-1)] = [1, -1] = 1·[1, -1]
- So λ = 1, v = [1, -1] is another eigenpair.

A 2×2 matrix has (at most) 2 eigenvalues. An N×N matrix has (at most) N eigenvalues.

**Why eigenvalues matter:** They reveal fundamental properties of the system the matrix represents. For a vibrating structure, each eigenvalue corresponds to a natural vibration frequency.

### 1.2 What Is a Generalized Eigenvalue Problem?

Instead of Av = λv, we have:

```
Kv = λMv
```

Two matrices now. This is harder because M "weights" the right side. You can think of it as: "find directions where K and M are proportional."

**Why two matrices?** In structural engineering:
- K encodes stiffness (how resistant the structure is to deformation)
- M encodes mass (how heavy each part is)
- The eigenvalues λ give ω² (squared angular frequency of vibration)
- The frequency in Hz is f = √λ / (2π)

### 1.3 What Is a Sparse Matrix?

A matrix is "sparse" if most of its entries are zero. For example, a 1000×1000 matrix where each row has at most 5 nonzero entries has 1,000,000 cells but only ~5,000 nonzeros. That's 0.5% full.

**Why sparsity matters for algorithms:** If a matrix is sparse, you can multiply it by a vector in O(nnz) time (where nnz = number of nonzeros), not O(N²). Classical eigensolvers exploit this heavily.

**Structural matrices are sparse because:** Each element (beam, column) in a building only connects to its neighbors. A node on floor 3 doesn't directly interact with a node on floor 47. So K and M have nonzeros only where elements are physically connected — typically O(N) nonzeros total for a 1D structure.

**Is K an adjacency matrix?** Not exactly, but it's closely related. Here's the distinction:

- An **adjacency matrix** has entry A_ij = 1 if nodes i and j are connected, 0 otherwise. It only encodes topology (who connects to whom).
- The **stiffness matrix K** has the same sparsity pattern as the adjacency matrix (nonzero in the same positions), but the values encode *how strongly* nodes are connected — the physical stiffness of the beam or spring between them.

So K is like a *weighted* adjacency matrix where the weights represent physical stiffness. Similarly, M (the mass matrix) has a similar sparsity pattern but encodes how mass is distributed.

**Graph Laplacian analogy:** If you know graph theory, K is actually much closer to a **graph Laplacian** than a plain adjacency matrix:
- Graph Laplacian: L = D - A (degree matrix minus adjacency matrix)
- Stiffness matrix: K_ii = sum of stiffnesses of all springs attached to node i; K_ij = negative stiffness of spring between i and j

The structure is identical: positive diagonal, negative off-diagonals, rows sum to zero (for a free-free structure). The eigenvalues of a graph Laplacian give connectivity information; the eigenvalues of K (weighted by M) give vibration frequencies. Same math, different physical interpretation.

**Concrete example (3-node chain with springs):**

```
    k₁      k₂      k₃
[wall]---[node 1]---[node 2]---[node 3]
          mass m₁    mass m₂    mass m₃
```

Stiffness matrix (after fixing the wall):
```
K = [ k₁+k₂   -k₂      0   ]
    [ -k₂    k₂+k₃   -k₃   ]
    [  0      -k₃      k₃   ]
```

This is tridiagonal (bandwidth = 1) because each node only connects to its immediate neighbors. Entry K_ij is nonzero only if there's a spring between nodes i and j. The diagonal K_ii is the total stiffness attached to node i. This is exactly the graph Laplacian structure with physical stiffness values as weights.

**Now the mass matrix M — two versions:**

The mass matrix encodes how mass is distributed across the structure. There are two standard ways to build it:

**Version 1: Lumped mass (diagonal)**

The simplest approach — concentrate all mass at the nodes:

```
M_lumped = [ m₁   0    0  ]
           [  0   m₂   0  ]
           [  0    0   m₃ ]
```

Each diagonal entry is just the mass of that floor. Off-diagonals are zero. This says: "floor 1 has mass m₁, floor 2 has mass m₂, etc. The masses don't interact."

**Physical meaning:** When floor 2 accelerates, only m₂ resists (F = ma). The mass of floor 1 doesn't care what floor 2 is doing. This is an approximation — it pretends all mass lives exactly at the nodes.

**Why CS people should care:** Lumped mass is diagonal, so L = √M is diagonal, so L⁻¹ is diagonal, so the whitened matrix A = L⁻¹KL⁻ᵀ preserves the sparsity of K. This keeps the matrix sparse, which directly reduces the quantum measurement cost (explained in detail in Part 6 — for now, just know that a sparser matrix = fewer measurements needed on the quantum computer).

**Version 2: Consistent mass (tridiagonal for 1D)**

A more physically accurate model that accounts for how mass is distributed *along* each beam/spring element, not just at the endpoints:

```
M_consistent = [ (m₁+m₂)/3    m₂/6        0       ]
               [   m₂/6     (m₂+m₃)/3    m₃/6     ]
               [    0         m₃/6      m₃/3      ]
```

(Exact coefficients depend on the element formulation; the 1/3 and 1/6 come from integrating shape functions over the element length.)

**Physical meaning of off-diagonal entries:** M_ij ≠ 0 means "when node j accelerates, node i also feels an inertial force." This happens because the beam between them has distributed mass — if you shake one end of a heavy beam, the other end feels it through the beam's own inertia.

**Analogy:** Think of two people holding opposite ends of a heavy rope. If person A jerks their end, person B feels a tug — not because of a spring force, but because the rope's mass creates inertial coupling. That coupling is what the off-diagonal entries in M encode.

**Why this matters for the paper:** Consistent mass is more accurate physics (especially for higher-frequency modes), but it makes M non-diagonal. When you Cholesky-factor a non-diagonal M and compute A = L⁻¹KL⁻ᵀ, the inverse L⁻¹ is dense, which makes A dense, which means more quantum measurements are needed (the full explanation of *why* denser matrices cost more to measure is in Part 6, Section 6.1).

**Summary of the two mass models:**

| Property | Lumped mass | Consistent mass |
|----------|------------|----------------|
| Structure | Diagonal | Tridiagonal (1D) or banded |
| Physical accuracy | Good for low modes, worse for high modes | Accurate for all modes |
| Cholesky factor L | Diagonal | Banded (but L⁻¹ is dense) |
| Whitened A = L⁻¹KL⁻ᵀ | Sparse (same pattern as K) | Dense (fills in) |
| Quantum measurement cost | Low | 5× higher |

### 1.4 How Do Classical Algorithms Solve This?

The Lanczos algorithm (and its relatives like ARPACK) find the k smallest eigenvalues of a sparse N×N system in approximately:

```
Time:  O(k · nnz(K) · iterations)
Space: O(N · k)
```

**What is `nnz`?** It stands for "number of nonzeros" — the count of entries in the matrix that are not zero. For a sparse N×N matrix, nnz << N². For a 1D structural chain, each node connects to at most 2 neighbors, so each row has at most ~3 nonzero entries, giving nnz = O(N). For a 2D mesh, nnz = O(N) as well (each node has a bounded number of neighbors). For a dense matrix, nnz = N².

**Why does nnz matter?** The core operation in Lanczos is matrix-vector multiplication: compute Kv for some vector v. For a dense matrix this costs O(N²). For a sparse matrix it costs O(nnz) because you only touch the nonzero entries. Since Lanczos does many such multiplications, the total cost is proportional to nnz, not N².

**What are "iterations" and why do we iterate?**

The Lanczos algorithm is *iterative* — it doesn't compute eigenvalues in one shot. Here's the intuition:

1. **Start with a random vector** v₁.
2. **Multiply by the matrix:** compute Kv₁. This gives a new vector that "leans toward" the dominant eigenvector.
3. **Orthogonalize** against previous vectors (so you don't just converge to the same eigenvector repeatedly).
4. **Repeat:** compute Kv₂, orthogonalize, get v₃, etc.

After t iterations, you have a set of t vectors {v₁, v₂, ..., vₜ} that span a "Krylov subspace." You then solve a small t×t eigenvalue problem (cheap, since t << N) to get approximate eigenvalues of the original N×N matrix.

**Why does this converge?** Each multiplication by K amplifies components along the dominant eigenvectors. After enough iterations, the subspace contains excellent approximations to the eigenvectors you want. The number of iterations needed depends on the "spectral gap" (how well-separated the eigenvalues are). For well-separated eigenvalues, convergence is fast (maybe 20-50 iterations). For clustered eigenvalues, you might need hundreds.

**Concrete example:** For a 1,000,000×1,000,000 sparse matrix with nnz = 5,000,000, finding k=10 eigenvalues might take 50 iterations:
- Each iteration: one matrix-vector multiply = O(5,000,000) operations
- Total: 50 × 5,000,000 = 250,000,000 operations
- On a modern laptop at ~10⁹ ops/sec: about 0.25 seconds

This is why classical eigensolvers are so hard to beat — they exploit sparsity directly and converge in very few iterations for well-separated eigenvalues.

For a sparse system where nnz = O(N), the total cost is roughly O(N·k·iterations). These algorithms are 40+ years old, extremely well-optimized, and run on laptops for N = 1,000,000.

### 1.5 What Is a Qubit?

A classical bit is 0 or 1. A qubit is a unit vector in a 2-dimensional complex vector space. You can write it as:

```
|ψ⟩ = α|0⟩ + β|1⟩
```

where α and β are complex numbers with |α|² + |β|² = 1.

**n qubits together** represent a unit vector in a 2ⁿ-dimensional space. So:
- 3 qubits → 8-dimensional vector
- 10 qubits → 1024-dimensional vector
- 20 qubits → 1,048,576-dimensional vector

This is the "exponential compression" that makes quantum computing exciting: you store an N-dimensional vector using only log₂(N) qubits.

### 1.6 What Is VQE (Variational Quantum Eigensolver)?

VQE is a quantum algorithm for finding the smallest eigenvalue of a matrix. Here's how it works:

1. **Parameterize a quantum state:** Build a circuit with adjustable parameters θ = (θ₁, θ₂, ..., θₚ). The circuit produces a state |ψ(θ)⟩.

2. **Measure the energy:** Compute E(θ) = ⟨ψ(θ)|H|ψ(θ)⟩ (the "expectation value" of matrix H in state ψ). This is always ≥ the smallest eigenvalue (by the variational principle).

3. **Optimize:** Use a classical optimizer (like gradient descent) to minimize E(θ).

4. **Converge:** When E(θ) stops decreasing, you've (hopefully) found the smallest eigenvalue and its eigenvector.

**CS analogy:** This is exactly like training a neural network. The circuit is the model, θ are the weights, E(θ) is the loss function, and you use gradient descent to minimize it. The key difference: the "model" is a quantum circuit that can represent 2ⁿ-dimensional vectors with only n·(L+1) parameters (where L is circuit depth).

### 1.7 What Is VQD (Variational Quantum Deflation)?

VQE finds eigenvalue #1 (the smallest). What about #2, #3, etc.?

**Deflation idea:** After finding eigenvector u₁, modify the matrix:

```
H_new = H + β · |u₁⟩⟨u₁|
```

This adds β to eigenvalue #1, pushing it above eigenvalue #2. Now eigenvalue #2 becomes the smallest, and VQE finds it.

**The notation |u₁⟩⟨u₁|:** This is a "projector" — a matrix that, when multiplied by any vector v, gives (u₁ᵀv)·u₁. It projects onto the direction of u₁. Adding β times this projector to H increases only the eigenvalue associated with u₁.

---

## Part 2: The Paper's Pipeline (What It Actually Does)

### 2.1 The Full Algorithm, Step by Step

```
INPUT:  M (mass matrix, N×N, sparse, symmetric positive definite)
        K (stiffness matrix, N×N, sparse, symmetric)
        m (number of eigenvalues to find)

Step 1: CHOLESKY FACTORIZATION
        Compute M = LLᵀ (L is lower-triangular)
        Cost: O(N) for 1D banded chains; O(N^{3/2}) for 2D sparse meshes (nested dissection); O(N³) dense worst case

Step 2: WHITENING (transform to standard eigenvalue problem)
        Compute A = L⁻¹ K L⁻ᵀ (via triangular solves, not explicit inverse)
        Now Ay = λy is equivalent to Kφ = λMφ where y = Lᵀφ
        WARNING: A can be DENSE even when K and M are sparse

Step 3: NORMALIZATION
        Compute Gershgorin upper bound: U_G = max_i(A_ii + Σ_{j≠i} |A_ij|)
        Choose optimization scale s_t (e.g., from a few Lanczos steps)
        Set Ã = A / s_t

Step 4: PADDING (make dimensions a power of 2)
        If N=6, need 3 qubits, so pad to 8×8:
        H = [Ã    0  ]
            [0   α·I₂]   where α = U_G/s_t + δ (δ > 0)
        The padding eigenvalues (= α) are guaranteed > all physical eigenvalues

Step 5: PAULI DECOMPOSITION
        Write H = Σ hₗ Pₗ where each Pₗ is a tensor product of Pauli matrices
        Count the nonzero terms N_P (this determines measurement cost)

Step 6: FOR EACH MODE r = 1, 2, ..., m:
        a. Set penalty: β_j = U_r - ε̂_j + δ_j + τ  for all j < r
        b. Form deflated Hamiltonian: H_r = H + Σ β_j |û_j⟩⟨û_j|
        c. Run VQE with multiple random restarts
        d. Check all diagnostics (frequency error, residual, variance, etc.)
        e. Recover physical mode: φ = L⁻ᵀ y (truncate padding, normalize)

OUTPUT: Eigenvalues λ₁, ..., λₘ and eigenvectors φ₁, ..., φₘ
        Frequencies f_i = √λ_i / (2π)
```

### 2.2 Why Each Step Exists

| Step | Why it's needed | What goes wrong without it |
|------|----------------|---------------------------|
| Cholesky | VQE needs ONE matrix, not two | Can't formulate the optimization problem |
| Whitening | Makes quantum orthogonality = physical orthogonality | Modes aren't physically meaningful |
| Normalization | Eigenvalues should be in [0,1] for good optimization | Gradients are tiny or huge; optimizer fails |
| Padding | Qubit count must be integer; 2ⁿ dimensions | Can't map to a quantum circuit |
| Pauli decomposition | Quantum computers measure in the Pauli basis | Can't estimate the energy |
| Deflation | VQE only finds the minimum; need higher modes too | Can only get eigenvalue #1 |

---

## Part 3: The Main Theorem (Deflation Threshold)

### 3.1 The Problem

You've found mode 1 with eigenvalue ε₁. You want mode 2 (eigenvalue ε₂ > ε₁). You add a penalty:

```
H_deflated = H + β · |u₁⟩⟨u₁|
```

**Question:** How large must β be?

### 3.2 Theorem 1 (Necessary and Sufficient Condition)

**Statement:** Mode r is the unique ground state (minimum eigenvalue) of the deflated operator if and only if:

```
β_j > ε_r - ε_j    for ALL j < r
```

**In words:** The penalty for each previously-found mode j must be larger than the gap between the target eigenvalue εᵣ and that mode's eigenvalue εⱼ.

### 3.3 Proof (Complete)

**Setup:** The original matrix H has eigenvectors |u₁⟩, |u₂⟩, ..., |uₙ⟩ with eigenvalues ε₁ ≤ ε₂ ≤ ... ≤ εₙ. The deflated matrix is:

```
H_r = H + Σⱼ₌₁ʳ⁻¹ βⱼ |uⱼ⟩⟨uⱼ|
```

**Key observation:** The eigenvectors don't change! Since |uⱼ⟩ are already eigenvectors of H, and |uⱼ⟩⟨uⱼ| only acts on |uⱼ⟩ itself:

- H_r |uₖ⟩ = H|uₖ⟩ + Σⱼ βⱼ |uⱼ⟩⟨uⱼ|uₖ⟩

The eigenvectors are orthonormal, so ⟨uⱼ|uₖ⟩ = 1 if j=k, else 0. Therefore:

- For k < r: H_r |uₖ⟩ = εₖ|uₖ⟩ + βₖ|uₖ⟩ = (εₖ + βₖ)|uₖ⟩
- For k ≥ r: H_r |uₖ⟩ = εₖ|uₖ⟩ (unchanged)

**The eigenvalues of H_r are:**
- Mode k < r: eigenvalue = εₖ + βₖ (boosted)
- Mode k ≥ r: eigenvalue = εₖ (unchanged)
- Padding states: eigenvalue = α (unchanged, above all physical)

**For |uᵣ⟩ to be the unique minimum:**
- Need εₖ + βₖ > εᵣ for all k < r → βₖ > εᵣ - εₖ ✓
- Need εₖ > εᵣ for all k > r → true by assumption (eigenvalues are ordered)
- Need α > εᵣ → true by padding construction

**Why "if and only if":** If any βⱼ ≤ εᵣ - εⱼ, then εⱼ + βⱼ ≤ εᵣ, meaning mode j has an eigenvalue ≤ εᵣ in the deflated operator. Then |uᵣ⟩ is NOT the unique minimum. □

### 3.4 Worked Example

Consider a 3-mode system with eigenvalues ε₁ = 0.1, ε₂ = 0.3, ε₃ = 0.7.

**Finding mode 2 (r=2):**
- Need β₁ > ε₂ - ε₁ = 0.3 - 0.1 = 0.2
- Any β₁ > 0.2 works. E.g., β₁ = 0.25.
- After deflation: mode 1 has eigenvalue 0.1 + 0.25 = 0.35 > 0.3 ✓
- Mode 2 still has eigenvalue 0.3 → it's the minimum ✓

**Finding mode 3 (r=3):**
- Need β₁ > ε₃ - ε₁ = 0.7 - 0.1 = 0.6
- Need β₂ > ε₃ - ε₂ = 0.7 - 0.3 = 0.4
- E.g., β₁ = 0.65, β₂ = 0.45.
- After deflation: mode 1 → 0.75, mode 2 → 0.75, mode 3 → 0.7 (minimum) ✓

### 3.5 Why You Can't Just Set β = ∞

**The problem with approximate projectors:** In practice, you don't know the exact eigenvectors u₁, u₂, etc. You have approximations û₁, û₂. The actual deflated operator is:

```
Ĥ_r = H + Σⱼ βⱼ |ûⱼ⟩⟨ûⱼ|
```

The error (perturbation) is:

```
‖Δ‖ = ‖Ĥ_r - H_r‖ ≤ Σⱼ βⱼ · ‖|ûⱼ⟩⟨ûⱼ| - |uⱼ⟩⟨uⱼ|‖ = Σⱼ βⱼ · sin(angle(ûⱼ, uⱼ))
```

**The perturbation is proportional to β.** Larger penalty → larger error amplification.

**Davis-Kahan theorem** then bounds the error in the recovered mode:

```
sin(angle(recovered mode, true mode)) ≤ 2‖Δ‖ / gap
```

where "gap" is the distance from the target eigenvalue to the nearest other eigenvalue in the deflated spectrum.

**CS analogy:** This is exactly like the condition number in numerical linear algebra. A system Ax = b with condition number κ amplifies input errors by factor κ. Here, β plays the role of the condition number — it amplifies projector errors. The optimal strategy: use the smallest β that exceeds the threshold.

**Another analogy:** Learning rate in SGD. Too small → doesn't converge (below threshold = wrong mode found). Too large → overshoots (amplifies noise from imperfect prior modes). The sweet spot is determined by the spectral gap.

### 3.6 Experimental Verification

The paper sweeps β from 0 to 1.25× the threshold for modes 2, 3, and 4. The result is a perfect step function:
- β/threshold < 1.0 → overlap with target mode = 0 (wrong mode found)
- β/threshold > 1.0 → overlap with target mode = 1 (correct mode found)
- Transition happens at exactly β/threshold = 1.0

This confirms the theorem is sharp (tight) — not just sufficient but also necessary.

---

## Part 4: The Preprocessing Problem (Why Sparsity Dies)

### 4.1 The Cholesky Whitening Step

**Goal:** Convert Kφ = λMφ into Ay = λy (standard form that VQE can handle).

**Method:**
1. Factor M = LLᵀ (Cholesky decomposition — L is lower-triangular)
2. Set A = L⁻¹ K L⁻ᵀ
3. Now Ay = λy with y = Lᵀφ

**Why this works (Proposition 1, proved in the paper):**

The key property is that the Euclidean dot product of the transformed vectors equals the M-weighted dot product of the original vectors:

```
yᵢᵀyⱼ = (Lᵀφᵢ)ᵀ(Lᵀφⱼ) = φᵢᵀLLᵀφⱼ = φᵢᵀMφⱼ
```

**Why this matters:** On a quantum computer, two states are "orthogonal" when their Euclidean dot product is zero. After whitening, quantum orthogonality = physical mass-orthogonality. The quantum computer's natural notion of "perpendicular" matches the physics.

### 4.2 The Sparsity Destruction Problem

**The critical issue:** Even when K and M are both sparse, A = L⁻¹KL⁻ᵀ can be DENSE.

**Why?** The inverse of a sparse matrix is generally dense. Think of it this way: if L is a lower-triangular matrix with a band structure, L⁻¹ has nonzeros everywhere below the diagonal because back-substitution propagates values.

**Concrete example from the paper (N=128):**

| Mass model | K sparsity | M sparsity | A density | Pauli terms |
|-----------|-----------|-----------|-----------|-------------|
| Lumped (M = diagonal) | tridiagonal | diagonal | 2.3% | 576 |
| Consistent (M = tridiagonal) | tridiagonal | tridiagonal | 30.8% | 2,950 |

**Why lumped mass preserves sparsity:** When M is diagonal, L = √M is also diagonal. Then L⁻¹ is diagonal. So A = L⁻¹KL⁻ᵀ = D⁻¹KD⁻¹ where D is diagonal — this just rescales rows and columns of K, preserving its sparsity pattern.

**Why consistent mass destroys sparsity:** When M is tridiagonal, L has a band structure but L⁻¹ is dense. The product L⁻¹KL⁻ᵀ fills in.

### 4.3 Why This Matters for Quantum Computing

The number of nonzero Pauli terms in the decomposition of A directly determines the measurement cost. Each Pauli term (or group of commuting terms) requires separate quantum circuit executions.

**The fundamental trade-off:**

The paper frames this as a "three-way" problem (Eq. trilemma), but it's really a direct two-way trade-off:

> **A more physically accurate mass matrix costs more to measure on a quantum computer.**

Sparsity of A is not an independent design variable — it's the deterministic *consequence* of your choice of M:

```
Choice of M  →  determines L  →  determines A = L⁻¹KL⁻ᵀ  →  determines Pauli count  →  determines measurement cost
```

- Lumped mass (diagonal M): Cheap to measure (576 Pauli terms at N=128), but physically less accurate for high-frequency modes
- Consistent mass (tridiagonal M): Better physics, but 5× more Pauli terms (2,950 at N=128)

The paper's data at multiple sizes:

| N | Lumped Pauli terms | Consistent Pauli terms | Ratio |
|---|---|---|---|
| 4 | 8 | 10 | 1.25× |
| 16 | 48 | 136 | 2.8× |
| 64 | 256 | 1,241 | 4.8× |
| 128 | 576 | 2,950 | 5.1× |

The consistent-mass Pauli count grows roughly as O(N²). This quadratic measurement cost destroys the logarithmic qubit advantage.

### 4.4 The Net Effect

```
Storage:     log₂(N) qubits     — exponential savings ✓
Measurement: O(N²) Pauli terms  — quadratic cost ✗
Classical:   O(N·k) operations  — linear cost ✓✓
```

The quantum approach saves exponentially on storage but pays quadratically on measurement. Classical Lanczos on a sparse system costs O(N·k). The quantum approach is asymptotically worse.

---

## Part 5: Circuit Depth and Expressibility

### 5.1 The Representation Capacity Problem

A quantum circuit with n qubits and depth L has n·(L+1) adjustable parameters.

The target (an arbitrary unit vector in 2ⁿ dimensions) has 2ⁿ - 1 degrees of freedom (it's a unit vector on a (2ⁿ - 1)-sphere).

**Necessary condition for the circuit to be expressive enough:**

```
n·(L+1) ≥ 2ⁿ - 1
```

Solving for L:

```
L ≥ (2ⁿ - 1)/n - 1
```

**Worked examples:**
- n=3 (8 dimensions): L ≥ (7/3) - 1 = 1.33 → need at least depth 2
- n=4 (16 dimensions): L ≥ (15/4) - 1 = 2.75 → need at least depth 3
- n=7 (128 dimensions): L ≥ (127/7) - 1 = 17.1 → need at least depth 18

**Key insight:** Required depth grows EXPONENTIALLY with qubit count (since 2ⁿ/n grows exponentially). Having log₂(N) qubits is not enough — you also need depth proportional to N/log(N).

### 5.2 Experimental Confirmation

The paper's exact-state VQE results:

| DOF (N) | Qubits (n) | Depth (L) | Parameters | Frequency error |
|---------|-----------|-----------|------------|-----------------|
| 8 | 3 | 1 | 6 | 27.3% |
| 8 | 3 | 2 | 9 | ~0% (10⁻¹³) |
| 16 | 4 | 1 | 8 | 81.8% |
| 16 | 4 | 2 | 12 | 5.2% |
| 16 | 4 | 3 | 16 | ~0% (10⁻⁹) |

**Interpretation:**
- 8 DOFs need 7 degrees of freedom → 9 parameters (depth 2) is enough ✓
- 16 DOFs need 15 degrees of freedom → 16 parameters (depth 3) is enough ✓
- 12 parameters (depth 2) for 16 DOFs is NOT enough → 5% error remains

The threshold is approximately: parameters ≥ degrees of freedom.

### 5.3 This Is a Capacity Result, Not an Optimization Result

The paper ran 20 random restarts at depth 1 for the 6-DOF problem. Every single one gave ~128% error (minimum: 127.99%, median: 130.20%, maximum: 130.20%).

**What this proves:** The failure is NOT due to the optimizer getting stuck in a local minimum. If it were, different random starts would give different errors. The failure is SYSTEMATIC — the circuit literally cannot represent the target vector at depth 1, regardless of how you optimize.

**CS analogy:** This is like proving that a 2-layer neural network with ReLU activations cannot represent a certain function, regardless of the weights. It's a representational impossibility, not a training difficulty.

### 5.4 Barren Plateaus (Vanishing Gradients)

As qubit count grows, the gradient variance of a randomly initialized circuit decays:

| Qubits | Gradient variance (uniform init) | Gradient variance (near-identity init) |
|--------|------|------|
| 2 | 2.62×10⁻² | 4.05×10⁻⁴ |
| 3 | 1.09×10⁻² | 8.87×10⁻⁴ |
| 4 | 4.51×10⁻³ | 1.36×10⁻⁴ |
| 5 | 3.47×10⁻³ | 3.26×10⁻⁵ |
| 6 | 2.22×10⁻³ | 1.53×10⁻⁵ |
| 7 | 1.59×10⁻³ | 3.14×10⁻⁵ |
| 8 | 2.04×10⁻³ | 2.23×10⁻⁵ |
| 9 | 9.82×10⁻⁴ | 1.97×10⁻⁵ |
| 10 | 1.05×10⁻³ | 8.76×10⁻⁶ |
| 11 | 1.12×10⁻³ | 2.22×10⁻⁵ |
| 12 | 9.82×10⁻⁴ | 1.52×10⁻⁵ |

The gradient shrinks for the first few qubits, then **flattens near 10⁻³** rather than continuing to fall. A barren plateau is *defined* by gradient variance decaying **exponentially** in qubit count; the test is whether `log V` falls linearly in `n_q` (exponential) or in `log n_q` (polynomial).

**What the extended study (2–12 qubits) finds:** a **power law fits far better than an exponential** — full-range R²=0.96 vs 0.83 (favouring the power law by ΔAIC≈15); over the small-system-free range (n_q≥4) the advantage narrows to R²=0.90 vs 0.85 (ΔAIC≈4, moderate). A linear fit to the tail (n_q≥7) gives a shallow slope of −0.12/qubit (95% CI [−0.23, −0.01]) — statistically nonzero but far below the exponential rate that defines a barren plateau. So the gradients keep shrinking, but at an **approximately polynomial, not exponential, rate.**

**CS analogy:** akin to the vanishing-gradient problem in deep nets before ResNets/batchnorm — but here, over the tested range, the decay looks polynomial rather than the exponential collapse a true barren plateau would show.

**Honest caveat from the paper:** n_q ≤ 12 is still far from asymptotic, so this does **not** rule out a barren plateau emerging at larger scale — it only shows no exponential signature over the accessible range.

---

## Part 6: Measurement and Sample Complexity

### 6.1 How You Measure on a Quantum Computer

You cannot directly read out a quantum state. Instead, you decompose the matrix H into a sum of "Pauli strings":

```
H = Σₗ hₗ Pₗ
```

where each Pₗ is a tensor product of 2×2 Pauli matrices:

```
I = [1 0; 0 1]    X = [0 1; 1 0]    Y = [0 -i; i 0]    Z = [1 0; 0 -1]
```

For n qubits, each Pₗ is a product like X⊗Z⊗I⊗Y (one Pauli per qubit). There are at most 4ⁿ possible strings, but most have coefficient hₗ ≈ 0.

**To estimate the energy E = ⟨ψ|H|ψ⟩:**
1. For each Pauli string Pₗ with nonzero hₗ:
   - Prepare the state |ψ⟩
   - Measure in the basis defined by Pₗ
   - Get outcome +1 or -1
   - Repeat many times to estimate the mean μₗ = ⟨ψ|Pₗ|ψ⟩
2. Compute E = Σₗ hₗ · μₗ

### 6.2 Shot Noise (Sample Complexity)

Each measurement of Pₗ gives a ±1 outcome (like a biased coin flip). With n shots:

```
Var(μ̂ₗ) = (1 - μₗ²) / n
```

This is just the variance of a sample mean from a Bernoulli-like distribution. The energy estimate variance is:

```
Var(Ê) = Σₗ hₗ² · (1 - μₗ²) / nₗ
```

**Key point:** This is standard Monte Carlo estimation. The error decreases as 1/√n. There is NO quantum speedup in the measurement phase.

### 6.3 Why the Lowest Frequency Is Hardest

The physical frequency is f = √(s·ε) / (2π). By the delta method (first-order Taylor expansion of a function of a random variable):

```
Var(f̂) ≈ (s / (8π²f))² · Var(ε̂)
```

The factor 1/f² means: **the smallest frequency has the largest relative error.**

**Intuition:** The square root function has a steep slope near zero. A small absolute error in ε (near zero) produces a large relative error in f = √ε.

**Numerical confirmation from the paper (at 10⁵ shots/term):**

| Mode | Frequency (Hz) | RMSE (%) |
|------|---------------|----------|
| 1 | 4.48 | 2.23 |
| 2 | 12.27 | 0.34 |
| 3 | 19.38 | 0.15 |
| 4 | 25.18 | 0.08 |

Mode 1 is **28× harder** to estimate than mode 4 at the same shot budget.

### 6.4 Can Grouped Measurements Help?

Yes, but not enough. "Qubit-wise commuting" (QWC) groups are sets of Pauli strings that can be measured simultaneously (in one circuit execution).

The paper computes an optimistic lower bound: if you perfectly group measurements and optimally allocate shots across groups, the best possible RMSE for mode 1 at the same total budget (1.9×10⁶ total shots) is **0.83%**.

This is better than the 2.23% from independent measurement, but still substantial. The grouping helps by a factor of ~2.7, but doesn't fundamentally solve the problem.

### 6.5 The Damage Detection Failure

**The downstream task:** Given modal data from a "healthy" structure and a "possibly damaged" structure, identify which element is damaged.

**Method:** Compute the "modal strain energy fraction" for each element e and each mode i:

```
η_{e,i} = φᵢᵀKₑφᵢ / φᵢᵀKφᵢ
```

The damage index for element e is:

```
D_e = Σᵢ wᵢ |η_{e,i}^damaged - η_{e,i}^baseline|
```

Declare the element with the largest D_e as damaged.

**Results at 10⁵ shots/term with 10% stiffness loss:**
- Correct localization rate: **26%** (random guessing = 16.7%)
- Median margin: **negative** (wrong element usually scores higher)

**How many shots do you actually need?**

The paper derives (using the delta method and hypothesis testing theory):

| Stiffness loss | Shots/term for detection (80% power) | Shots/term for localization (80% power) |
|---|---|---|
| 2% | 1.32×10⁷ | 1.27×10⁹ |
| 5% | 2.04×10⁶ | 2.11×10⁸ |
| 10% | 4.82×10⁵ | 5.71×10⁷ |
| 20% | 1.10×10⁵ | 1.83×10⁷ |

For 10% damage with 648 measurement settings: **5.71×10⁷ × 648 ≈ 3.7×10¹⁰ total circuit shots**.

**CS interpretation:** This is like saying a classifier needs 37 billion training samples to achieve 80% accuracy on a binary decision. The sample complexity makes the approach impractical.

---

## Part 7: The Normalization and Padding Scheme

### 7.1 Why Normalize?

The eigenvalues of A can be very large (e.g., 43,000 in the test cases). Quantum optimization works best when eigenvalues are in [0, 1]. So divide by a scale factor s_t.

### 7.2 The Gershgorin Bound

**Theorem (Gershgorin, 1931):** For any symmetric matrix A, every eigenvalue lies within at least one "Gershgorin disc." The largest eigenvalue satisfies:

```
λ_max ≤ U_G = max_i (A_ii + Σ_{j≠i} |A_ij|)
```

**Proof sketch:** If Av = λv, look at the row i where |v_i| is largest. Then:
```
λ v_i = A_ii v_i + Σ_{j≠i} A_ij v_j
|λ - A_ii| ≤ Σ_{j≠i} |A_ij| · |v_j/v_i| ≤ Σ_{j≠i} |A_ij|
```
So λ ≤ A_ii + Σ_{j≠i} |A_ij| for some i. The max over i gives the bound.

**Cost:** O(N²) — just sum each row. No eigensolve needed.

**How tight is it?** From the paper's experiments, the ratio U_G/λ_max is between 1.014 and 1.110 for these 1D chains. So the bound is only 1-11% loose.

### 7.3 Why Pad?

If N=6, you need ⌈log₂(6)⌉ = 3 qubits, giving a state space of dimension 2³ = 8. You must embed the 6×6 matrix into an 8×8 matrix.

**The padding scheme:**

```
H = [A/s_t    0    ]     (6×6 block in top-left)
    [  0    α·I₂   ]     (2×2 identity scaled by α in bottom-right)
```

where α = U_G/s_t + δ (strictly above the normalized physical bound).

**Safety guarantee:** The padding eigenvalues equal α, which is strictly above U_G/s_t ≥ λ_max/s_t. So padding states always have higher energy than any physical state. The optimizer (which minimizes energy) can never accidentally converge to a padding state.

**CS analogy:** This is like adding sentinel values to a sorted array that are guaranteed larger than all real values. A min-search will never return a sentinel.

### 7.4 The Two-Scale Idea

Previous work used a single scale s for both normalization and padding safety. This paper separates them:

1. **Optimization scale s_t:** Makes gradients well-conditioned. Can be estimated cheaply (a few Lanczos iterations). Doesn't need to be ≥ λ_max.

2. **Certified bound U_G:** Guarantees padding safety. Must be ≥ λ_max. Doesn't need to be tight.

**Why separate?** If you use U_G (which may be conservative) as the optimization scale, you compress the physical spectrum near zero, making gradients tiny. If you use s_t (which may be < λ_max) for padding, you might not safely exclude padding states. The two-scale scheme avoids both problems.

---

## Part 8: Validation and Diagnostics

### 8.1 The Diagnostic Vector (Why No Single Score)

The paper defines a per-mode diagnostic vector:

```
D_i = (ε_{f,i}, ρ_i, σ²_{H,i}, ℓ_{p,i}, o_i, 1 - MAC_ii)
```

| Component | Formula | What it measures | Failure mode it catches |
|-----------|---------|-----------------|------------------------|
| ε_{f,i} | \|f_quantum - f_exact\| / f_exact | Frequency accuracy | Wrong eigenvalue |
| ρ_i | ‖Kφ - λMφ‖ / (‖Kφ‖ + λ‖Mφ‖) | Residual of the eigenequation | Not actually an eigenvector |
| σ²_{H,i} | ⟨H²⟩ - ⟨H⟩² | Hamiltonian variance | Superposition of eigenstates, not a pure eigenstate |
| ℓ_{p,i} | Probability in padding dimensions | Padding leakage | Answer leaked into fake dimensions |
| o_i | max overlap with prior modes | Mode contamination | Deflation didn't fully separate modes |
| 1 - MAC_ii | 1 - (cosine similarity)² with classical mode | Shape mismatch | Right eigenvalue but wrong eigenvector |

**Why not combine into one number?** A good value in one component can hide catastrophic failure in another. Example: a state with 99% overlap with mode 1 and 1% leakage into padding might have a reasonable energy (low ε_f) but terrible shape accuracy (low MAC). A single score would average these out.

**CS analogy:** This is like reporting precision, recall, F1, AUC, and calibration separately for a classifier. "95% accuracy" can hide complete failure on the minority class.

### 8.2 Modal Assurance Criterion (MAC)

```
MAC_ij = |φᵢᵀMφⱼ|² / [(φᵢᵀMφᵢ)(φⱼᵀMφⱼ)]
```

This is cosine similarity squared in the M-weighted inner product:
- MAC = 1 → vectors are parallel (perfect match)
- MAC = 0 → vectors are orthogonal (completely different)

For mass-normalized modes (φᵀMφ = 1), this simplifies to MAC_ij = |φᵢᵀMφⱼ|².

### 8.3 Mode Tracking Over Time

**Problem:** In structural health monitoring, you solve the eigenvalue problem every day. Today's mode 1 might not correspond to yesterday's mode 1 (modes can swap order if frequencies are close).

**Solution (Hungarian algorithm):** 
1. Solve each day independently (no temporal penalty in the optimizer)
2. Build a cost matrix: c_ij = how different is today's mode j from yesterday's mode i?
3. Find the minimum-cost assignment (bipartite matching)

**Why NOT add a tracking penalty to the optimizer:** If you penalize the optimizer for deviating from yesterday's answer, you bias today's solution toward yesterday's. This masks damage (a damaged mode would be pulled toward its healthy predecessor).

**CS analogy:** This is "tracking by detection" in multi-object tracking (computer vision). You run a detector independently on each frame, then match detections across frames. You never let the tracker influence the detector.

---

## Part 9: The Complete Resource Accounting

### 9.1 The Resource Vector

The paper reports (and refuses to combine into a single number):

```
R = (n_q, N_P, N_G, N_shot, N_eval, D_circ, N_2q,
     C_asm, C_chol, C_bound, C_Pauli, C_prep)
```

| Symbol | Meaning | Example value (6-DOF) |
|--------|---------|----------------------|
| n_q | Qubits | 3 |
| N_P | Nonzero Pauli terms | 20 |
| N_G | Measurement groups (QWC) | 8 |
| N_shot | Total shots per mode | 10⁵ × 19 = 1.9×10⁶ |
| N_eval | Optimizer function evaluations | ~450 per restart × 8 restarts |
| D_circ | Circuit depth | 2 (for exact recovery) |
| N_2q | Two-qubit gates per circuit | 6 (CNOT ring × depth) |
| C_asm | Classical assembly cost | O(N) for sparse structures |
| C_chol | Cholesky factorization cost | O(N) banded (1D chains); O(N^{3/2}) 2D; O(N³) dense |
| C_bound | Spectral bound computation | O(N²) for Gershgorin |
| C_Pauli | Pauli decomposition cost | O(4ⁿ · N) naively |
| C_prep | State preparation cost | O(N) gates worst case |

### 9.2 Why No Single "Efficiency" Metric

Each resource has different scaling, different hardware constraints, and different costs:
- Qubits are limited by hardware (currently ~1000 noisy qubits)
- Circuit depth is limited by decoherence time
- Shots are limited by wall-clock time
- Classical preprocessing is limited by CPU/memory

Combining them requires arbitrary weighting (how many shots equal one qubit?). The paper argues this is like combining FLOPs, memory, latency, and energy into one "performance" number — it's meaningless without a specific cost model.

---

## Part 10: Experimental Design Philosophy

### 10.1 Separating Error Sources

The paper runs five experiments, each isolating a different failure mode:

| Experiment | What it uses | What it isolates |
|-----------|-------------|-----------------|
| A: Deflation threshold | Exact eigenvalues, exact states | Is the theorem correct? |
| B: Hamiltonian scaling | Exact Cholesky, count Pauli terms | How does preprocessing cost grow? |
| C: Depth/restarts | Exact state vectors (no noise) | Is failure from capacity or optimization? |
| D: Finite-shot frequency | Exact states, sampled measurements | How does shot noise affect output? |
| E: Finite-shot damage | Exact states, sampled measurements | Can you make a downstream decision? |

**The principle:** Don't run everything on noisy hardware and get one combined error number. Instead:
- Tests A-C use exact state vectors → isolates ansatz/optimizer error
- Tests D-E use exact states but sampled measurements → isolates shot noise
- No hardware noise anywhere → that's explicitly declared as a separate study

**CS analogy:** This is like the ML practice of:
1. First, test your model on clean data (check architecture is correct)
2. Then, add label noise (check robustness to noise)
3. Then, add distribution shift (check generalization)

Rather than training on noisy, shifted data from the start and not knowing which errors come from where.

---

## Part 11: What the Paper Proves vs. What It Doesn't

### 11.1 Proven (Exact Guarantees)

1. **Cholesky whitening preserves eigenvalues** and converts M-orthogonality to Euclidean orthogonality. (Proposition 1, proved by direct calculation.)

2. **Safe padding excludes artificial states** from the physical low spectrum for any positive optimization scale. (Proposition 2, proved by construction: padding eigenvalue = α > U_G/s_t ≥ λ_max/s_t.)

3. **The deflation threshold is necessary and sufficient.** (Theorem 1, proved by explicit eigenvalue calculation of the deflated operator.)

4. **Approximate projectors degrade the solution proportionally to β × sin(angle error).** (Derived from Davis-Kahan theorem, a standard result in matrix perturbation theory.)

### 11.2 Not Proven (Open Questions)

1. **Global convergence of the optimizer** — The cost function is nonconvex. No guarantee that gradient descent finds the global minimum. (This is the same as the open question of whether SGD finds global minima in deep learning.)

2. **Barren plateaus** — Still cannot prove/disprove an *asymptotic* barren plateau, but the extended 2–12 qubit study finds the gradient variance decaying *polynomially*, not exponentially (power law beats exponential, ΔAIC≈15 full-range), i.e. no barren-plateau signature over the accessible range; 12 qubits is still not asymptotic.

3. **Any quantum advantage** — The paper explicitly states: "These results do not establish quantum advantage."

4. **Pauli scaling law** — The N² reference line is empirical. The paper does NOT prove whether the Pauli term count is O(N²), O(N log N), or something else.

5. **Formal comparison to Lanczos** — The paper says classical is "overwhelmingly preferable" but never writes a formal T_quantum vs T_classical expression.

### 11.3 The Paper's Actual Conclusion

> "For the present benchmarks, classical methods are overwhelmingly preferable."

And:

> "These results do not establish quantum advantage. They identify the mathematical conditions and resource bottlenecks that a future structural quantum algorithm must overcome."

This is deliberately weaker than an impossibility result. It's saying: "We can't prove it's impossible, but we've shown exactly how expensive it is, and right now it's far worse than classical."

---

## Part 12: The Fundamental Bottlenecks (Summary)

### 12.1 Where the Costs Come From

```
TOTAL QUANTUM COST ≈ 
    Classical preprocessing (banded Cholesky O(N) for 1D chains)
  + Whitening + Pauli decomposition (densifies A → O(N²) terms)
  + Per mode: shots × Pauli_terms × optimizer_iterations
  + Mode recovery (O(N²) for back-transformation)
```

### 12.2 Comparison to Classical

```
CLASSICAL LANCZOS COST ≈ O(k · nnz(K) · iterations)
```

For a sparse 1D chain with nnz = O(N), finding k modes costs O(N·k·iter).

The quantum approach (corrected — preprocessing is *not* the dominant term for 1D chains):
- The banded Cholesky is only O(N); it is the *whitening + Pauli construction* that hurts — it densifies A and produces O(N²) Pauli terms. Merely building and storing that operator is Θ(N²), already more than the entire classical solve for large N, before a single shot.
- Then measurement: O(N²) Pauli terms × O(1/ε²) shots per term = O(N²/ε²) per mode
- Classical: O(N·k) total

### 12.3 The Five Bottlenecks, Ranked

1. **Whitening-induced densification (operator construction):** the banded Cholesky itself is only O(N), but whitening densifies the operator and building its Pauli representation is Θ(N²) in both time and space — already more than the entire classical solve for large N, before any measurement. (An earlier draft mislabeled this as O(N³) Cholesky; that dense cost does not apply to the 1D banded chains the paper studies.)

2. **Pauli term count (measurement cost):** Grows as ~O(N²) for consistent mass. Each term needs separate measurements. Destroys the log(N) qubit advantage.

3. **Circuit depth:** Must grow with N (not just log N). Required depth ≈ N/log(N). Each layer adds noise on real hardware.

4. **Shot noise for low-frequency modes:** The 1/f² amplification means the most important mode (the fundamental frequency) is the hardest to estimate. No quantum speedup in measurement.

5. **Downstream decision (damage detection):** Even with perfect eigenvalues, detecting small damage requires ~10¹⁰ total circuit shots. The signal-to-noise ratio is too low.

---

## Part 13: What Would Need to Change for Quantum to Win

The paper identifies specific conditions that a future quantum algorithm must satisfy:

1. **Avoid dense whitening:** Use a generalized Rayleigh quotient (measure both K and M separately) or domain decomposition to preserve sparsity. But this introduces new problems (denominator estimation, M-metric orthogonality enforcement).

2. **Reduce Pauli terms:** Use alternative encodings, exploit structure in the matrix, or use classical shadows (a measurement technique that estimates many observables from fewer circuit runs).

3. **Solve the depth problem:** Find ansatz architectures that are expressive enough with depth O(poly(log N)) rather than O(N/log N). This is an open research question.

4. **Reduce shot requirements:** Use correlated measurement strategies, variance reduction techniques, or reformulate the damage detection problem to require less precision.

5. **Find problems where N is enormous but classical sparsity doesn't help:** If nnz(K) = O(N²) (dense classical problem), then classical Lanczos costs O(N²·k) and the quantum approach might become competitive. But structural matrices are typically sparse.

---

## Part 14: Glossary of Domain-Specific Terms

| Term | Domain | Meaning |
|------|--------|---------|
| DOF (degree of freedom) | Engineering | One independent displacement coordinate (e.g., horizontal motion of one floor) |
| Mode shape | Engineering | Eigenvector — the pattern of deformation at a natural frequency |
| Natural frequency | Engineering | √eigenvalue / (2π) — a frequency at which the structure resonates |
| Stiffness matrix K | Engineering | Encodes how resistant the structure is to deformation |
| Mass matrix M | Engineering | Encodes how heavy each part of the structure is |
| Lumped mass | Engineering | Diagonal M — mass concentrated at nodes (simpler, less accurate) |
| Consistent mass | Engineering | Non-diagonal M — mass distributed along elements (more accurate) |
| Qubit | Quantum | A quantum bit — a unit vector in 2D complex space |
| Ansatz | Quantum | The parameterized circuit family (like a neural network architecture) |
| Hamiltonian | Quantum/Physics | The matrix whose eigenvalues you want to find |
| Pauli string | Quantum | A tensor product of {I,X,Y,Z} matrices — the measurement basis |
| Shot | Quantum | One execution of the quantum circuit + measurement |
| Barren plateau | Quantum | Exponentially vanishing gradients in quantum circuits |
| VQE | Quantum | Variational Quantum Eigensolver — finds minimum eigenvalue |
| VQD | Quantum | Variational Quantum Deflation — finds excited states sequentially |
| MAC | Engineering | Modal Assurance Criterion — cosine similarity² of mode shapes |
| Gershgorin disc | Linear algebra | A circle in the complex plane guaranteed to contain an eigenvalue |
| Cholesky factorization | Linear algebra | Decomposing a positive-definite matrix as M = LLᵀ |
| Davis-Kahan theorem | Linear algebra | Bounds eigenvector perturbation given eigenvalue gap and matrix perturbation |
| Delta method | Statistics | First-order variance propagation through a nonlinear function |
| Hungarian algorithm | Combinatorics | Optimal bipartite matching in O(n³) |

---

## Part 15: Reading the Code

### 15.1 `reproduce_study.py` — The Core Proof of Concept

This script implements the entire pipeline for a 6-storey building:

```python
# Build the structural model (6 floors, each with mass and stiffness)
model = build_shear_building()  # Returns M, K, element stiffnesses

# Classical preprocessing
A, L = mass_normalized_operator(M, K)  # Cholesky whitening
scale = gershgorin_upper(A)            # Compute U_G
H = safe_pad(A/scale, target_dim=8)    # Pad 6×6 → 8×8

# Quantum simulation (no actual quantum hardware)
result = solve_vqd(H, n_modes=4, depth=2, ...)  # VQD with L-BFGS-B optimizer

# Recovery
Phi = recover_modes(result.states, L, physical_dim=6, M=model.M)
```

The `RealAmplitudeCircuit` class simulates a quantum circuit classically:
- Ry gates (single-qubit rotations) parameterized by angles θ
- CNOT gates (two-qubit entangling gates) in a ring topology
- The state is computed as a 2ⁿ-dimensional vector (exact, no noise)

### 15.2 `revision_study.py` — The Scaling Studies

Extends the analysis to:
- Chains with N = 4, 8, 16, 32, 64, 128 DOFs
- Both lumped and consistent mass
- Pauli decomposition with optimized algorithm (avoids building 2ⁿ×2ⁿ Kronecker products)
- Finite-shot sampling (simulates quantum measurement noise)
- Gradient variance measurement (barren plateau diagnostic)

### 15.3 `minor_revision_study.py` — Measurement Floor and Shot Requirements

Computes:
- The theoretical best-case RMSE with grouped (QWC) measurements
- The delta-method shot requirements for damage detection and localization

### 15.4 `validate_outputs.py` — Integrity Checks

Verifies:
- All `\cite{}` keys in the LaTeX have corresponding entries in `references.bib`
- Specific numerical values match expected results (e.g., exactly 2950 Pauli terms for consistent mass at N=128)

---

## Part 16: Open Problems (Potential CS Contributions)

### 16.1 Prove the Pauli Scaling Law

**Question:** For a 1D chain with consistent mass and Cholesky whitening, is the number of nonzero Pauli terms Θ(N²)?

**What's known:** The paper shows 10 → 136 → 1241 → 2950 for N = 4, 16, 64, 128. An independent recompute (`scratchpad/verify_complexity.py`) fits a *local* log-log slope of **≈1.4–1.5** over this range — i.e. **sub-quadratic**, so O(N²) is a worst-case *upper bound* (4^{n_q} ≤ 4N²), not a tight law. Whether the true asymptotic exponent is 2 (constant limiting density) or 2−α (density decaying as 1/N^α) is still unproven.

**Approach:** The Pauli term count equals the number of nonzero entries in the matrix A when expressed in the computational basis (roughly). Since consistent-mass whitening fills A to ~30% density at N=128, and density × N² = nonzeros, the Pauli count should be O(density × N²). If density converges to a constant, it's Θ(N²). If density decays as 1/N^α, it's O(N^{2-α}).

### 16.2 Write the End-to-End Complexity

**Status: now done in the paper** — this is `tab:complexity` (Table 1) in Section VI, with a worked 128-DOF example (`tab:concrete`). Corrected sketch (banded 1D chain):
```
T_quantum(N, k, ε) = O(N)                     [banded Cholesky, 1D chains]
                   + O(N² log N)              [whitening + Pauli decomposition]
                   + k · N_P · shots · iters   [VQE per mode, N_P = O(N²)]
                   + k · O(N²)                 [recovery]
```

**Compare to:**
```
T_classical(N, k, ε) = O(k · nnz(K) · log(1/ε))  [Lanczos]
```

For sparse K with nnz = O(N): T_classical = O(N·k·log(1/ε)).

### 16.3 Prove a Conditional Impossibility

**Candidate theorem:** "Under binary amplitude encoding with Cholesky whitening of a consistent-mass 1D chain, the total measurement cost for k modes to accuracy ε is Ω(N²·k/ε²), which exceeds the classical cost O(N·k·log(1/ε)) for all N where nnz(K) = O(N)."

This would be a formal separation result (conditional on the encoding choice).

### 16.4 Find the Crossover Point

**Question:** Is there any N (and any problem structure) where the quantum approach becomes competitive?

**Candidate:** Dense classical problems where nnz(K) = O(N²). Then classical Lanczos costs O(N²·k), and if the quantum Pauli count is also O(N²) but with a smaller constant (due to grouping), there might be a crossover. But this requires the preprocessing to not dominate.

---

## Part 17: One-Paragraph Summary

A quantum computer can store a structural eigenvector in log(N) qubits, but actually computing that eigenvector requires: (1) classical preprocessing that is itself cheap for a banded 1D chain (O(N) Cholesky) but destroys the sparsity via whitening — densifying the operator so that merely building its O(N²) Pauli representation can cost more than solving the problem classically; (2) circuit depth that must grow proportionally to N/log(N), not just log(N); (3) a deflation penalty for finding multiple eigenvalues that must be set precisely — the paper proves the exact necessary-and-sufficient threshold and shows that overshooting amplifies errors from imperfect prior modes; (4) O(N²) Pauli measurements in the worst case for consistent mass; and (5) approximately 3.7×10¹⁰ total circuit shots to make a practical engineering decision (damage localization) at 10% damage severity. Classical Lanczos solves the same problem in O(N·k) time for sparse systems. The quantum approach is currently far worse, and this paper identifies exactly where each bottleneck lies, with theorems and controlled experiments, rather than making vague claims about potential future advantage.
