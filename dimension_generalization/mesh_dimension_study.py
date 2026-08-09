"""Dimension-generalization study: 2D and 3D scalar mesh extensions of the
paper's 1D shear-chain model.

Purpose: the manuscript's Cholesky-fill and Pauli-density results are
demonstrated only on 1D chains, the *best case* for sparsity. This script
tests whether the same effect (whitening destroys sparsity; nonzero Pauli
terms grow with N) holds -- and worsens -- as spatial dimension increases,
and validates the manuscript's subspace/MAC degeneracy diagnostics (Eq. mac,
Eq. subspace) on a genuinely degenerate case, which no 1D chain can produce.

Design choices:
  - The mesh is a SCALAR grid (one DOF per node, nearest-neighbor springs),
    the direct generalization of the 1D chain -- not the VECTOR truss/beam
    elements already sitting unused in Appendix A. This isolates the effect
    of spatial dimension from the confound of DOFs-per-node.
  - All dimension-agnostic linear algebra (Cholesky whitening, Pauli
    decomposition, QWC grouping, Gershgorin bound) is imported directly from
    revision_study.py, not reimplemented, so these numbers are produced by
    the same code that generated the paper's 1D Table 2.
  - Two different boundary conditions are used for two different purposes:
      * grid_MK: fixes an entire "wall" layer along axis 0 (generalizing the
        1D chain's fixed base), chosen so free-DOF counts land on exact
        powers of 2. Used for the scaling sweep. Properties are smoothly
        heterogeneous (asymmetric in the grid axes) to avoid accidental
        symmetry degeneracies inflating or deflating the sparsity count.
      * square_uniform_MK: pins only the (0,0) corner (a single Dirichlet
        constraint that preserves the i<->j reflection symmetry of a square
        grid) with UNIFORM properties, deliberately chosen to produce
        genuine degenerate mode pairs. Used only for the degeneracy demo.

Written in a separate subdirectory to keep the project root uncluttered.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parent
PARENT = ROOT.parent
DATA = ROOT / "data"
FIG = ROOT / "figures"
DATA.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

sys.path.insert(0, str(PARENT))
from revision_study import mass_whiten, pauli_coefficients, qwc_groups, gersh  # noqa: E402

K0 = 2.2e8
M0 = 2.0e4


def _local_blocks(k_val: float, m_val: float, mass_type: str):
    ke2 = k_val * np.array([[1.0, -1.0], [-1.0, 1.0]])
    if mass_type == "consistent":
        me2 = m_val / 6.0 * np.array([[2.0, 1.0], [1.0, 2.0]])
    elif mass_type == "lumped":
        me2 = m_val / 2.0 * np.eye(2)
    else:
        raise ValueError(mass_type)
    return ke2, me2


def _assemble_grid_KM(shape: tuple[int, ...], mass_type: str, heterogeneous: bool
                       ) -> tuple[np.ndarray, np.ndarray]:
    """Nearest-neighbor scalar spring/mass network on a d-dimensional grid.
    ``shape`` fixes the FULL node grid (no boundary conditions applied here).
    """
    dims = len(shape)
    n = int(np.prod(shape))
    idx = np.arange(n).reshape(shape)
    K = np.zeros((n, n))
    M = np.zeros((n, n))

    if heterogeneous:
        coords = np.indices(shape, dtype=float)  # (dims, *shape)
        # Asymmetric per-axis weights: deliberately breaks any i<->j grid
        # symmetry so the scaling sweep never hits an accidental degeneracy.
        weights = np.array([2.0 ** d for d in range(dims)])
        denom = sum(w * max(s - 1, 1) for w, s in zip(weights, shape))
        radial = sum(w * coords[d] for d, w in enumerate(weights)) / denom
    else:
        radial = np.zeros(shape)

    for d in range(dims):
        lo = [slice(None)] * dims
        hi = [slice(None)] * dims
        lo[d] = slice(0, shape[d] - 1)
        hi[d] = slice(1, shape[d])
        a_ids = idx[tuple(lo)].ravel()
        b_ids = idx[tuple(hi)].ravel()
        r_a = radial[tuple(lo)].ravel()
        r_b = radial[tuple(hi)].ravel()
        for a, b, ra, rb in zip(a_ids, b_ids, r_a, r_b):
            r = 0.5 * (ra + rb)
            k_val = K0 * (1.0 - 0.35 * r)
            m_val = (M0 / dims) * (1.0 - 0.20 * r)
            ke2, me2 = _local_blocks(k_val, m_val, mass_type)
            ids = (int(a), int(b))
            for p in range(2):
                for q in range(2):
                    K[ids[p], ids[q]] += ke2[p, q]
                    M[ids[p], ids[q]] += me2[p, q]
    return M, K


def grid_MK(free_shape: tuple[int, ...], mass_type: str = "consistent",
            heterogeneous: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """Free-DOF (M, K) for a d-dimensional grid fixed on the axis-0=0 wall
    layer -- the direct generalization of the 1D chain's fixed base.
    ``free_shape`` gives the FREE grid size; total free DOFs = prod(free_shape).
    """
    full_shape = (free_shape[0] + 1,) + tuple(free_shape[1:])
    M_full, K_full = _assemble_grid_KM(full_shape, mass_type, heterogeneous)
    idx_full = np.arange(M_full.shape[0]).reshape(full_shape)
    free_ids = idx_full[1:].ravel()
    return M_full[np.ix_(free_ids, free_ids)], K_full[np.ix_(free_ids, free_ids)]


def square_uniform_MK(n_side: int, mass_type: str = "consistent"
                       ) -> tuple[np.ndarray, np.ndarray]:
    """Uniform square grid, free everywhere except a single pinned corner
    (0,0). Preserves the i<->j reflection symmetry needed for genuine
    degenerate mode pairs; used only for the degeneracy demonstration.
    """
    M_full, K_full = _assemble_grid_KM((n_side, n_side), mass_type, heterogeneous=False)
    free_ids = np.arange(1, M_full.shape[0])  # drop only the (0,0) corner
    return M_full[np.ix_(free_ids, free_ids)], K_full[np.ix_(free_ids, free_ids)]


def _density_nnz(A: np.ndarray) -> tuple[int, float]:
    thr = 1e-12 * np.max(np.abs(A))
    nnz = int(np.count_nonzero(np.abs(A) > thr))
    return nnz, nnz / A.size


def scaling_sweep() -> pd.DataFrame:
    cases = (
        [("2D", (4, 4)), ("2D", (8, 8)), ("2D", (16, 16)), ("2D", (32, 32))]
        + [("3D", (4, 4, 4)), ("3D", (4, 4, 8)), ("3D", (8, 8, 16))]
    )
    rows = []
    for dim_label, shape in cases:
        n = int(np.prod(shape))
        nq = int(round(np.log2(n)))
        assert 2 ** nq == n, f"shape {shape} does not give a power-of-2 N={n}"
        for mass_type in ("lumped", "consistent"):
            M, K = grid_MK(shape, mass_type)
            A, _, _, _ = mass_whiten(M, K)
            lam_max = float(eigh(A, eigvals_only=True, subset_by_index=[n - 1, n - 1])[0])
            sG = gersh(A)
            nnz, density = _density_nnz(A)
            H = A / lam_max
            pdf = pauli_coefficients(H, tol=1e-10)
            groups = qwc_groups(pdf) if n <= 64 else np.nan
            rows.append({
                "dimension": dim_label, "shape": str(shape), "mass_type": mass_type,
                "dof": n, "qubits": nq, "gershgorin_ratio": sG / lam_max,
                "nnz_A": nnz, "density_A": density,
                "pauli_terms": int(len(pdf)), "qwc_groups": groups,
            })
            print(f"[{dim_label} {shape} {mass_type:10s}] N={n:4d} n_q={nq} "
                  f"density={density:.4f} pauli_terms={len(pdf):5d} "
                  f"gersh_ratio={sG/lam_max:.4f} qwc={groups}")
    df = pd.DataFrame(rows)
    df.to_csv(DATA / "dimension_scaling.csv", index=False)
    return df


def _comparison_figure(df_hd: pd.DataFrame) -> None:
    df_1d = pd.read_csv(PARENT / "revision_data" / "hamiltonian_scaling.csv")
    df_1d = df_1d[df_1d.mass_type == "consistent"].copy()
    df_1d["dimension"] = "1D"

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    markers = {"1D": "o", "2D": "s", "3D": "^"}
    for dim_label, marker in markers.items():
        if dim_label == "1D":
            g = df_1d.sort_values("dof")
        else:
            g = df_hd[(df_hd.dimension == dim_label) & (df_hd.mass_type == "consistent")].sort_values("dof")
        axes[0].loglog(g["dof"], g["density_A"], marker=marker, label=f"{dim_label} (consistent mass)")
        axes[1].loglog(g["dof"], g["pauli_terms"], marker=marker, label=f"{dim_label} (consistent mass)")
    axes[0].set_xlabel("Structural DOFs, $N$")
    axes[0].set_ylabel("Density of $\\mathbf{A}$")
    axes[1].set_xlabel("Structural DOFs, $N$")
    axes[1].set_ylabel("Nonzero Pauli terms")
    for ax in axes:
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "dimension_scaling.pdf", bbox_inches="tight")
    fig.savefig(FIG / "dimension_scaling.png", dpi=250, bbox_inches="tight")
    plt.close(fig)


def degeneracy_demo(n_side: int = 4) -> pd.DataFrame:
    """Validate Eq. mac / Eq. subspace on a genuine degenerate pair from a
    symmetric square grid -- the first case in the paper able to produce one.
    """
    M, K = square_uniform_MK(n_side, mass_type="consistent")
    lam, Phi = eigh(K, M)  # generalized eig; Phi columns are M-orthonormal
    gaps = np.diff(lam)
    i = int(np.argmin(np.abs(gaps)))  # tightest adjacent pair
    lam_a, lam_b = lam[i], lam[i + 1]
    rel_gap = abs(lam_a - lam_b) / max(abs(lam_a), abs(lam_b), 1e-30)
    phi_a, phi_b = Phi[:, i], Phi[:, i + 1]
    print(f"[degeneracy] n_side={n_side}, mode pair ({i},{i+1}): "
          f"lambda_a={lam_a:.6e}, lambda_b={lam_b:.6e}, relative gap={rel_gap:.3e}")

    def mac(u: np.ndarray, v: np.ndarray) -> float:
        num = (u @ (M @ v)) ** 2
        den = (u @ (M @ u)) * (v @ (M @ v))
        return float(num / den)

    rows = []
    for theta_deg in range(0, 91, 5):
        th = np.deg2rad(theta_deg)
        phi_a_rot = np.cos(th) * phi_a + np.sin(th) * phi_b
        phi_b_rot = -np.sin(th) * phi_a + np.cos(th) * phi_b
        mac_a = mac(phi_a, phi_a_rot)
        mac_b = mac(phi_b, phi_b_rot)
        Yc = np.column_stack([phi_a, phi_b])
        Yq = np.column_stack([phi_a_rot, phi_b_rot])
        gram = Yc.T @ (M @ Yq)
        sv = np.linalg.svd(gram, compute_uv=False)
        sv = np.clip(sv, -1.0, 1.0)
        eps_sub = float(np.sqrt(np.sum(1.0 - sv ** 2)))
        rows.append({"theta_deg": theta_deg, "mac_a": mac_a, "mac_b": mac_b, "eps_sub": eps_sub})
    df = pd.DataFrame(rows)
    df.to_csv(DATA / "degeneracy_demo.csv", index=False)
    print(df.to_string(index=False))

    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax.plot(df.theta_deg, df.mac_a, "o-", label=r"MAC$(\phi_a,\phi_a')$ (naive, per-mode)")
    ax.plot(df.theta_deg, df.eps_sub, "s-", label=r"$\varepsilon_{\mathrm{sub}}$ (subspace angle)")
    ax.set_xlabel(r"Rotation within degenerate subspace, $\theta$ (deg)")
    ax.set_ylabel("Diagnostic value")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "degeneracy_demo.pdf", bbox_inches="tight")
    fig.savefig(FIG / "degeneracy_demo.png", dpi=250, bbox_inches="tight")
    plt.close(fig)
    return df


def main() -> None:
    df_scale = scaling_sweep()
    _comparison_figure(df_scale)
    degeneracy_demo()


if __name__ == "__main__":
    main()
