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

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt

# IEEE prohibits Type 3 fonts in submitted PDFs. Matplotlib's default pdf.fonttype is 3;
# 42 emits TrueType instead. Must be set before any figure is created.
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
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


def _assemble_grid_KM(shape: tuple[int, ...], mass_type: str, heterogeneous: bool,
                      axis_scale: tuple[float, ...] | None = None,
                      probe: tuple[int, float] | None = None
                      ) -> tuple[np.ndarray, np.ndarray]:
    """Nearest-neighbor scalar spring/mass network on a d-dimensional grid.
    ``shape`` fixes the FULL node grid (no boundary conditions applied here).

    ``axis_scale[d]`` multiplies the stiffness of every axis-``d`` spring; used by
    near_degenerate_sweep() to break the square grid's i<->j reflection symmetry by
    a controlled amount. ``probe=(edge_index, fraction)`` reduces one spring's
    stiffness by ``fraction``, the fixed model discrepancy whose effect on the modes
    is measured. Both default to None, reproducing the original assembly exactly.
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

    edge = 0
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
            if axis_scale is not None:
                k_val *= axis_scale[d]
            if probe is not None and edge == probe[0]:
                k_val *= (1.0 - probe[1])
            edge += 1
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


def square_uniform_MK(n_side: int, mass_type: str = "consistent",
                      axis_scale: tuple[float, ...] | None = None,
                      probe: tuple[int, float] | None = None
                      ) -> tuple[np.ndarray, np.ndarray]:
    """Uniform square grid, free everywhere except a single pinned corner
    (0,0). Preserves the i<->j reflection symmetry needed for genuine
    degenerate mode pairs; used only for the degeneracy demonstration.
    """
    M_full, K_full = _assemble_grid_KM((n_side, n_side), mass_type, heterogeneous=False,
                                       axis_scale=axis_scale, probe=probe)
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


ORDERING_CASES = (("2D", (16, 16)), ("3D", (8, 8, 8)))   # N = 256 and N = 512
ORDERING_SEED = 20260623


def _pattern(A: np.ndarray, rel: float = 1e-12) -> np.ndarray:
    return np.abs(A) > rel * np.max(np.abs(A))


def _bandwidth(A: np.ndarray) -> int:
    i, j = np.nonzero(_pattern(A))
    return int(np.max(np.abs(i - j)))


def ordering_study() -> pd.DataFrame:
    """Does a fill-reducing ordering help or hurt the measurement cost?

    Classical sparse practice applies a fill-reducing permutation before factorizing.
    This asks what that same step does to the Pauli representation. Three symmetric
    permutations of the structural graph of K are compared at fixed N, fixed spatial
    dimension and fixed physics (a symmetric permutation changes no eigenvalue):

      natural  -- the locality-preserving row-major indexing used everywhere else;
      rcm      -- reverse Cuthill-McKee, the fill-reducing ordering;
      random   -- a fixed-seed scramble, included as the worst-case bracket.

    Reported per ordering: the Cholesky fill nnz(L) of M, the density of the whitened
    A, and the nonzero Pauli count. The lumped-mass rows are the controlling case:
    there M is diagonal, so nnz(L) = N and the density of A are both invariant under
    permutation, and any change in N_P is attributable to Pauli locality alone.
    """
    from scipy.linalg import cholesky
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import reverse_cuthill_mckee

    rows = []
    for dim_label, shape in ORDERING_CASES:
        for mass_type in ("lumped", "consistent"):
            M_nat, K_nat = grid_MK(shape, mass_type)
            n = M_nat.shape[0]
            nq = int(round(np.log2(n)))
            perms = {
                "natural": np.arange(n),
                "rcm": np.asarray(reverse_cuthill_mckee(
                    csr_matrix(_pattern(K_nat).astype(np.int8)), symmetric_mode=True)),
                "random": np.random.default_rng(ORDERING_SEED).permutation(n),
            }
            for name, p in perms.items():
                M, K = M_nat[np.ix_(p, p)], K_nat[np.ix_(p, p)]
                L = cholesky(M, lower=True, check_finite=False)
                A, _, _, _ = mass_whiten(M, K)
                lam = float(np.max(np.abs(eigh(A, eigvals_only=True))))
                pdf = pauli_coefficients(A / lam, tol=1e-10)
                rows.append({
                    "dimension": dim_label, "shape": str(shape), "mass_type": mass_type,
                    "dof": n, "qubits": nq, "ordering": name,
                    "bandwidth_K": _bandwidth(K),
                    "nnz_L": int(np.count_nonzero(_pattern(L))),
                    "density_A": float(np.count_nonzero(_pattern(A))) / A.size,
                    "pauli_terms": int(len(pdf)),
                    "pauli_bound_lumped": (nq + 2) * 2 ** (nq - 1),
                    "pauli_bound_realsym": 2 ** (2 * nq - 1) + 2 ** (nq - 1),
                })
                print(f"[ordering] {dim_label} {mass_type:10s} {name:7s} "
                      f"nnz_L={rows[-1]['nnz_L']:6d} density_A={rows[-1]['density_A']:.4f} "
                      f"N_P={len(pdf):7d}")
    df = pd.DataFrame(rows)
    base = df[df.ordering == "natural"].set_index(["dimension", "mass_type"])["pauli_terms"]
    df["pauli_ratio_to_natural"] = [
        r.pauli_terms / base.loc[(r.dimension, r.mass_type)] for r in df.itertuples()]
    df.to_csv(DATA / "ordering_study.csv", index=False)
    return df


N_SIDE_NEARDEG = 3          # 3x3 grid, corner pinned -> N = 8 = 2^3 free DOFs
PROBE_FRACTION = 0.01       # fixed 1% stiffness loss on one spring
CANONICAL_EDGE = 4          # interior axis-0 spring joining nodes (1,1) and (2,1)
ASYMMETRIES = (0.001, 0.002, 0.005, 0.01, 0.02, 0.05)


def near_degenerate_sweep() -> pd.DataFrame:
    """The informative degeneracy case: a *near*-degenerate pair.

    degeneracy_demo() rotates an exactly degenerate pair inside its own subspace,
    where MAC = cos^2(theta) and eps_sub = 0 hold analytically, so it confirms the
    implementation rather than discovering behaviour. Here the pair is split by a
    physical asymmetry (axis-0 springs stiffened by `delta`) and a *fixed* small
    model discrepancy is applied at every gap. Holding the discrepancy constant
    while the gap varies isolates the effect of the gap alone on each diagnostic
    and on the Davis-Kahan certificate.

    N = 8 = 2^3, so no spectral padding is needed and G_r reduces to the physical
    gap eps_{r+1} - eps_r. The reflection i<->j is an exact symmetry of both K and
    M and acts as diag(-1, +1) on the pair, so the degeneracy is symmetry-protected
    rather than accidental. Every probe edge is swept, so the reported behaviour
    cannot be an artefact of where the discrepancy was placed.
    """
    n = N_SIDE_NEARDEG
    M0_, K0_ = square_uniform_MK(n, "consistent")
    lam0 = eigh(K0_, M0_, eigvals_only=True)
    r = int(np.argmin(np.diff(lam0) / lam0[1:]))
    n_edges = 2 * n * (n - 1)
    print(f"[near-deg] N={len(lam0)}  pair=({r},{r+1})  "
          f"unperturbed rel gap={(lam0[r+1]-lam0[r])/lam0[r+1]:.2e}")

    rows = []
    for delta in ASYMMETRIES:
        for edge in range(n_edges):
            M_, K_ = square_uniform_MK(n, "consistent", axis_scale=(1.0 + delta, 1.0))
            _, Kp = square_uniform_MK(n, "consistent", axis_scale=(1.0 + delta, 1.0),
                                      probe=(edge, PROBE_FRACTION))
            A, _, _, _ = mass_whiten(M_, K_)
            Ap, _, _, _ = mass_whiten(M_, Kp)      # probe leaves M, hence L, unchanged
            scale = float(np.max(eigh(A, eigvals_only=True)))
            H, Hp = A / scale, Ap / scale
            eps, U = eigh(H)
            _, Uh = eigh(Hp)
            norm_delta = float(np.max(np.abs(eigh(Hp - H, eigvals_only=True))))

            G_r = float(eps[r + 1] - eps[r])       # no padding at N=8, so alpha does not bind
            g_pair = float(min(eps[r] - eps[r - 1], eps[r + 2] - eps[r + 1]))
            Y, Yh = U[:, r:r + 2], Uh[:, r:r + 2]
            sv = np.clip(np.linalg.svd(Y.T @ Yh, compute_uv=False), -1.0, 1.0)
            mac_r = float((U[:, r] @ Uh[:, r]) ** 2)
            rows.append({
                "delta": delta, "probe_edge": edge,
                "canonical": edge == CANONICAL_EDGE,
                "rel_gap": float((eps[r + 1] - eps[r]) / eps[r + 1]),
                "norm_delta": norm_delta, "G_r": G_r, "g_pair": g_pair,
                "mac_r": mac_r,
                "mac_r1": float((U[:, r + 1] @ Uh[:, r + 1]) ** 2),
                "eps_sub": float(np.sqrt(np.sum(1.0 - sv ** 2))),
                "sin_individual": float(np.sqrt(max(1.0 - mac_r, 0.0))),
                "sin_subspace": float(np.sqrt(max(1.0 - sv.min() ** 2, 0.0))),
                "dk_individual": 2.0 * norm_delta / G_r,
                "dk_subspace": 2.0 * norm_delta / g_pair,
                # Weyl's precondition ||Delta|| < g/2 is exactly the condition
                # 2||Delta||/g < 1, i.e. the bound is valid precisely when non-vacuous.
                "weyl_ok_individual": bool(norm_delta < G_r / 2),
                "weyl_ok_subspace": bool(norm_delta < g_pair / 2),
            })
    df = pd.DataFrame(rows)
    df.to_csv(DATA / "near_degenerate_sweep.csv", index=False)
    can = df[df.canonical]
    print(can[["delta", "rel_gap", "mac_r", "eps_sub",
               "dk_individual", "dk_subspace"]].to_string(index=False))
    return df


def _fit_exponents(df_hd: pd.DataFrame) -> dict:
    """Power-law exponents of the consistent-mass Pauli count, per dimension.

    These are quoted in Sec. VIII-B alongside the 1D value from revision_study, and are
    persisted here so that every reported exponent has a generator.
    """
    out = {}
    df_1d = pd.read_csv(PARENT / "revision_data" / "hamiltonian_scaling.csv")
    g1 = df_1d[df_1d.mass_type == "consistent"].sort_values("dof")
    series = {"1D": (g1.dof.to_numpy(float), g1.pauli_terms.to_numpy(float))}
    for dim in ("2D", "3D"):
        g = df_hd[(df_hd.dimension == dim) & (df_hd.mass_type == "consistent")].sort_values("dof")
        series[dim] = (g.dof.to_numpy(float), g.pauli_terms.to_numpy(float))
    for dim, (n, y) in series.items():
        out[f"pauli_exponent_{dim}_full_range"] = float(np.polyfit(np.log(n), np.log(y), 1)[0])
        out[f"n_points_{dim}"] = int(len(n))
    m = series["1D"][0] >= 16
    out["pauli_exponent_1D_n_ge_16"] = float(
        np.polyfit(np.log(series["1D"][0][m]), np.log(series["1D"][1][m]), 1)[0])
    (DATA / "dimension_fit_statistics.json").write_text(json.dumps(out, indent=2))
    for k, v in out.items():
        print(f"[dim-fit] {k} = {v}")
    return out


def main() -> None:
    df_scale = scaling_sweep()
    _comparison_figure(df_scale)
    _fit_exponents(df_scale)
    degeneracy_demo()
    near_degenerate_sweep()
    ordering_study()

    # Mirror the generated figures into the directory the manuscript reads, exactly as
    # revision_study.main() does. Without this the copies under ../figures are manual and
    # silently go stale: they were still Aug-9 Type-3 PDFs after the generators had moved on.
    manuscript_figures = PARENT / "figures"
    manuscript_figures.mkdir(exist_ok=True)
    for source in FIG.glob("*"):
        if source.suffix.lower() in {".pdf", ".png"}:
            (manuscript_figures / source.name).write_bytes(source.read_bytes())


if __name__ == "__main__":
    main()
