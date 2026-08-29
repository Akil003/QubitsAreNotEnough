"""Major-revision numerical studies for the structural VQE manuscript.

The script adds:
1. theorem-driven VQD penalty threshold tests;
2. Gershgorin tightness and Hamiltonian/Pauli scaling for lumped and consistent mass;
3. Cholesky preprocessing cost and density growth;
4. trainability diagnostics through random-initialization gradient variance;
5. larger 8-, 16-, and 32-DOF state-vector VQE benchmarks;
6. finite-shot modal-frequency and modal-strain-energy damage-localization studies;
7. stylized depolarizing and readout-bias error budgets.

No quantum SDK is required.
"""
from __future__ import annotations

import itertools
import json
import math
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.linalg import cholesky, eigh, solve_triangular

from reproduce_study import (RealAmplitudeCircuit, build_shear_building, gershgorin_upper,
                             mac_matrix, mass_normalized_operator, recover_modes,
                             relative_residual, safe_pad, solve_vqd)

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "revision_figures"
DATA = ROOT / "revision_data"
FIG.mkdir(exist_ok=True)
DATA.mkdir(exist_ok=True)
# Shared by finite_shot_frequency, damage_shot_study and stylized_noise_budget, in that
# order, which is the order main() calls them in. Results are deterministic for that
# order but depend on it: calling the damage study first shifts its accuracies by a few
# percentage points. Studies added later should take their own generator, as
# restart_study, theory_beta_study, gradient_trainability and the oracle damage variant do.
RNG = np.random.default_rng(20260623)


def chain_model(n: int, mass_type: str = "lumped", damage_story: int | None = None,
                damage_fraction: float = 0.0) -> tuple[np.ndarray, np.ndarray, list[np.ndarray]]:
    """Fixed-base axial/shear chain with n free translational DOFs.

    The consistent-mass option is assembled from two-node element matrices and is
    intentionally included because mass whitening can densify the transformed operator.
    """
    # Smoothly varying properties prevent exact Toeplitz symmetries from understating counts.
    xi = np.linspace(0.0, 1.0, n)
    k = 2.2e8 * (1.0 - 0.35 * xi)
    me = 2.0e4 * (1.0 - 0.20 * xi)
    if damage_story is not None:
        k = k.copy()
        k[damage_story] *= 1.0 - damage_fraction

    Kfull = np.zeros((n + 1, n + 1))
    Mfull = np.zeros((n + 1, n + 1))
    Ke_full: list[np.ndarray] = []
    for e in range(n):
        ids = [e, e + 1]
        ke2 = k[e] * np.array([[1.0, -1.0], [-1.0, 1.0]])
        if mass_type == "consistent":
            me2 = me[e] / 6.0 * np.array([[2.0, 1.0], [1.0, 2.0]])
        elif mass_type == "lumped":
            me2 = me[e] / 2.0 * np.eye(2)
        else:
            raise ValueError(mass_type)
        Kef = np.zeros_like(Kfull)
        for a in range(2):
            for b in range(2):
                Kfull[ids[a], ids[b]] += ke2[a, b]
                Mfull[ids[a], ids[b]] += me2[a, b]
                Kef[ids[a], ids[b]] += ke2[a, b]
        Ke_full.append(Kef)
    # Eliminate fixed base coordinate.
    K = Kfull[1:, 1:]
    M = Mfull[1:, 1:]
    Ke = [x[1:, 1:] for x in Ke_full]
    return M, K, Ke


def chain_MK(n: int, mass_type: str = "consistent") -> tuple[np.ndarray, np.ndarray]:
    """M and K identical to ``chain_model(n, mass_type)`` but without the
    ``O(N^3)``-memory per-element matrix list, so studies that need only the
    assembled operators (e.g. the trainability sweep) can reach larger n.
    """
    xi = np.linspace(0.0, 1.0, n)
    k = 2.2e8 * (1.0 - 0.35 * xi)
    me = 2.0e4 * (1.0 - 0.20 * xi)
    Kfull = np.zeros((n + 1, n + 1))
    Mfull = np.zeros((n + 1, n + 1))
    for e in range(n):
        ids = (e, e + 1)
        ke2 = k[e] * np.array([[1.0, -1.0], [-1.0, 1.0]])
        if mass_type == "consistent":
            me2 = me[e] / 6.0 * np.array([[2.0, 1.0], [1.0, 2.0]])
        elif mass_type == "lumped":
            me2 = me[e] / 2.0 * np.eye(2)
        else:
            raise ValueError(mass_type)
        for a in range(2):
            for b in range(2):
                Kfull[ids[a], ids[b]] += ke2[a, b]
                Mfull[ids[a], ids[b]] += me2[a, b]
    return Mfull[1:, 1:], Kfull[1:, 1:]


def mass_whiten(M: np.ndarray, K: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, float]:
    t0 = time.perf_counter()
    L = cholesky(M, lower=True, check_finite=False)
    t_chol = time.perf_counter() - t0
    t0 = time.perf_counter()
    left = solve_triangular(L, K, lower=True, check_finite=False)
    A = solve_triangular(L, left.T, lower=True, check_finite=False).T
    A = 0.5 * (A + A.T)
    t_transform = time.perf_counter() - t0
    return A, L, t_chol, t_transform


def gersh(A: np.ndarray) -> float:
    return float(np.max(np.diag(A) + np.sum(np.abs(A), axis=1) - np.abs(np.diag(A))))


def pad(H: np.ndarray, dim: int, alpha: float) -> np.ndarray:
    out = np.zeros((dim, dim))
    out[: H.shape[0], : H.shape[1]] = H
    if dim > H.shape[0]:
        out[H.shape[0]:, H.shape[0]:] = alpha * np.eye(dim - H.shape[0])
    return out


def label_masks(labels: tuple[str, ...]) -> tuple[int, int, int, int]:
    """Return x mask, z mask, nonidentity mask, and Y count for MSB-first labels."""
    x = z = non = ny = 0
    n = len(labels)
    for q, p in enumerate(labels):
        bit = 1 << (n - 1 - q)
        if p in ("X", "Y"):
            x |= bit
        if p in ("Z", "Y"):
            z |= bit
        if p != "I":
            non |= bit
        if p == "Y":
            ny += 1
    return x, z, non, ny


def pauli_coefficients(H: np.ndarray, tol: float = 1e-11) -> pd.DataFrame:
    """Pauli decomposition using signed permutation action, avoiding Kronecker matrices."""
    d = H.shape[0]
    n = int(round(math.log2(d)))
    rows = []
    basis = np.arange(d, dtype=np.int64)
    # There are 4^n Pauli strings but only 2^n distinct z-masks, and the parity vector
    # depends on the z-mask alone, so caching it removes a 2^n-fold repetition of the
    # inner Python loop. Values are unchanged.
    parity_cache: dict[int, np.ndarray] = {}
    for labels in itertools.product("IXYZ", repeat=n):
        xmask, zmask, nonmask, ny = label_masks(labels)
        target = basis ^ xmask
        parity = parity_cache.get(zmask)
        if parity is None:
            parity = np.array([int((int(i) & zmask).bit_count() & 1) for i in basis])
            parity_cache[zmask] = parity
        phase = (1j ** ny) * (1.0 - 2.0 * parity)
        coeff = np.sum(phase * H[basis, target]) / d
        if abs(coeff) > tol:
            rows.append({
                "pauli": "".join(labels),
                "coefficient": float(np.real_if_close(coeff, tol=1000).real),
                "xmask": xmask,
                "zmask": zmask,
                "nonmask": nonmask,
                "weight": int(nonmask.bit_count()),
            })
    # kind="stable" is required, not cosmetic. Each Pauli term draws one binomial from
    # a shared generator in sample_prepared, so the row order fixes which term consumes
    # which draw. Several benchmarks here contain exactly tied |coefficient| values
    # (the six-DOF Hamiltonian has two groups of four), and the default quicksort breaks
    # those ties differently across environments, which silently changes every Monte
    # Carlo result. A stable sort over the deterministic itertools.product generation
    # order makes the term list, and hence the sampling, reproducible.
    return pd.DataFrame(rows).sort_values("coefficient", key=lambda s: np.abs(s),
                                          ascending=False, kind="stable")


def qwc_groups(df: pd.DataFrame) -> int:
    """Greedy qubit-wise commuting grouping."""
    groups: list[list[tuple[int, int, int]]] = []
    for row in df.itertuples(index=False):
        item = (int(row.xmask), int(row.zmask), int(row.nonmask))
        placed = False
        for group in groups:
            ok = True
            for x2, z2, non2 in group:
                conflict = item[2] & non2 & ((item[0] ^ x2) | (item[1] ^ z2))
                if conflict:
                    ok = False
                    break
            if ok:
                group.append(item)
                placed = True
                break
        if not placed:
            groups.append([item])
    return len(groups)


def pauli_expectation(state: np.ndarray, pauli: str) -> float:
    labels = tuple(pauli)
    xmask, zmask, _, ny = label_masks(labels)
    basis = np.arange(state.size, dtype=np.int64)
    target = basis ^ xmask
    parity = np.array([int((int(i) & zmask).bit_count() & 1) for i in basis])
    phase = (1j ** ny) * (1.0 - 2.0 * parity)
    # <psi|P|psi> = sum_j conj(psi_{j xor x}) phase_j psi_j
    val = np.sum(np.conjugate(state[target]) * phase * state[basis])
    return float(np.real_if_close(val).real)



def prepare_operator(df: pd.DataFrame, state: np.ndarray, readout_p: float = 0.0,
                     depolarizing_p: float = 0.0) -> list[tuple[float, float, bool]]:
    prepared = []
    for row in df.itertuples(index=False):
        pstr = row.pauli
        is_identity = set(pstr) == {"I"}
        mu = 1.0 if is_identity else pauli_expectation(state, pstr)
        if not is_identity:
            mu *= (1.0 - depolarizing_p)
            mu *= (1.0 - 2.0 * readout_p) ** int(row.weight)
        prepared.append((float(row.coefficient), float(np.clip(mu, -1.0, 1.0)), is_identity))
    return prepared


def sample_prepared(prepared: list[tuple[float, float, bool]], shots_per_term: int,
                    rng: np.random.Generator) -> float:
    total = 0.0
    for coeff, mu, identity in prepared:
        if identity:
            mhat = 1.0
        else:
            plus = rng.binomial(shots_per_term, (1.0 + mu) / 2.0)
            mhat = 2.0 * plus / shots_per_term - 1.0
        total += coeff * mhat
    return float(total)


def theorem_penalty_study() -> pd.DataFrame:
    M, K, _ = chain_model(6, "lumped")
    A, _, _, _ = mass_whiten(M, K)
    lam, Y = eigh(A)
    s = lam[-1]
    lamn = lam / s
    H = pad(A / s, 8, 1.10)
    # Extend exact physical modes to padded space.
    Yp = np.zeros((8, 6))
    Yp[:6, :] = Y
    rows = []
    for r in range(1, 4):
        threshold = float(max(lamn[r] - lamn[:r]))
        beta_grid = np.linspace(0.0, max(1.25 * threshold, 0.05), 121)
        for beta in beta_grid:
            Hd = H.copy()
            for j in range(r):
                Hd += beta * np.outer(Yp[:, j], Yp[:, j])
            vals, vecs = eigh(Hd)
            v = vecs[:, 0]
            overlap = float(abs(v @ Yp[:, r]) ** 2)
            recovered = int(np.argmax(np.abs(Yp.T @ v) ** 2) + 1)
            rows.append({
                "target_mode": r + 1,
                "beta": beta,
                "threshold": threshold,
                "beta_over_threshold": beta / threshold if threshold > 0 else np.nan,
                "target_overlap": overlap,
                "recovered_mode": recovered,
                "ground_value": vals[0],
            })
    df = pd.DataFrame(rows)
    df.to_csv(DATA / "penalty_threshold.csv", index=False)

    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    for r, g in df.groupby("target_mode"):
        ax.plot(g["beta_over_threshold"], g["target_overlap"], label=f"Target mode {r}")
    ax.axvline(1.0, linestyle="--", linewidth=1.2, label="Theoretical threshold")
    ax.set_xlim(0, 1.25)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel(r"Common deflation weight / exact threshold, $\beta/\beta_{\min}$")
    ax.set_ylabel("Squared overlap with target mode")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "penalty_threshold.pdf", bbox_inches="tight")
    fig.savefig(FIG / "penalty_threshold.png", dpi=250, bbox_inches="tight")
    plt.close(fig)

    # Perturbed-projector stability study for target mode 4.
    r = 3
    threshold = float(max(lamn[r] - lamn[:r]))
    beta = threshold + 0.05
    exact_def = H.copy()
    for j in range(r):
        exact_def += beta * np.outer(Yp[:, j], Yp[:, j])
    exact_vals = eigh(exact_def, eigvals_only=True)
    target_val = lamn[r]
    nearest = int(np.argmin(np.abs(exact_vals - target_val)))
    separation = float(np.min(np.abs(np.delete(exact_vals, nearest) - target_val)))
    perturb_rows = []
    # The grid is refined below g_r/2 so that most points satisfy the Weyl precondition
    # of the Davis-Kahan bound; the coarser tail is retained to show where it lapses.
    angles = np.concatenate([np.linspace(0.0, 0.035, 8), np.linspace(0.04, 0.20, 9)])
    for angle in angles:
        Hhat = H.copy()
        delta = np.zeros_like(H)
        for j in range(r):
            # Mix with a different physical state while preserving normalization. Note
            # that for j = 0 the partner (j + r) % 6 == r is the TARGET mode itself, so
            # the contamination lies along the direction deflation is meant to protect.
            # This is deliberately adversarial and explains the large amplification.
            q = Yp[:, (j + r) % 6]
            vhat = np.cos(angle) * Yp[:, j] + np.sin(angle) * q
            term_hat = beta * np.outer(vhat, vhat)
            term_exact = beta * np.outer(Yp[:, j], Yp[:, j])
            Hhat += term_hat
            delta += term_hat - term_exact
        vals, vecs = eigh(Hhat)
        v = vecs[:, 0]
        sin_theta = float(np.sqrt(max(0.0, 1.0 - abs(v @ Yp[:, r]) ** 2)))
        perturb = float(np.linalg.norm(delta, 2))
        # Eq. (dkbound): sin(theta) <= 2||Delta_r|| / g_r, valid only while
        # ||Delta_r|| < g_r / 2. Outside that range the expression exceeds one and the
        # guarantee is vacuous, which the weyl_condition_met flag records.
        bound = min(1.0, 2.0 * perturb / max(separation, 1e-15))
        perturb_rows.append({"projector_angle_rad": angle, "sin_theta": sin_theta,
                             "operator_perturbation": perturb, "davis_kahan_bound": bound,
                             "weyl_condition_met": bool(perturb < 0.5 * separation),
                             "spectral_separation": separation})
    pd.DataFrame(perturb_rows).to_csv(DATA / "deflation_perturbation.csv", index=False)
    return df


def scaling_and_pauli_study() -> pd.DataFrame:
    rows = []
    for mass_type in ("lumped", "consistent"):
        for n in (4, 8, 16, 32, 64, 128):
            M, K, _ = chain_model(n, mass_type)
            A, _, tc, tt = mass_whiten(M, K)
            lam_max = float(eigh(A, eigvals_only=True, subset_by_index=[n - 1, n - 1])[0])
            sG = gersh(A)
            H = A / lam_max
            t0 = time.perf_counter()
            pdf = pauli_coefficients(H, tol=1e-10)
            tpauli = time.perf_counter() - t0
            groups = qwc_groups(pdf) if n <= 64 else np.nan
            rows.append({
                "mass_type": mass_type,
                "dof": n,
                "qubits": int(round(math.log2(n))),
                "lambda_max": lam_max,
                "gershgorin_bound": sG,
                "gershgorin_ratio": sG / lam_max,
                "nnz_A": int(np.count_nonzero(np.abs(A) > 1e-12 * np.max(np.abs(A)))),
                "density_A": float(np.count_nonzero(np.abs(A) > 1e-12 * np.max(np.abs(A))) / A.size),
                "pauli_terms": int(len(pdf)),
                "qwc_groups": groups,
                "cholesky_seconds": tc,
                "transform_seconds": tt,
                "pauli_seconds": tpauli,
            })
    df = pd.DataFrame(rows)
    df.to_csv(DATA / "hamiltonian_scaling.csv", index=False)

    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    for mt, g in df.groupby("mass_type"):
        ax.loglog(g["dof"], g["pauli_terms"], marker="o", label=f"{mt.capitalize()} mass")
    ax.loglog(df["dof"].unique(), df["dof"].unique() ** 2, linestyle=":", label=r"$N^2$ reference")
    ax.set_xlabel("Structural DOFs, $N$")
    ax.set_ylabel("Nonzero Pauli terms")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "pauli_scaling.pdf", bbox_inches="tight")
    fig.savefig(FIG / "pauli_scaling.png", dpi=250, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    for mt, g in df.groupby("mass_type"):
        ax.semilogx(g["dof"], g["gershgorin_ratio"], marker="o", label=f"{mt.capitalize()} mass")
    ax.axhline(1.0, linestyle="--", linewidth=1.0)
    ax.set_xlabel("Structural DOFs, $N$")
    ax.set_ylabel(r"Gershgorin looseness, $s_G/\lambda_{\max}$")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "gershgorin_tightness.pdf", bbox_inches="tight")
    fig.savefig(FIG / "gershgorin_tightness.png", dpi=250, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    for mt, g in df.groupby("mass_type"):
        ax.loglog(g["dof"], g["density_A"], marker="o", label=f"{mt.capitalize()} mass")
    ax.set_xlabel("Structural DOFs, $N$")
    ax.set_ylabel("Density of mass-whitened operator")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "whitening_density.pdf", bbox_inches="tight")
    fig.savefig(FIG / "whitening_density.png", dpi=250, bbox_inches="tight")
    plt.close(fig)
    return df


def preprocessing_timing() -> pd.DataFrame:
    rows = []
    for mass_type in ("lumped", "consistent"):
        for n in (32, 64, 128, 256, 512):
            trials = []
            for _ in range(3):
                M, K, _ = chain_model(n, mass_type)
                A, L, tc, tt = mass_whiten(M, K)
                trials.append((tc, tt, np.count_nonzero(L), np.count_nonzero(np.abs(A) > 1e-12*np.max(np.abs(A)))))
            arr = np.array(trials, dtype=float)
            rows.append({"mass_type": mass_type, "dof": n,
                         "cholesky_seconds": float(np.median(arr[:, 0])),
                         "transform_seconds": float(np.median(arr[:, 1])),
                         "nnz_L": int(np.median(arr[:, 2])),
                         "nnz_A": int(np.median(arr[:, 3]))})
    df = pd.DataFrame(rows)
    df.to_csv(DATA / "preprocessing_timing.csv", index=False)
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    for mt, g in df.groupby("mass_type"):
        ax.loglog(g["dof"], g["cholesky_seconds"] + g["transform_seconds"], marker="o",
                  label=f"{mt.capitalize()} mass")
    ax.set_xlabel("Structural DOFs, $N$")
    ax.set_ylabel("Dense whitening wall time (s)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "preprocessing_scaling.pdf", bbox_inches="tight")
    fig.savefig(FIG / "preprocessing_scaling.png", dpi=250, bbox_inches="tight")
    plt.close(fig)
    return df


def gradient_trainability() -> pd.DataFrame:
    # Use a dedicated seed so this diagnostic is reproducible independently of
    # the Monte Carlo studies executed before or after it. A separate bootstrap
    # stream leaves the main sampling stream (hence the point estimates) untouched.
    rng = np.random.default_rng(20260623)
    boot = np.random.default_rng(20260624)
    rows = []
    for nqubits in range(2, 13):
        n = 2 ** nqubits
        M, K = chain_MK(n, "consistent")
        A, _, _, _ = mass_whiten(M, K)
        H = A / eigh(A, eigvals_only=True, subset_by_index=[n - 1, n - 1])[0]
        depth = 2
        circ = RealAmplitudeCircuit(nqubits, depth)
        for init in ("uniform", "near_identity"):
            grads = []
            energies = []
            for _ in range(80):
                if init == "uniform":
                    theta = rng.uniform(-np.pi, np.pi, circ.n_parameters)
                else:
                    theta = rng.normal(0.0, 0.10, circ.n_parameters)
                shift = np.zeros_like(theta)
                shift[0] = np.pi / 2
                yp = circ.state(theta + shift)
                ep = float(yp @ H @ yp)
                ym = circ.state(theta - shift)
                em = float(ym @ H @ ym)
                grad = 0.5 * (ep - em)
                grads.append(grad)
                psi = circ.state(theta)
                energies.append(float(psi @ H @ psi))
            # Bootstrap CI for the gradient variance (2000 resamples of the 80 draws).
            g_arr = np.asarray(grads)
            bidx = boot.integers(0, g_arr.size, size=(2000, g_arr.size))
            bvar = g_arr[bidx].var(axis=1, ddof=1)
            ci_lo, ci_hi = (float(x) for x in np.quantile(bvar, [0.025, 0.975]))
            rows.append({"qubits": nqubits, "dof": n, "initialization": init,
                         "gradient_variance": float(np.var(grads, ddof=1)),
                         "gradient_variance_ci_lo": ci_lo,
                         "gradient_variance_ci_hi": ci_hi,
                         "se_log_variance": float(np.std(np.log(bvar), ddof=1)),
                         "mean_abs_gradient": float(np.mean(np.abs(grads))),
                         "energy_variance": float(np.var(energies, ddof=1))})
    df = pd.DataFrame(rows)
    df.to_csv(DATA / "gradient_trainability.csv", index=False)

    # Model discrimination on the uniform-init series.
    # Exponential model: log V = a + b*n_q ; polynomial model: log V = a + p*log(n_q).
    # Same parameter count (2), so R^2 and AIC rank identically; AIC quantifies the gap.
    uni = df[df.initialization == "uniform"].sort_values("qubits")
    nq = uni["qubits"].to_numpy(dtype=float)
    logv = np.log(uni["gradient_variance"].to_numpy())

    def _fit(x: np.ndarray, y: np.ndarray) -> dict:
        coeff, cov = np.polyfit(x, y, 1, cov=True)
        resid = y - np.polyval(coeff, x)
        rss = float(np.sum(resid ** 2))
        m = y.size
        return {"slope": float(coeff[0]), "slope_se": float(np.sqrt(cov[0, 0])),
                "r2": float(1.0 - rss / np.sum((y - y.mean()) ** 2)),
                "aic": float(m * np.log(rss / m) + 2 * 2)}

    def _report(mask: np.ndarray, label: str) -> dict:
        e = _fit(nq[mask], logv[mask])
        p = _fit(np.log(nq[mask]), logv[mask])
        print(f"[gradient] {label} (n={int(mask.sum())}): "
              f"exp R^2={e['r2']:.3f} rate={e['slope']:.3f}/qubit AIC={e['aic']:.1f} | "
              f"poly R^2={p['r2']:.3f} exponent={p['slope']:.3f} AIC={p['aic']:.1f} | "
              f"dAIC(exp-poly)={e['aic'] - p['aic']:+.1f}")
        return {"exp_r2": e["r2"], "exp_slope": e["slope"], "poly_r2": p["r2"],
                "poly_exponent": p["slope"], "delta_aic_exp_minus_poly": e["aic"] - p["aic"]}

    stats: dict[str, float] = {}
    for key, m, lab in (("full_range", nq >= 2, "full range n_q=2-12"),
                        ("nq_ge_4", nq >= 4, "restricted n_q=4-12")):
        stats.update({f"gradient_{key}_{k}": v for k, v in _report(m, lab).items()})
    # Direct flatness test: is the tail (n_q>=7) exponential slope distinguishable from 0?
    tail = _fit(nq[nq >= 7], logv[nq >= 7])
    lo, hi = tail["slope"] - 1.96 * tail["slope_se"], tail["slope"] + 1.96 * tail["slope_se"]
    print(f"[gradient] tail slope n_q=7-12: {tail['slope']:.3f}/qubit "
          f"(95% CI [{lo:.3f}, {hi:.3f}]) -> {'consistent with 0 (flat)' if lo < 0 < hi else 'nonzero'}")
    stats.update({"gradient_tail_slope": tail["slope"],
                  "gradient_tail_ci_lo": lo, "gradient_tail_ci_hi": hi})
    # Sec. VII-E also reports the gradient variance normalized by the spread of the
    # objective itself, which is the statistic that is actually flat.
    ratio = _fit(nq, np.log(uni["gradient_variance"].to_numpy()
                            / uni["energy_variance"].to_numpy()))
    stats.update({"gradient_normalized_slope": ratio["slope"],
                  "gradient_normalized_ci_lo": ratio["slope"] - 1.96 * ratio["slope_se"],
                  "gradient_normalized_ci_hi": ratio["slope"] + 1.96 * ratio["slope_se"],
                  "gradient_normalized_r2": ratio["r2"]})
    print(f"[gradient] normalized Var[g]/Var[E] slope: {ratio['slope']:+.4f}/qubit "
          f"(95% CI [{stats['gradient_normalized_ci_lo']:+.4f}, "
          f"{stats['gradient_normalized_ci_hi']:+.4f}], R^2={ratio['r2']:.3f})")

    e_full, p_full = _fit(nq, logv), _fit(np.log(nq), logv)
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    for init, g in df.groupby("initialization"):
        g = g.sort_values("qubits")
        gv = g["gradient_variance"].to_numpy()
        yerr = np.vstack([gv - g["gradient_variance_ci_lo"].to_numpy(),
                          g["gradient_variance_ci_hi"].to_numpy() - gv])
        ax.errorbar(g["qubits"], gv, yerr=np.clip(yerr, 0, None), marker="o",
                    capsize=2, label=init.replace("_", " ").title())
    ax.semilogy(nq, np.exp(np.polyval([e_full["slope"], np.polyfit(nq, logv, 1)[1]], nq)), "--",
                color="0.45", label=fr"Exponential fit ($R^2={e_full['r2']:.2f}$)")
    ax.semilogy(nq, np.exp(np.polyval(np.polyfit(np.log(nq), logv, 1), np.log(nq))), ":",
                color="0.45", label=fr"Polynomial fit ($R^2={p_full['r2']:.2f}$)")
    ax.set_yscale("log")
    ax.set_xlabel("Qubits")
    ax.set_ylabel("Variance of a parameter-shift gradient")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "gradient_trainability.pdf", bbox_inches="tight")
    fig.savefig(FIG / "gradient_trainability.png", dpi=250, bbox_inches="tight")
    plt.close(fig)
    return df, stats


def trainability_depth_sweep() -> pd.DataFrame:
    """Gradient variance across circuit depth, averaged over parameters.

    The main sweep holds L=2 and differentiates theta_0, which is the least
    plateau-prone parameter in the shallowest useful circuit -- so it could not
    exhibit gradient concentration even if concentration were present. This sweep
    varies the depth and picks the differentiated parameter uniformly at random per
    sample, so the reported variance is a parameter-averaged quantity at deeper
    circuits. It uses its own generator and does not touch the shared stream.
    """
    rng = np.random.default_rng(20260626)
    rows = []
    for nqubits in range(2, 11):
        n = 2 ** nqubits
        M, K = chain_MK(n, "consistent")
        A, _, _, _ = mass_whiten(M, K)
        H = A / eigh(A, eigvals_only=True, subset_by_index=[n - 1, n - 1])[0]
        for depth in (2, 4, 8):
            circ = RealAmplitudeCircuit(nqubits, depth)
            grads = []
            for _ in range(80):
                theta = rng.uniform(-np.pi, np.pi, circ.n_parameters)
                a = int(rng.integers(circ.n_parameters))   # parameter-averaged
                shift = np.zeros_like(theta)
                shift[a] = np.pi / 2
                yp = circ.state(theta + shift)
                ym = circ.state(theta - shift)
                grads.append(0.5 * (float(yp @ H @ yp) - float(ym @ H @ ym)))
            rows.append({"qubits": nqubits, "depth": depth,
                         "parameters": circ.n_parameters,
                         "gradient_variance": float(np.var(grads, ddof=1)),
                         "mean_abs_gradient": float(np.mean(np.abs(grads)))})
        print(f"[depth-sweep] n_q={nqubits}: " + "  ".join(
            f"L={r['depth']} var={r['gradient_variance']:.3e}"
            for r in rows[-3:]))
    df = pd.DataFrame(rows)
    df.to_csv(DATA / "trainability_depth_sweep.csv", index=False)
    for depth, g in df.groupby("depth"):
        g = g.sort_values("qubits")
        slope = float(np.polyfit(g.qubits.to_numpy(float),
                                 np.log(g.gradient_variance.to_numpy()), 1)[0])
        print(f"[depth-sweep] L={depth}: exponential slope {slope:+.3f}/qubit "
              f"(barren-plateau reference -0.693)")
    return df


def larger_vqe_benchmarks() -> pd.DataFrame:
    rows = []
    for n in (8, 16):
        M, K, _ = chain_model(n, "lumped")
        A, _, _, _ = mass_whiten(M, K)
        exact = eigh(A, eigvals_only=True)
        s = exact[-1]
        H = A / s
        for depth in (1, 2, 3):
            result = solve_vqd(H, n_modes=1, depth=depth, physical_dim=n,
                               gamma=0.0, beta=1.05, starts=2, seed=31+n+depth,
                               maxiter=160)
            lam_q = result.eigenvalues_scaled[0] * s
            freq_err = 100.0 * abs(math.sqrt(lam_q) - math.sqrt(exact[0])) / math.sqrt(exact[0])
            rows.append({"dof": n, "qubits": int(math.log2(n)), "depth": depth,
                         "parameters": int(math.log2(n)) * (depth + 1),
                         "frequency_error_pct": freq_err,
                         "variance": result.variances[0],
                         "iterations": result.iterations[0]})
    df = pd.DataFrame(rows)
    df.to_csv(DATA / "larger_vqe_benchmarks.csv", index=False)
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    for n, g in df.groupby("dof"):
        ax.semilogy(g["depth"], np.maximum(g["frequency_error_pct"], 1e-14), marker="o", label=f"{n} DOFs")
    ax.set_xlabel("Entangling depth")
    ax.set_ylabel("Fundamental-frequency error (%)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "larger_vqe_error.pdf", bbox_inches="tight")
    fig.savefig(FIG / "larger_vqe_error.png", dpi=250, bbox_inches="tight")
    plt.close(fig)
    return df


def restart_study() -> pd.DataFrame:
    """Optimizer-restart sensitivity for the fundamental mode of the shear building.

    Twenty independent uniform initializations per depth, one L-BFGS-B run each, so
    that systematic ansatz error is separated from optimizer luck. Initialization k
    is drawn from default_rng(k) inside solve_vqd, so every row is reproducible on
    its own and the study does not touch the shared module generator.
    """
    model = build_shear_building()
    lam, _ = eigh(model.K, model.M)
    f1 = math.sqrt(lam[0]) / (2 * math.pi)
    A, _ = mass_normalized_operator(model.M, model.K)
    s = gershgorin_upper(A)
    H = safe_pad(A / s, target_dim=8, alpha=1.25)
    rows = []
    for depth in (1, 2):
        for seed in range(20):
            res = solve_vqd(H, n_modes=1, depth=depth, physical_dim=6,
                            gamma=0.0, starts=1, seed=seed)
            lam_q = float(res.eigenvalues_scaled[0]) * s
            rows.append({
                "depth": depth,
                "seed": seed,
                "frequency_error_pct":
                    100.0 * abs(math.sqrt(max(lam_q, 0.0)) / (2 * math.pi) - f1) / f1,
                "iterations": int(res.iterations[0]),
                "variance": float(res.variances[0]),
            })
    df = pd.DataFrame(rows)
    df.to_csv(DATA / "optimizer_restarts.csv", index=False)
    for depth, g in df.groupby("depth"):
        e = g.frequency_error_pct
        print(f"[restarts] depth {depth}: min={e.min():.6g} median={e.median():.6g} "
              f"max={e.max():.6g}")
    return df


def theory_beta_study() -> pd.DataFrame:
    """Execute Algorithm 1's certified deflation rule end to end.

    Uses the quantum-only Gershgorin bound policy U_r = U_G / s_t, which requires no
    auxiliary classical eigensolve. For target r the penalty on previously extracted
    mode j is

        beta_j = U_r - eps_hat_j + delta_j + tau,

    where eps_hat_j is the energy actually returned for mode j and delta_j its
    residual radius sqrt(sigma_H^2) -- the exact-state-vector case of the
    residual-plus-shot rule, with no sampling term. Since s_t = U_G here, U_r = 1 and
    beta_j > eps_r - eps_j holds for every physical eigenvalue, so the threshold of
    Theorem 1 is satisfied by construction rather than by tuning.
    """
    model = build_shear_building()
    lam, Phi = eigh(model.K, model.M)
    freqs = np.sqrt(lam) / (2 * np.pi)
    A, L = mass_normalized_operator(model.M, model.K)
    s = gershgorin_upper(A)
    H = safe_pad(A / s, target_dim=8, alpha=1.25)
    U_r = gershgorin_upper(A) / s          # = 1.0 under the Gershgorin policy
    tau = 1e-3

    used_betas: dict[int, list[float]] = {}
    current_depth = 0

    def schedule(r: int, j: int, energies: list[float], variances: list[float]) -> float:
        delta_j = math.sqrt(max(variances[j], 0.0))
        beta_j = U_r - energies[j] + delta_j + tau
        used_betas.setdefault(current_depth, []).append(beta_j)
        return beta_j

    rows = []
    for depth in (1, 2):
        current_depth = depth
        res = solve_vqd(H, n_modes=4, depth=depth, physical_dim=6, gamma=0.0,
                        starts=4, seed=17, beta_schedule=schedule)
        lam_q = res.eigenvalues_scaled * s
        freq_q = np.sqrt(np.clip(lam_q, 0.0, None)) / (2 * np.pi)
        Phi_q = recover_modes(res.states, L, 6, model.M)
        MAC = mac_matrix(Phi[:, :4], Phi_q[:, :4], model.M)
        for i in range(4):
            rows.append({
                "depth": depth, "mode": i + 1,
                "f_exact": freqs[i], "f_q": freq_q[i],
                "err": 100.0 * abs(freq_q[i] - freqs[i]) / freqs[i],
                "mac": MAC[i, i],
                "var": res.variances[i],
                "overlap": res.max_overlaps[i],
                "leak": res.padding_leakage[i],
                "res": relative_residual(model.K, model.M, Phi_q[:, i], lam_q[i]),
            })
    df = pd.DataFrame(rows)
    df["beta_min"] = df.depth.map({d: min(v) for d, v in used_betas.items()})
    df["beta_max"] = df.depth.map({d: max(v) for d, v in used_betas.items()})
    df.to_csv(DATA / "six_dof_theory_beta.csv", index=False)
    for depth, g in df.groupby("depth"):
        print(f"[theory-beta] depth {depth}: max freq err={g.err.max():.6g}% "
              f"min MAC={g.mac.min():.12f} max overlap={g.overlap.max():.3e}")
    for d, v in sorted(used_betas.items()):
        print(f"[theory-beta] depth {d} certified penalties span "
              f"{min(v):.4f} to {max(v):.4f} (U_r={U_r:.4f}, tau={tau})")
    return df


def exact_modal_data(n: int = 6, mass_type: str = "lumped", damage_story: int | None = None,
                     damage_fraction: float = 0.0):
    M, K, Ke = chain_model(n, mass_type, damage_story, damage_fraction)
    A, L, _, _ = mass_whiten(M, K)
    lam, Y = eigh(A)
    s = lam[-1]
    # Certified padding (Prop. 2): alpha must exceed U = U_G / s_t. The literal 1.10
    # already satisfies this for the benchmarks used here (U = 1.0920 for the six-DOF
    # lumped chain), so no reported number changes; the max() keeps the certificate
    # valid if the model changes -- the N=4 chain, for instance, has U = 1.1102.
    alpha = max(1.10, gersh(A) / s + 1e-3)
    H = pad(A / s, 2 ** math.ceil(math.log2(n)), alpha)
    Yp = np.zeros((H.shape[0], n))
    Yp[:n, :] = Y
    return M, K, Ke, A, L, lam, Yp, s, H


def finite_shot_frequency() -> pd.DataFrame:
    _, _, _, _, _, lam, Yp, s, H = exact_modal_data()
    pdf = pauli_coefficients(H, tol=1e-12)
    rows = []
    for shots in (100, 1000, 10000, 100000):
        for mode in range(4):
            estimates = []
            prepared = prepare_operator(pdf, Yp[:, mode])
            for _ in range(400):
                ehat = sample_prepared(prepared, shots, RNG)
                fhat = math.sqrt(max(s * ehat, 0.0)) / (2 * math.pi)
                estimates.append(fhat)
            f0 = math.sqrt(lam[mode]) / (2 * math.pi)
            estimates = np.array(estimates)
            rows.append({"shots_per_pauli": shots, "mode": mode + 1, "f_exact": f0,
                         "mean_f": float(np.mean(estimates)), "std_f": float(np.std(estimates, ddof=1)),
                         "rmse_pct": float(100*np.sqrt(np.mean((estimates-f0)**2))/f0),
                         "bias_pct": float(100*(np.mean(estimates)-f0)/f0)})
    df = pd.DataFrame(rows)
    df.to_csv(DATA / "finite_shot_frequency.csv", index=False)
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    for mode, g in df.groupby("mode"):
        ax.loglog(g["shots_per_pauli"], g["rmse_pct"], marker="o", label=f"Mode {mode}")
    ax.set_xlabel("Shots per Pauli term")
    ax.set_ylabel("Frequency RMSE (%)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "finite_shot_frequency.pdf", bbox_inches="tight")
    fig.savefig(FIG / "finite_shot_frequency.png", dpi=250, bbox_inches="tight")
    plt.close(fig)
    return df


def transformed_element_operators(Ke: list[np.ndarray], L: np.ndarray, s: float, dim: int) -> list[np.ndarray]:
    out = []
    for ke in Ke:
        left = solve_triangular(L, ke, lower=True, check_finite=False)
        Ae = solve_triangular(L, left.T, lower=True, check_finite=False).T
        Ae = 0.5 * (Ae + Ae.T) / s
        out.append(pad(Ae, dim, 0.0))
    return out


def damage_shot_study(convention: str = "baseline",
                      rng: np.random.Generator | None = None,
                      make_figure: bool = True) -> pd.DataFrame:
    """Finite-shot modal-strain-energy damage localization.

    ``convention="baseline"`` builds the element operators measured on the damaged
    structure from the BASELINE element stiffnesses Ke0. This is the standard
    modal-strain-energy convention (Stubbs 1995, Cornwell 1999) and the only
    deployable one, since L^-1 Ke^d L^-T cannot be assembled without already knowing
    the damage. It is what the manuscript reports.

    ``convention="oracle"`` instead uses the damaged element stiffnesses. That
    presumes the quantity being inferred and is not implementable; it is retained
    only to reproduce the sensitivity case quoted in the manuscript's Remark 3,
    where the explicit (1 - d_e) numerator factor is shown to cancel most of the
    localization signal.
    """
    if convention not in ("baseline", "oracle"):
        raise ValueError(convention)
    # The oracle variant uses a dedicated generator so it is reproducible on its own
    # and does not perturb the shared stream that produces the reported results.
    if rng is None:
        rng = RNG if convention == "baseline" else np.random.default_rng(20260625)
    base = exact_modal_data()
    _, _, Ke0, _, L0, lam0, Y0, s0, H0 = base
    elem0 = transformed_element_operators(Ke0, L0, s0, H0.shape[0])
    elem0_p = [pauli_coefficients(x, tol=1e-12) for x in elem0]
    H0p = pauli_coefficients(H0, tol=1e-12)
    rows = []
    for severity in (0.02, 0.05, 0.10, 0.20):
        damaged = exact_modal_data(damage_story=2, damage_fraction=severity)
        _, _, Ked, _, Ld, lamd, Yd, sd, Hd = damaged
        # Numerator operators for the damaged epoch. Under "baseline" only the
        # normalization s_d and the state Yd are damaged quantities; the denominator
        # <H_d> is the measured eigenvalue, so the ratio stays consistently scaled.
        Ke_numerator = Ke0 if convention == "baseline" else Ked
        elemd = transformed_element_operators(Ke_numerator, Ld, sd, Hd.shape[0])
        elemd_p = [pauli_coefficients(x, tol=1e-12) for x in elemd]
        Hdp = pauli_coefficients(Hd, tol=1e-12)
        prep_H0 = [prepare_operator(H0p, Y0[:, mode]) for mode in range(4)]
        prep_Hd = [prepare_operator(Hdp, Yd[:, mode]) for mode in range(4)]
        prep_e0 = [[prepare_operator(elem0_p[e], Y0[:, mode]) for e in range(6)] for mode in range(4)]
        prep_ed = [[prepare_operator(elemd_p[e], Yd[:, mode]) for e in range(6)] for mode in range(4)]
        for shots in (100, 1000, 10000, 100000):
            correct = 0
            margins = []
            for _ in range(200):
                eta0 = np.zeros((6, 4))
                etad = np.zeros((6, 4))
                for mode in range(4):
                    den0 = sample_prepared(prep_H0[mode], shots, rng)
                    dend = sample_prepared(prep_Hd[mode], shots, rng)
                    for e in range(6):
                        eta0[e, mode] = sample_prepared(prep_e0[mode][e], shots, rng) / max(den0, 1e-10)
                        etad[e, mode] = sample_prepared(prep_ed[mode][e], shots, rng) / max(dend, 1e-10)
                di = 0.25 * np.sum(np.abs(etad - eta0), axis=1)   # w_i = 1/4, sum_i w_i = 1
                pred = int(np.argmax(di))
                correct += int(pred == 2)
                others = np.delete(di, 2)
                margins.append((di[2] - np.max(others)) / max(abs(di[2]), 1e-12))
            rows.append({"damage_fraction": severity, "shots_per_pauli": shots,
                         "localization_accuracy": correct / 200.0,
                         "median_normalized_margin": float(np.median(margins)),
                         "p10_margin": float(np.quantile(margins, 0.10))})
    df = pd.DataFrame(rows)
    suffix = "" if convention == "baseline" else "_oracle_convention"
    df.to_csv(DATA / f"damage_shot_localization{suffix}.csv", index=False)
    if not make_figure:
        return df
    pivot = df.pivot(index="damage_fraction", columns="shots_per_pauli", values="localization_accuracy")
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    im = ax.imshow(pivot.values, vmin=0, vmax=1, aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(pivot.columns)), [f"{x:g}" for x in pivot.columns])
    ax.set_yticks(range(len(pivot.index)), [f"{100*x:.0f}%" for x in pivot.index])
    ax.set_xlabel("Shots per Pauli term")
    ax.set_ylabel("Stiffness loss")
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = pivot.values[i, j]
            ax.text(j, i, f"{val:.2f}", ha="center", va="center")
    fig.colorbar(im, ax=ax, label="Correct localization probability")
    fig.tight_layout()
    fig.savefig(FIG / "damage_shot_localization.pdf", bbox_inches="tight")
    fig.savefig(FIG / "damage_shot_localization.png", dpi=250, bbox_inches="tight")
    plt.close(fig)
    return df


def stylized_noise_budget() -> pd.DataFrame:
    _, _, _, _, _, lam, Yp, s, H = exact_modal_data()
    pdf = pauli_coefficients(H, tol=1e-12)
    mode = 0
    f0 = math.sqrt(lam[mode]) / (2*math.pi)
    rows = []
    for depol in (0.0, 0.001, 0.005, 0.01, 0.02):
        for readout in (0.0, 0.002, 0.01, 0.02):
            vals = []
            prepared = prepare_operator(pdf, Yp[:, mode], readout_p=readout, depolarizing_p=depol)
            for _ in range(300):
                e = sample_prepared(prepared, 10000, RNG)
                vals.append(math.sqrt(max(s*e, 0.0))/(2*math.pi))
            vals = np.array(vals)
            rows.append({"depolarizing_probability": depol, "readout_flip_probability": readout,
                         "frequency_bias_pct": 100*(np.mean(vals)-f0)/f0,
                         "frequency_rmse_pct": 100*np.sqrt(np.mean((vals-f0)**2))/f0,
                         "frequency_std_pct": 100*np.std(vals, ddof=1)/f0})
    df = pd.DataFrame(rows)
    df.to_csv(DATA / "stylized_noise_budget.csv", index=False)
    return df


def _scaling_exponents(scaling: pd.DataFrame) -> dict:
    """Power-law exponents of the consistent-mass Pauli count quoted in the text."""
    g = scaling[scaling.mass_type == "consistent"].sort_values("dof")
    n, y = g.dof.to_numpy(float), g.pauli_terms.to_numpy(float)
    fit = lambda m: float(np.polyfit(np.log(n[m]), np.log(y[m]), 1)[0])
    return {"pauli_exponent_full_range": fit(n >= 4),
            "pauli_exponent_n_ge_16": fit(n >= 16),
            "pauli_local_exponent_last_doubling": float(np.log2(y[-1] / y[-2]))}


def _wilson_intervals(damage: pd.DataFrame, trials: int = 200) -> dict:
    """Wilson score intervals for the localization accuracies quoted in the text."""
    z = 1.959963984540054
    out = {}
    for sev, shots in ((0.10, 100000), (0.10, 10000), (0.20, 100000), (0.02, 100000)):
        row = damage[(damage.damage_fraction == sev) & (damage.shots_per_pauli == shots)]
        if row.empty:
            continue
        ph = float(row.localization_accuracy.iloc[0])
        d = 1.0 + z * z / trials
        c = (ph + z * z / (2 * trials)) / d
        h = z * math.sqrt(ph * (1 - ph) / trials + z * z / (4 * trials * trials)) / d
        tag = f"wilson_{int(sev * 100)}pct_{shots}"
        out[f"{tag}_lo"] = 100.0 * max(c - h, 0.0)
        out[f"{tag}_hi"] = 100.0 * min(c + h, 1.0)
    return out


def main() -> None:
    # Recompute every dataset used by the revised manuscript. Stochastic
    # studies use the module-level fixed seed, while the gradient diagnostic
    # uses its own fixed seed so that execution order does not change results.
    theorem_penalty_study()
    scaling = scaling_and_pauli_study()
    timing = preprocessing_timing()
    gradient, gradient_stats = gradient_trainability()
    depth_sweep = trainability_depth_sweep()   # own generator; shared RNG untouched
    larger = larger_vqe_benchmarks()
    restarts = restart_study()        # own per-seed generators; shared RNG untouched
    theory = theory_beta_study()      # own generator inside solve_vqd
    shots = finite_shot_frequency()
    damage = damage_shot_study()
    # Remark 3 sensitivity case; dedicated generator, so the shared stream above and
    # therefore every reported number is unaffected by running it.
    damage_shot_study("oracle", make_figure=False)
    noise = stylized_noise_budget()

    summary = {
        "gershgorin_ratio_range": [float(scaling.gershgorin_ratio.min()), float(scaling.gershgorin_ratio.max())],
        "max_pauli_terms": int(scaling.pauli_terms.max()),
        "consistent_mass_density_at_128": float(scaling[(scaling.mass_type == 'consistent') & (scaling.dof == 128)].density_A.iloc[0]),
        "lumped_mass_density_at_128": float(scaling[(scaling.mass_type == 'lumped') & (scaling.dof == 128)].density_A.iloc[0]),
        "gradient_variance_ratio_7q": float(
            gradient[(gradient.qubits == 7) & (gradient.initialization == 'near_identity')].gradient_variance.iloc[0] /
            gradient[(gradient.qubits == 7) & (gradient.initialization == 'uniform')].gradient_variance.iloc[0]
        ),
        "best_16dof_vqe_error_pct": float(larger[larger.dof == 16].frequency_error_pct.min()),
        "damage_accuracy_10pct_10000shots": float(damage[(damage.damage_fraction == 0.10) & (damage.shots_per_pauli == 10000)].localization_accuracy.iloc[0]),
        "mode1_shot_rmse_10000": float(shots[(shots["mode"] == 1) & (shots.shots_per_pauli == 10000)].rmse_pct.iloc[0]),
        "timing_max_seconds": float((timing.cholesky_seconds + timing.transform_seconds).max()),
        "noise_mode1_max_rmse_pct": float(noise.frequency_rmse_pct.max()),
        "restart_depth1_min_error_pct": float(restarts[restarts.depth == 1].frequency_error_pct.min()),
        "restart_depth1_median_error_pct": float(restarts[restarts.depth == 1].frequency_error_pct.median()),
        "restart_depth1_max_error_pct": float(restarts[restarts.depth == 1].frequency_error_pct.max()),
        "depth_sweep_slope_L2": float(np.polyfit(
            depth_sweep[depth_sweep.depth == 2].qubits.to_numpy(float),
            np.log(depth_sweep[depth_sweep.depth == 2].gradient_variance.to_numpy()), 1)[0]),
        "depth_sweep_slope_L8": float(np.polyfit(
            depth_sweep[depth_sweep.depth == 8].qubits.to_numpy(float),
            np.log(depth_sweep[depth_sweep.depth == 8].gradient_variance.to_numpy()), 1)[0]),
        "restart_depth2_median_error_pct": float(restarts[restarts.depth == 2].frequency_error_pct.median()),
        "certified_beta_min_depth2": float(theory[theory.depth == 2].beta_min.iloc[0]),
        "certified_beta_max_depth2": float(theory[theory.depth == 2].beta_max.iloc[0]),
        "certified_beta_max_depth1": float(theory[theory.depth == 1].beta_max.iloc[0]),
        **gradient_stats,
        **_scaling_exponents(scaling),
        **_wilson_intervals(damage),
    }
    (DATA / "revision_summary.json").write_text(json.dumps(summary, indent=2))

    # Synchronize the generated major-revision figures with the directory read
    # by the LaTeX manuscript. Original exact-state figures created by
    # reproduce_study.py remain untouched.
    manuscript_figures = ROOT / "figures"
    manuscript_figures.mkdir(exist_ok=True)
    for source in FIG.glob("*"):
        if source.suffix.lower() in {".pdf", ".png"}:
            (manuscript_figures / source.name).write_bytes(source.read_bytes())

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
