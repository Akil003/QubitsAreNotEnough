"""Reproducible state-vector proof-of-concept for the paper.

No quantum SDK is required. The script emulates a real-amplitude, hardware-efficient
variational circuit, solves a padded mass-normalized structural eigenproblem with VQD,
and evaluates modal and damage-sensitive quantities.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.linalg import cholesky, eigh, solve_triangular
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
DATA = ROOT / "data"
FIG.mkdir(exist_ok=True)
DATA.mkdir(exist_ok=True)

np.set_printoptions(precision=8, suppress=True)


@dataclass
class ModalModel:
    M: np.ndarray
    K: np.ndarray
    Ke: list[np.ndarray]
    masses: np.ndarray
    story_stiffness: np.ndarray


def build_shear_building(damage_story: int | None = None, damage_fraction: float = 0.0) -> ModalModel:
    """Create a six-storey shear-building benchmark in SI units.

    Floor masses decrease modestly with height. Story stiffnesses also decrease with
    height. ``damage_story`` is one-based in the paper but zero-based here.
    """
    masses = np.array([220, 210, 200, 190, 180, 170.0]) * 1e3  # kg
    stiffness = np.array([220, 210, 200, 190, 180, 170.0]) * 1e6  # N/m
    if damage_story is not None:
        stiffness = stiffness.copy()
        stiffness[damage_story] *= 1.0 - damage_fraction

    n = len(masses)
    M = np.diag(masses)
    K = np.zeros((n, n))
    Ke: list[np.ndarray] = []

    for s, ks in enumerate(stiffness):
        e = np.zeros((n, n))
        e[s, s] += ks
        if s > 0:
            e[s - 1, s - 1] += ks
            e[s, s - 1] -= ks
            e[s - 1, s] -= ks
        K += e
        Ke.append(e)
    return ModalModel(M=M, K=K, Ke=Ke, masses=masses, story_stiffness=stiffness)


def mass_normalized_operator(M: np.ndarray, K: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    L = cholesky(M, lower=True)
    left = solve_triangular(L, K, lower=True)
    A = solve_triangular(L, left.T, lower=True).T
    A = 0.5 * (A + A.T)
    return A, L


def gershgorin_upper(A: np.ndarray) -> float:
    return float(max(A[i, i] + np.sum(np.abs(A[i, :])) - abs(A[i, i]) for i in range(A.shape[0])))


def safe_pad(Abar: np.ndarray, target_dim: int, alpha: float = 1.25) -> np.ndarray:
    if target_dim < Abar.shape[0]:
        raise ValueError("target_dim must not be smaller than the physical dimension")
    H = np.zeros((target_dim, target_dim))
    n = Abar.shape[0]
    H[:n, :n] = Abar
    if target_dim > n:
        H[n:, n:] = alpha * np.eye(target_dim - n)
    return H


def cnot_matrix(control: int, target: int, n_qubits: int) -> np.ndarray:
    dim = 2**n_qubits
    U = np.zeros((dim, dim))
    for idx in range(dim):
        bits = [(idx >> (n_qubits - 1 - j)) & 1 for j in range(n_qubits)]
        if bits[control]:
            bits[target] ^= 1
        out = 0
        for bit in bits:
            out = (out << 1) | bit
        U[out, idx] = 1.0
    return U


def ry(theta: float) -> np.ndarray:
    return np.array(
        [[np.cos(theta / 2.0), -np.sin(theta / 2.0)],
         [np.sin(theta / 2.0), np.cos(theta / 2.0)]]
    )


def kron_all(mats: list[np.ndarray]) -> np.ndarray:
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


class RealAmplitudeCircuit:
    def __init__(self, n_qubits: int, depth: int):
        self.n_qubits = n_qubits
        self.depth = depth
        self.n_parameters = n_qubits * (depth + 1)
        edges = [(q, q + 1) for q in range(n_qubits - 1)]
        if n_qubits > 2:
            edges.append((n_qubits - 1, 0))
        self.entanglers = [cnot_matrix(c, t, n_qubits) for c, t in edges]

    def state(self, theta: np.ndarray) -> np.ndarray:
        psi = np.zeros(2**self.n_qubits)
        psi[0] = 1.0
        p = 0
        for _ in range(self.depth):
            Urot = kron_all([ry(theta[p + q]) for q in range(self.n_qubits)])
            p += self.n_qubits
            psi = Urot @ psi
            for Ucnot in self.entanglers:
                psi = Ucnot @ psi
        Urot = kron_all([ry(theta[p + q]) for q in range(self.n_qubits)])
        return Urot @ psi


@dataclass
class VQDResult:
    eigenvalues_scaled: np.ndarray
    states: np.ndarray
    variances: np.ndarray
    max_overlaps: np.ndarray
    padding_leakage: np.ndarray
    iterations: np.ndarray


def solve_vqd(
    H: np.ndarray,
    n_modes: int,
    depth: int,
    physical_dim: int,
    gamma: float = 0.0,
    beta: float = 3.0,
    starts: int = 8,
    seed: int = 12,
    maxiter: int = 450,
) -> VQDResult:
    n_qubits = int(round(np.log2(H.shape[0])))
    circuit = RealAmplitudeCircuit(n_qubits, depth)
    H2 = H @ H
    rng = np.random.default_rng(seed)
    states: list[np.ndarray] = []
    eigenvalues: list[float] = []
    variances: list[float] = []
    overlaps: list[float] = []
    leakage: list[float] = []
    iterations: list[int] = []

    for mode in range(n_modes):
        best: tuple[float, float, float, float, float, int, np.ndarray] | None = None
        for _ in range(starts):
            x0 = rng.uniform(-np.pi, np.pi, circuit.n_parameters)

            def objective(x: np.ndarray) -> float:
                y = circuit.state(x)
                energy = float(y @ H @ y)
                variance = float(y @ H2 @ y - energy**2)
                deflation = sum(float(y @ p) ** 2 for p in states)
                return energy + gamma * max(variance, 0.0) + beta * deflation

            res = minimize(
                objective,
                x0,
                method="L-BFGS-B",
                options={"maxiter": maxiter, "ftol": 1e-13, "gtol": 1e-9, "maxls": 60},
            )
            y = circuit.state(res.x)
            energy = float(y @ H @ y)
            variance = float(max(y @ H2 @ y - energy**2, 0.0))
            max_overlap = max([float(y @ p) ** 2 for p in states] or [0.0])
            pad = float(np.sum(y[physical_dim:] ** 2))
            record = (objective(res.x), energy, variance, max_overlap, pad, int(res.nit), y)
            if best is None or record[0] < best[0]:
                best = record
        assert best is not None
        states.append(best[-1])
        eigenvalues.append(best[1])
        variances.append(best[2])
        overlaps.append(best[3])
        leakage.append(best[4])
        iterations.append(best[5])

    return VQDResult(
        eigenvalues_scaled=np.array(eigenvalues),
        states=np.column_stack(states),
        variances=np.array(variances),
        max_overlaps=np.array(overlaps),
        padding_leakage=np.array(leakage),
        iterations=np.array(iterations),
    )


def recover_modes(states: np.ndarray, L: np.ndarray, physical_dim: int, M: np.ndarray) -> np.ndarray:
    Y = states[:physical_dim, :]
    Phi = solve_triangular(L.T, Y, lower=False)
    for i in range(Phi.shape[1]):
        norm = np.sqrt(Phi[:, i].T @ M @ Phi[:, i])
        Phi[:, i] /= norm
    return Phi


def mac_matrix(A: np.ndarray, B: np.ndarray, W: np.ndarray | None = None) -> np.ndarray:
    if W is None:
        W = np.eye(A.shape[0])
    out = np.zeros((A.shape[1], B.shape[1]))
    for i in range(A.shape[1]):
        for j in range(B.shape[1]):
            num = abs(A[:, i].T @ W @ B[:, j]) ** 2
            den = (A[:, i].T @ W @ A[:, i]) * (B[:, j].T @ W @ B[:, j])
            out[i, j] = float(num / den)
    return out


def modal_strain_energy_fractions(Phi: np.ndarray, K: np.ndarray, Ke: list[np.ndarray]) -> np.ndarray:
    eta = np.zeros((len(Ke), Phi.shape[1]))
    for e, Ke_i in enumerate(Ke):
        for i in range(Phi.shape[1]):
            den = float(Phi[:, i].T @ K @ Phi[:, i])
            eta[e, i] = float(Phi[:, i].T @ Ke_i @ Phi[:, i]) / den
    return eta


def pauli_decomposition(H: np.ndarray, tol: float = 1e-12) -> pd.DataFrame:
    paulis = {
        "I": np.eye(2),
        "X": np.array([[0, 1], [1, 0]], dtype=complex),
        "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
        "Z": np.array([[1, 0], [0, -1]], dtype=complex),
    }
    n = int(round(np.log2(H.shape[0])))
    rows = []
    import itertools
    for labels in itertools.product(paulis.keys(), repeat=n):
        P = kron_all([paulis[x] for x in labels])
        coeff = np.trace(P.conj().T @ H) / (2**n)
        if abs(coeff) > tol:
            rows.append({"pauli": "".join(labels), "coefficient": float(np.real_if_close(coeff).real)})
    return pd.DataFrame(rows).sort_values("coefficient", key=lambda x: np.abs(x), ascending=False)


def relative_residual(K: np.ndarray, M: np.ndarray, phi: np.ndarray, lam: float) -> float:
    r = K @ phi - lam * M @ phi
    den = np.linalg.norm(K @ phi) + abs(lam) * np.linalg.norm(M @ phi)
    return float(np.linalg.norm(r) / den)


def run_case(model: ModalModel, name: str, starts: int = 4) -> dict:
    lam, Phi = eigh(model.K, model.M)
    freqs = np.sqrt(lam) / (2 * np.pi)
    A, L = mass_normalized_operator(model.M, model.K)
    scale = gershgorin_upper(A)
    H = safe_pad(A / scale, target_dim=8, alpha=1.25)

    standard_shallow = solve_vqd(H, 4, depth=1, physical_dim=6, gamma=0.0, starts=starts, seed=17)
    residual_shallow = solve_vqd(H, 4, depth=1, physical_dim=6, gamma=1.0, starts=starts, seed=17)
    standard_deep = solve_vqd(H, 4, depth=2, physical_dim=6, gamma=0.0, starts=starts, seed=17)

    methods = {
        "VQD L=1": standard_shallow,
        "MR-VQD L=1": residual_shallow,
        "VQD L=2": standard_deep,
    }
    out = {
        "lam_exact": lam,
        "freq_exact": freqs,
        "Phi_exact": Phi,
        "A": A,
        "L": L,
        "scale": scale,
        "H": H,
        "methods": {},
    }

    for label, result in methods.items():
        lam_q = result.eigenvalues_scaled * scale
        freq_q = np.sqrt(np.clip(lam_q, 0.0, None)) / (2 * np.pi)
        Phi_q = recover_modes(result.states, L, 6, model.M)
        MAC = mac_matrix(Phi[:, :4], Phi_q[:, :4], model.M)
        residuals = np.array([relative_residual(model.K, model.M, Phi_q[:, i], lam_q[i]) for i in range(4)])
        out["methods"][label] = {
            "result": result,
            "lam": lam_q,
            "freq": freq_q,
            "Phi": Phi_q,
            "MAC": MAC,
            "residual": residuals,
            "freq_error_pct": 100 * np.abs(freq_q - freqs[:4]) / freqs[:4],
        }

    # Export method table.
    rows = []
    for label, d in out["methods"].items():
        for i in range(4):
            rows.append({
                "case": name,
                "method": label,
                "mode": i + 1,
                "f_exact_hz": freqs[i],
                "f_quantum_hz": d["freq"][i],
                "frequency_error_pct": d["freq_error_pct"][i],
                "MAC_diagonal": d["MAC"][i, i],
                "relative_residual": d["residual"][i],
                "hamiltonian_variance": d["result"].variances[i],
                "padding_leakage": d["result"].padding_leakage[i],
                "iterations": d["result"].iterations[i],
            })
    pd.DataFrame(rows).to_csv(DATA / f"{name}_modal_results.csv", index=False)
    return out


def main() -> None:
    baseline_model = build_shear_building()
    damaged_model = build_shear_building(damage_story=2, damage_fraction=0.10)  # third story
    baseline = run_case(baseline_model, "baseline")
    damaged = run_case(damaged_model, "damaged")

    # Pauli decomposition.
    pauli = pauli_decomposition(baseline["H"])
    pauli.to_csv(DATA / "pauli_decomposition.csv", index=False)

    # Modal strain-energy damage index using first four modes.
    eta0_exact = modal_strain_energy_fractions(baseline["Phi_exact"][:, :4], baseline_model.K, baseline_model.Ke)
    etad_exact = modal_strain_energy_fractions(damaged["Phi_exact"][:, :4], damaged_model.K, damaged_model.Ke)
    DI_exact = np.sum(np.abs(etad_exact - eta0_exact), axis=1)

    phi0_q = baseline["methods"]["VQD L=2"]["Phi"][:, :4]
    phid_q = damaged["methods"]["VQD L=2"]["Phi"][:, :4]
    eta0_q = modal_strain_energy_fractions(phi0_q, baseline_model.K, baseline_model.Ke)
    etad_q = modal_strain_energy_fractions(phid_q, damaged_model.K, damaged_model.Ke)
    DI_q = np.sum(np.abs(etad_q - eta0_q), axis=1)

    damage_df = pd.DataFrame({
        "story": np.arange(1, 7),
        "exact_damage_index": DI_exact,
        "quantum_damage_index": DI_q,
    })
    damage_df.to_csv(DATA / "damage_index.csv", index=False)

    # Sensitivity matrix S_ie = -phi_i^T K_e phi_i / lambda_i.
    S = np.zeros((4, 6))
    for i in range(4):
        for e in range(6):
            S[i, e] = -(baseline["Phi_exact"][:, i].T @ baseline_model.Ke[e] @ baseline["Phi_exact"][:, i]) / baseline["lam_exact"][i]
    pd.DataFrame(S, index=[f"Mode {i}" for i in range(1, 5)], columns=[f"Story {i}" for i in range(1, 7)]).to_csv(DATA / "damage_sensitivity.csv")

    # Figure 1: workflow.
    fig, ax = plt.subplots(figsize=(12, 3.2))
    ax.axis("off")
    boxes = [
        (0.02, "Finite-element model\n$\\mathbf{M},\\mathbf{K}$"),
        (0.22, "Mass normalization\n$\\mathbf{A}=\\mathbf{L}^{-1}\\mathbf{K}\\mathbf{L}^{-T}$"),
        (0.43, "Safe padding and\nPauli decomposition"),
        (0.64, "VQE / VQD / SSVQE\nwith residual control"),
        (0.84, "Modal validation and\ndamage indicators"),
    ]
    for x, text in boxes:
        ax.text(x + 0.07, 0.52, text, ha="center", va="center", fontsize=11,
                bbox=dict(boxstyle="round,pad=0.55", fc="white", ec="black"), transform=ax.transAxes)
    for x in [0.18, 0.39, 0.60, 0.80]:
        ax.annotate("", xy=(x + 0.035, 0.52), xytext=(x - 0.015, 0.52),
                    arrowprops=dict(arrowstyle="->", lw=1.5), xycoords=ax.transAxes)
    fig.tight_layout()
    fig.savefig(FIG / "workflow.pdf", bbox_inches="tight")
    fig.savefig(FIG / "workflow.png", dpi=250, bbox_inches="tight")
    plt.close(fig)

    # Figure 2: frequency error by method/mode.
    labels = list(baseline["methods"].keys())
    x = np.arange(4)
    width = 0.24
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    for j, label in enumerate(labels):
        ax.bar(x + (j - 1) * width, baseline["methods"][label]["freq_error_pct"], width, label=label)
    ax.set_yscale("log")
    ax.set_xticks(x, [f"Mode {i}" for i in range(1, 5)])
    ax.set_ylabel("Absolute frequency error (%)")
    ax.legend(frameon=False)
    ax.grid(True, axis="y", which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG / "frequency_error.pdf", bbox_inches="tight")
    fig.savefig(FIG / "frequency_error.png", dpi=250, bbox_inches="tight")
    plt.close(fig)

    # Figure 3: MAC heatmap for expressive ansatz.
    MAC = baseline["methods"]["VQD L=2"]["MAC"]
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    im = ax.imshow(MAC, vmin=0, vmax=1, cmap="viridis")
    for i in range(4):
        for j in range(4):
            ax.text(j, i, f"{MAC[i, j]:.3f}", ha="center", va="center", color="white" if MAC[i,j] < 0.55 else "black")
    ax.set_xticks(range(4), [f"Q{i}" for i in range(1, 5)])
    ax.set_yticks(range(4), [f"C{i}" for i in range(1, 5)])
    ax.set_xlabel("Variational mode")
    ax.set_ylabel("Classical mode")
    fig.colorbar(im, ax=ax, label="Mass-weighted MAC")
    fig.tight_layout()
    fig.savefig(FIG / "mac_heatmap.pdf", bbox_inches="tight")
    fig.savefig(FIG / "mac_heatmap.png", dpi=250, bbox_inches="tight")
    plt.close(fig)

    # Figure 4: damage index.
    fig, ax = plt.subplots(figsize=(8.0, 4.5))
    x = np.arange(1, 7)
    ax.plot(x, DI_exact, marker="o", label="Classical modes")
    ax.plot(x, DI_q, marker="s", linestyle="--", label="Variational modes")
    ax.axvline(3, linestyle=":", linewidth=1.5, label="True damaged story")
    ax.set_xticks(x)
    ax.set_xlabel("Story spring")
    ax.set_ylabel("Four-mode strain-energy change index")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG / "damage_index.pdf", bbox_inches="tight")
    fig.savefig(FIG / "damage_index.png", dpi=250, bbox_inches="tight")
    plt.close(fig)

    # Figure 5: sensitivity heatmap.
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    im = ax.imshow(np.abs(S), aspect="auto", cmap="magma")
    ax.set_xticks(range(6), [str(i) for i in range(1, 7)])
    ax.set_yticks(range(4), [str(i) for i in range(1, 5)])
    ax.set_xlabel("Story spring")
    ax.set_ylabel("Mode")
    fig.colorbar(im, ax=ax, label=r"$|\partial\ln\lambda_i/\partial d_e|$")
    fig.tight_layout()
    fig.savefig(FIG / "sensitivity_heatmap.pdf", bbox_inches="tight")
    fig.savefig(FIG / "sensitivity_heatmap.png", dpi=250, bbox_inches="tight")
    plt.close(fig)

    # Figure 6: Pauli coefficient magnitude.
    fig, ax = plt.subplots(figsize=(8.2, 4.5))
    show = pauli.head(20).copy()
    ax.bar(np.arange(len(show)), np.abs(show["coefficient"]))
    ax.set_xticks(np.arange(len(show)), show["pauli"], rotation=60, ha="right")
    ax.set_ylabel("Absolute Pauli coefficient")
    ax.set_xlabel("Pauli string")
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG / "pauli_spectrum.pdf", bbox_inches="tight")
    fig.savefig(FIG / "pauli_spectrum.png", dpi=250, bbox_inches="tight")
    plt.close(fig)

    # Compact JSON summary for LaTeX/table generation.
    summary = {
        "exact_frequencies_hz": baseline["freq_exact"][:4].tolist(),
        "n_nonzero_pauli_terms": int(len(pauli)),
        "damage_story": 3,
        "damage_fraction": 0.10,
        "damage_index_peak_exact_story": int(np.argmax(DI_exact) + 1),
        "damage_index_peak_quantum_story": int(np.argmax(DI_q) + 1),
        "methods": {},
    }
    for label, d in baseline["methods"].items():
        summary["methods"][label] = {
            "frequency_error_pct": d["freq_error_pct"].tolist(),
            "residual": d["residual"].tolist(),
            "variance": d["result"].variances.tolist(),
            "padding_leakage": d["result"].padding_leakage.tolist(),
            "mac_diagonal": np.diag(d["MAC"]).tolist(),
        }
    (DATA / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
