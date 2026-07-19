"""Lightweight integrity checks for the manuscript's generated outputs."""
from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent


def close(actual: float, expected: float, rtol: float = 5e-4) -> None:
    if not np.isclose(actual, expected, rtol=rtol, atol=0.0):
        raise AssertionError(f"Expected {expected}, obtained {actual}")


def validate_citations() -> None:
    tex = (ROOT / "Qubits_Are_Not_Enough_TQE.tex").read_text(encoding="utf-8")
    bib = (ROOT / "references.bib").read_text(encoding="utf-8")
    cited: set[str] = set()
    for match in re.finditer(r"\\cite[pt]?\{([^}]+)\}", tex):
        cited.update(key.strip() for key in match.group(1).split(","))
    available = set(re.findall(r"@\w+\{([^,]+),", bib))
    missing = sorted(cited - available)
    if missing:
        raise AssertionError(f"Missing bibliography keys: {missing}")


def validate_numbers() -> None:
    scaling = pd.read_csv(ROOT / "revision_data" / "hamiltonian_scaling.csv")
    shot = pd.read_csv(ROOT / "revision_data" / "finite_shot_frequency.csv")
    damage = pd.read_csv(ROOT / "revision_data" / "damage_shot_localization.csv")
    floor = pd.read_csv(ROOT / "minor_revision_data" / "measurement_floor.csv")
    required = pd.read_csv(ROOT / "minor_revision_data" / "damage_shot_requirements.csv")

    consistent = scaling[(scaling.mass_type == "consistent") & (scaling.dof == 128)].iloc[0]
    lumped = scaling[(scaling.mass_type == "lumped") & (scaling.dof == 128)].iloc[0]
    assert int(consistent.pauli_terms) == 2950
    assert int(lumped.pauli_terms) == 576
    close(float(consistent.density_A), 0.3077392578125, rtol=1e-12)
    close(float(lumped.density_A), 0.0233154296875, rtol=1e-12)

    mode1 = shot[(shot["mode"] == 1) & (shot.shots_per_pauli == 100_000)].iloc[0]
    close(float(mode1.rmse_pct), 2.233719, rtol=2e-4)

    damage10 = damage[(damage.damage_fraction == 0.10) & (damage.shots_per_pauli == 100_000)].iloc[0]
    close(float(damage10.localization_accuracy), 0.26, rtol=1e-12)

    qwc = floor[(floor["mode"] == 1) & (floor.shots_per_term_independent == 100_000)].iloc[0]
    close(float(qwc.qwc_optimal_floor_rmse_pct), 0.8287785573, rtol=1e-9)

    req10 = required[required.damage_fraction == 0.10].iloc[0]
    close(float(req10.shots_per_term_detection_80pct_power), 4.820845582e5, rtol=1e-8)
    close(float(req10.shots_per_term_localization_80pct_power), 5.708457450e7, rtol=1e-8)

    summary = json.loads((ROOT / "minor_revision_data" / "minor_revision_summary.json").read_text())
    close(float(summary["qwc_floor_mode1_at_100k_pct"]), 0.8287785573, rtol=1e-9)


def validate_complexity() -> None:
    """Tie the manuscript complexity table (tab:complexity) and the illustrative
    128-DOF example (tab:concrete) to the reproduction code.

    Only deterministic quantities are asserted (qubit count, exact Pauli-term
    count, whitened density, and local scaling-exponent sanity bands). Wall-clock
    and Lanczos matvec counts are machine/library dependent and are deliberately
    not gated here.
    """
    sys.path.insert(0, str(ROOT))
    import revision_study as rs  # safe: only path setup + RNG seed run on import
    from scipy.linalg import eigh

    def whiten(n: int, mass_type: str) -> np.ndarray:
        M, K, _ = rs.chain_model(n, mass_type)
        A, _, _, _ = rs.mass_whiten(M, K)
        return A

    def pauli_terms(A: np.ndarray, n: int) -> int:
        lam_max = float(eigh(A, eigvals_only=True, subset_by_index=[n - 1, n - 1])[0])
        return int(len(rs.pauli_coefficients(A / lam_max, tol=1e-10)))

    def nnz_density(A: np.ndarray) -> tuple[int, float]:
        nnz = int(np.count_nonzero(np.abs(A) > 1e-12 * np.max(np.abs(A))))
        return nnz, nnz / A.size

    # Illustrative 128-DOF, k=4 consistent-mass example (tab:concrete).
    n = 128
    if int(round(math.log2(n))) != 7:
        raise AssertionError("N=128 must encode in 7 qubits")
    A128 = whiten(n, "consistent")
    n_pauli = pauli_terms(A128, n)
    if n_pauli != 2950:
        raise AssertionError(f"consistent-mass Pauli count expected 2950, got {n_pauli}")
    _, dens = nnz_density(A128)
    close(dens, 0.3077392578125, rtol=1e-9)

    # The recomputed values must match the committed scaling table.
    scaling = pd.read_csv(ROOT / "revision_data" / "hamiltonian_scaling.csv")
    stored = scaling[(scaling.mass_type == "consistent") & (scaling.dof == 128)].iloc[0]
    if int(stored.pauli_terms) != n_pauli:
        raise AssertionError("recomputed Pauli count disagrees with hamiltonian_scaling.csv")
    close(float(stored.density_A), dens, rtol=1e-12)

    # Illustrative circuit-execution estimate is order 1e11 (arithmetic only).
    execs = 1e5 * n_pauli * 300 * 4
    if not (1e11 <= execs <= 1e12):
        raise AssertionError(f"execution estimate outside 1e11-1e12 band: {execs:.2e}")

    # Asymptotic-row sanity: local scaling exponents over the tested range.
    ns = [16, 32, 64, 128]
    nnz_cons, pauli_cons, nnz_lump = [], [], []
    for nn in ns:
        Ac, Al = whiten(nn, "consistent"), whiten(nn, "lumped")
        nnz_cons.append(nnz_density(Ac)[0])
        nnz_lump.append(nnz_density(Al)[0])
        pauli_cons.append(pauli_terms(Ac, nn))

    def slope(ys: list) -> float:
        return float(np.polyfit(np.log(ns), np.log(ys), 1)[0])

    s_nnz, s_pauli, s_lump = slope(nnz_cons), slope(pauli_cons), slope(nnz_lump)
    # Consistent-mass densification is superlinear but within the O(N^2) upper bound.
    if not (1.0 < s_nnz < 2.05):
        raise AssertionError(f"nnz(A) consistent slope out of band: {s_nnz:.3f}")
    if not (1.0 < s_pauli < 2.05):
        raise AssertionError(f"Pauli consistent slope out of band: {s_pauli:.3f}")
    # Lumped mass stays near-linear (near-diagonal operator).
    if not (s_lump < 1.3):
        raise AssertionError(f"lumped nnz(A) slope unexpectedly high: {s_lump:.3f}")


def main() -> None:
    validate_citations()
    validate_numbers()
    validate_complexity()
    print("All citation, numerical, and complexity integrity checks passed.")


if __name__ == "__main__":
    main()
