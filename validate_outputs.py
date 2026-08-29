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

    # Stochastic checks. These use wide bands on purpose: each Pauli term draws from a
    # single shared Generator, so the result depends on the order of the term list, and
    # the six-DOF Hamiltonian has exactly tied |coefficient| groups whose sort order is
    # not guaranteed across environments. Permuting only those ties moves the values by
    # up to ~12%. See the reproducibility note in README.md. The deterministic checks
    # below and in validate_complexity() remain tight.
    mode1 = shot[(shot["mode"] == 1) & (shot.shots_per_pauli == 100_000)].iloc[0]
    close(float(mode1.rmse_pct), 2.447855, rtol=0.15)

    damage10 = damage[(damage.damage_fraction == 0.10) & (damage.shots_per_pauli == 100_000)].iloc[0]
    accuracy = float(damage10.localization_accuracy)
    if not 0.70 <= accuracy <= 0.95:
        raise AssertionError(
            f"10% damage localization at 1e5 shots/term expected in [0.70, 0.95], got {accuracy}"
        )

    qwc = floor[(floor["mode"] == 1) & (floor.shots_per_term_independent == 100_000)].iloc[0]
    close(float(qwc.qwc_optimal_floor_rmse_pct), 0.8287785573, rtol=1e-9)

    req10 = required[required.damage_fraction == 0.10].iloc[0]
    close(float(req10.shots_per_term_detection_80pct_power), 8.748548529e4, rtol=1e-8)
    close(float(req10.shots_per_term_localization_80pct_power), 4.826152558e5, rtol=1e-8)

    summary = json.loads((ROOT / "minor_revision_data" / "minor_revision_summary.json").read_text())
    close(float(summary["qwc_floor_mode1_at_100k_pct"]), 0.8287785573, rtol=1e-9)

    # Restart study (Sec. VII-D). Deterministic per seed, so these are tight. Depth 1
    # is bimodal between two ansatz-limited attractors; depth 2 reaches machine noise.
    restarts = pd.read_csv(ROOT / "revision_data" / "optimizer_restarts.csv")
    d1 = restarts[restarts.depth == 1].frequency_error_pct
    d2 = restarts[restarts.depth == 2].frequency_error_pct
    if len(d1) != 20 or len(d2) != 20:
        raise AssertionError("expected 20 restarts per depth")
    close(float(d1.min()), 127.9856553, rtol=1e-6)
    close(float(d1.median()), 130.1995103, rtol=1e-6)
    if not float(d2.median()) < 1e-10:
        raise AssertionError(f"depth-2 restart median not at machine noise: {d2.median()}")

    # Algorithm 1 executed with certified per-pair penalties (Appendix B).
    theory = pd.read_csv(ROOT / "revision_data" / "six_dof_theory_beta.csv")
    deep = theory[theory.depth == 2]
    if not (float(deep.err.max()) < 1e-12 and float(deep.mac.min()) > 1 - 1e-12
            and float(deep.overlap.max()) < 4e-15 and float(deep.leak.max()) < 1e-15):
        raise AssertionError("certified-penalty run fails the acceptance set quoted in Appendix B")

    # Derived quantities quoted in the text. Each previously existed only in prose; these
    # assertions are what keeps the manuscript and the generators from drifting apart.
    rev = json.loads((ROOT / "revision_data" / "revision_summary.json").read_text())
    dim = json.loads((ROOT / "dimension_generalization" / "data"
                      / "dimension_fit_statistics.json").read_text())
    quoted = [
        # (value, printed in the manuscript as, tolerance)
        (rev["pauli_exponent_full_range"], 1.66, 5e-3),          # Sec. VIII-B
        (rev["pauli_exponent_n_ge_16"], 1.48, 5e-3),             # Sec. VIII-B
        (rev["pauli_local_exponent_last_doubling"], 1.25, 5e-3),  # Sec. VI-D
        (dim["pauli_exponent_2D_full_range"], 1.69, 5e-3),       # Sec. VIII-B
        (dim["pauli_exponent_3D_full_range"], 1.79, 5e-3),       # Sec. VIII-B
        (rev["gradient_normalized_slope"], -0.036, 5e-4),        # Sec. VII-E
        (rev["gradient_normalized_ci_lo"], -0.087, 5e-4),
        (rev["gradient_normalized_ci_hi"], 0.016, 5e-4),
        (rev["gradient_normalized_r2"], 0.17, 5e-3),
        (rev["gradient_tail_slope"], -0.118, 5e-4),              # Sec. VII-E
        (rev["certified_beta_min_depth2"], 0.68, 5e-3),          # Appendix B
        (rev["certified_beta_max_depth2"], 0.98, 5e-3),
        (rev["certified_beta_max_depth1"], 1.02, 5e-3),
        (rev["wilson_10pct_100000_lo"], 78.8, 0.05),             # Sec. VII-G
        (rev["wilson_10pct_100000_hi"], 88.9, 0.05),
        (rev["wilson_10pct_10000_lo"], 33.9, 0.05),
        (rev["wilson_10pct_10000_hi"], 47.4, 0.05),
        (rev["depth_sweep_slope_L2"], -0.42, 5e-3),               # Sec. VII-E
        (rev["depth_sweep_slope_L8"], -0.73, 5e-3),
    ]
    # The headline of the depth sweep: decay steepens with depth and reaches the
    # barren-plateau rate by L=8. Assert the ordering, not just the values.
    if not (rev["depth_sweep_slope_L8"] < rev["depth_sweep_slope_L2"] < 0):
        raise AssertionError("depth sweep no longer shows decay steepening with depth")
    if rev["depth_sweep_slope_L8"] > -math.log(2):
        raise AssertionError(
            f"L=8 slope {rev['depth_sweep_slope_L8']:.3f} no longer reaches the "
            f"barren-plateau rate {-math.log(2):.3f}"
        )
    for got, printed, tol in quoted:
        if abs(float(got) - printed) > tol:
            raise AssertionError(
                f"generated {got} does not round to the manuscript's {printed} (tol {tol})"
            )


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


def validate_regeneration(full: bool = False) -> None:
    """Re-execute the generators and diff against the committed CSVs.

    This is the only check here that can detect a *code* regression.
    ``validate_numbers`` compares committed files against literals baked into this
    script, so it verifies file integrity and nothing else -- ``MANIFEST.sha256``
    already does that better. The studies below are deterministic and reproduce
    bit-for-bit, so they are compared with ``rtol=0``.

    Generators write through the module-level ``DATA``/``FIG`` paths, which are
    redirected to a temporary directory here so the committed package is never
    touched. ``full=True`` adds the two Monte Carlo studies (~2 min); they are
    deterministic too, but only when the shared generator is at the state
    ``main()`` leaves it in, which is why the order below is fixed.
    """
    import shutil
    import tempfile

    sys.path.insert(0, str(ROOT))
    import matplotlib
    matplotlib.use("Agg")
    import numpy as _np
    import revision_study as rs
    import minor_revision_study as mrs

    def diff(name: str, produced: Path, committed: Path) -> None:
        a = pd.read_csv(produced / name)
        b = pd.read_csv(committed / name)
        if list(a.columns) != list(b.columns):
            raise AssertionError(f"{name}: column mismatch")
        for c in a.columns:
            if c.endswith("_seconds"):        # wall clock is machine-specific
                continue
            if not a[c].equals(b[c]):
                bad = int((a[c] != b[c]).sum())
                raise AssertionError(
                    f"{name}: column '{c}' differs from committed data in {bad} row(s); "
                    f"first produced={a[c].iloc[0]!r} committed={b[c].iloc[0]!r}"
                )

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        rev, fig, minor = tmp / "rev", tmp / "fig", tmp / "minor"
        for d in (rev, fig, minor):
            d.mkdir()
        saved = (rs.DATA, rs.FIG, mrs.DATA)
        rs.DATA, rs.FIG, mrs.DATA = rev, fig, minor
        try:
            rs.theorem_penalty_study()
            diff("penalty_threshold.csv", rev, ROOT / "revision_data")
            diff("deflation_perturbation.csv", rev, ROOT / "revision_data")

            rs.scaling_and_pauli_study()
            diff("hamiltonian_scaling.csv", rev, ROOT / "revision_data")

            rs.restart_study()
            diff("optimizer_restarts.csv", rev, ROOT / "revision_data")

            rs.theory_beta_study()
            diff("six_dof_theory_beta.csv", rev, ROOT / "revision_data")

            mrs.measurement_floor_study()
            diff("measurement_floor.csv", minor, ROOT / "minor_revision_data")

            mrs.damage_shot_requirement_study("baseline")
            diff("damage_shot_requirements.csv", minor, ROOT / "minor_revision_data")
            mrs.damage_shot_requirement_study("oracle")
            diff("damage_shot_requirements_oracle_convention.csv", minor,
                 ROOT / "minor_revision_data")

            if full:
                rs.trainability_depth_sweep()          # own generator, ~35 s
                diff("trainability_depth_sweep.csv", rev, ROOT / "revision_data")
                # Fixed order: finite_shot_frequency is main()'s first consumer of the
                # shared generator and damage_shot_study its second.
                rs.RNG = _np.random.default_rng(20260623)
                rs.finite_shot_frequency()
                diff("finite_shot_frequency.csv", rev, ROOT / "revision_data")
                rs.damage_shot_study()
                diff("damage_shot_localization.csv", rev, ROOT / "revision_data")
        finally:
            rs.DATA, rs.FIG, mrs.DATA = saved


def main() -> None:
    full = "--full" in sys.argv
    validate_citations()
    validate_numbers()
    validate_complexity()
    validate_regeneration(full=full)
    scope = "including Monte Carlo" if full else "deterministic only; use --full for Monte Carlo"
    print(f"All citation, numerical, complexity, and regeneration checks passed ({scope}).")


if __name__ == "__main__":
    main()
