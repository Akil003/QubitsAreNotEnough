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

    # sec:genquotient -- generalized Rayleigh quotient representation cost.
    gq = pd.read_csv(ROOT / "revision_data" / "generalized_quotient_comparison.csv")
    gq_l = gq[gq.mass_type == "lumped"].set_index("dof")
    gq_c = gq[gq.mass_type == "consistent"].set_index("dof")
    # Cross-check: the A columns must agree with hamiltonian_scaling / tab:scaling.
    for dof in (16, 64, 128):
        assert int(gq_l.loc[dof, "pauli_A"]) == int(
            scaling[(scaling.mass_type == "lumped") & (scaling.dof == dof)].pauli_terms.iloc[0])
        assert int(gq_c.loc[dof, "pauli_A"]) == int(
            scaling[(scaling.mass_type == "consistent") & (scaling.dof == dof)].pauli_terms.iloc[0])
    # Values quoted in the text and in tab:genquotient.
    assert int(gq_l.loc[128, "pauli_KM"]) == 503 and int(gq_l.loc[128, "pauli_A"]) == 576
    assert int(gq_c.loc[128, "pauli_KM"]) == 750 and int(gq_c.loc[128, "pauli_A"]) == 2950
    close(float(gq_l.loc[128, "ratio_KM_over_A"]), 503 / 576, rtol=1e-12)
    close(float(gq_c.loc[128, "ratio_KM_over_A"]), 750 / 2950, rtol=1e-12)
    assert int(gq_c.loc[64, "qwc_KM"]) == 64 and int(gq_c.loc[64, "qwc_A"]) == 221
    # The qualitative claim: under lumped mass the ratio crosses one only above N=32.
    assert (gq_l.loc[[4, 8, 16, 32], "ratio_KM_over_A"] > 1).all()
    assert (gq_l.loc[[64, 128], "ratio_KM_over_A"] < 1).all()
    # The claim that M opens no group beyond those K already requires.
    for frame in (gq_l, gq_c):
        sub = frame.loc[[4, 8, 16, 32, 64]]
        assert (sub["qwc_KM"] == sub["qwc_K"]).all()
    # sec:genquotient claims prop:lumped's counting bound covers the unwhitened pair too,
    # since K and M are real symmetric with nearest-neighbour couplings in the same indexing.
    for _, row in gq.iterrows():
        bound = (int(row.qubits) + 2) * 2 ** (int(row.qubits) - 1)
        assert int(row.pauli_K) <= bound, (row.dof, row.pauli_K, bound)
        assert int(row.pauli_M) <= bound, (row.dof, row.pauli_M, bound)
    # Lumped settings are exactly neutral; consistent settings strictly favour (K, M).
    grouped = [4, 8, 16, 32, 64]
    assert (gq_l.loc[grouped, "qwc_KM"] == gq_l.loc[grouped, "qwc_A"]).all()
    assert (gq_c.loc[[8, 16, 32, 64], "qwc_KM"] < gq_c.loc[[8, 16, 32, 64], "qwc_A"]).all()

    # G8: every reported localization proportion carries a Wilson interval, and each
    # interval quoted in the text must match the persisted one.
    rev_sum = json.loads((ROOT / "revision_data" / "revision_summary.json").read_text())
    for key, quoted in (
        ("wilson_10pct_100000", (78.8, 88.9)),      # 84.5%
        ("wilson_10pct_10000", (33.9, 47.4)),       # 40.5%
        ("wilson_10pct_1000", (17.7, 29.3)),        # 23.0%
        ("wilson_2pct_100000", (25.9, 38.8)),       # 32.0%
        ("wilson_oracle_10pct_100000", (19.5, 31.4)),  # 25.0%, Remark 3
    ):
        close(rev_sum[f"{key}_lo"], quoted[0], rtol=6e-3)
        close(rev_sum[f"{key}_hi"], quoted[1], rtol=6e-3)
    # The oracle convention must still give 25.0% at 1e5 -- the value Remark 3 quotes.
    oracle = pd.read_csv(ROOT / "revision_data"
                         / "damage_shot_localization_oracle_convention.csv")
    orow = oracle[(oracle.damage_fraction == 0.10) & (oracle.shots_per_pauli == 100_000)]
    close(float(orow.localization_accuracy.iloc[0]), 0.25, rtol=1e-12)

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
    # B3: the manuscript asserts the positivity clip and the denominator floor are
    # inactive at the two highest budgets. That is what makes those rows quotable.
    for _, r in shot.iterrows():
        if r.shots_per_pauli >= 10_000 and float(r.positivity_clip_fraction) != 0.0:
            raise AssertionError(
                f"positivity clip active at {int(r.shots_per_pauli)} shots, mode {int(r['mode'])}"
            )
    if abs(float(shot[(shot.shots_per_pauli == 100) & (shot["mode"] == 1)]
                 .positivity_clip_fraction.iloc[0]) - 0.238) > 0.002:
        raise AssertionError("mode-1 clip fraction at 1e2 shots no longer matches the quoted 23.8%")
    for _, r in damage.iterrows():
        if r.shots_per_pauli >= 10_000 and float(r.denominator_clip_fraction) != 0.0:
            raise AssertionError(
                f"denominator floor active at {int(r.shots_per_pauli)} shots, "
                f"{r.damage_fraction:.0%} damage"
            )

    # B2: bootstrap intervals on the depth-sweep slopes, as quoted in Sec. VII-E.
    ci = pd.read_csv(ROOT / "revision_data" / "trainability_depth_sweep_ci.csv")
    for depth, lo_q, hi_q in ((2, -0.49, -0.36), (4, -0.67, -0.57), (8, -0.77, -0.69)):
        r = ci[(ci.depth == depth) & (ci.range == "all")].iloc[0]
        if abs(float(r.ci_lo) - lo_q) > 0.005 or abs(float(r.ci_hi) - hi_q) > 0.005:
            raise AssertionError(
                f"L={depth} CI [{r.ci_lo:.3f}, {r.ci_hi:.3f}] does not match the quoted "
                f"[{lo_q}, {hi_q}]"
            )
    # The qualitative claim: L=2 excludes the O(2^-n_q) reference, L=8 does not.
    ref = -math.log(2)
    l2 = ci[(ci.depth == 2) & (ci.range == "all")].iloc[0]
    l8 = ci[(ci.depth == 8) & (ci.range == "all")].iloc[0]
    if not float(l2.ci_lo) > ref:
        raise AssertionError("L=2 interval no longer excludes the O(2^-n_q) reference")
    if not (float(l8.ci_lo) <= ref <= float(l8.ci_hi)):
        raise AssertionError("L=8 interval no longer contains the O(2^-n_q) reference")

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


def validate_lumped_pauli_count() -> None:
    """Check Proposition 4's equality condition over ALL admissible pairs.

    The count (n_q+2)2^(n_q-1) decomposes as 2^n_q from the x=0 diagonal sector plus
    n_q*2^(n_q-1) from the coupling masks. The diagonal sector is not automatic: a
    constant diagonal kills 2^n_q - 1 of it. This check therefore counts survivors in
    both sectors separately, and confirms the heterogeneous meshes attain the bound
    while the uniform ones fall to 3*2^(n_q-1) - 1.
    """
    sys.path.insert(0, str(ROOT))
    import revision_study as rs
    from scipy.linalg import eigh

    def survivors(A: np.ndarray) -> tuple[int, int]:
        d = A.shape[0]
        basis = np.arange(d, dtype=np.int64)
        thr = 1e-12 * np.abs(A).max()
        masks = sorted({i ^ j for i in range(d) for j in range(d) if abs(A[i, j]) > thr})
        diag = off = 0
        for x in masks:
            v = A[basis, basis ^ x]
            scale = 1e-10 * max(1.0, float(np.abs(v).max()))
            for z in range(d):
                if bin(x & z).count("1") % 2:
                    continue
                par = np.array([bin(int(i) & z).count("1") % 2 for i in basis])
                if abs(float(np.sum(((-1.0) ** par) * v))) > scale:
                    if x == 0:
                        diag += 1
                    else:
                        off += 1
        return diag, off

    def uniform_chain(n: int) -> tuple[np.ndarray, np.ndarray]:
        K = np.zeros((n + 1, n + 1))
        M = np.zeros((n + 1, n + 1))
        for e in range(n):
            ke = 2.2e8 * np.array([[1.0, -1.0], [-1.0, 1.0]])
            for a in range(2):
                for b in range(2):
                    K[e + a, e + b] += ke[a, b]
                    M[e + a, e + b] += 1.0e4 * (1.0 if a == b else 0.0)
        return M[1:, 1:], K[1:, 1:]

    for nq in (3, 4, 5):
        n = 2 ** nq
        cases = {"heterogeneous": rs.chain_model(n, "lumped")[:2], "uniform": uniform_chain(n)}
        for label, (M, K) in cases.items():
            A, _, _, _ = rs.mass_whiten(M, K)
            H = A / float(eigh(A, eigvals_only=True)[-1])
            diag, off = survivors(H)
            reported = len(rs.pauli_coefficients(H, tol=1e-10))
            if diag + off != reported:
                raise AssertionError(
                    f"{label} n_q={nq}: admissible-pair count {diag}+{off} != N_P {reported}"
                )
            if label == "heterogeneous":
                if diag != 2 ** nq:
                    raise AssertionError(
                        f"heterogeneous n_q={nq}: x=0 sector {diag}, expected full {2**nq}"
                    )
                if reported != (nq + 2) * 2 ** (nq - 1):
                    raise AssertionError(f"heterogeneous n_q={nq}: bound not attained")
            else:
                if diag != 1:
                    raise AssertionError(
                        f"uniform n_q={nq}: x=0 sector {diag}, expected 1 (constant diagonal)"
                    )
                if reported != 3 * 2 ** (nq - 1) - 1:
                    raise AssertionError(
                        f"uniform n_q={nq}: expected {3*2**(nq-1)-1}, got {reported}"
                    )


def validate_ordering() -> None:
    """sec:ordering -- fill-reducing orderings trade one sparsity for the other."""
    df = pd.read_csv(ROOT / "dimension_generalization" / "data" / "ordering_study.csv")
    assert len(df) == 12, len(df)
    idx = df.set_index(["dimension", "mass_type", "ordering"])

    # The natural ordering attains prop:lumped exactly, and is the minimum everywhere.
    for dim in ("2D", "3D"):
        nat = idx.loc[(dim, "lumped", "natural")]
        assert int(nat.pauli_terms) == int(nat.pauli_bound_lumped), (dim, nat.pauli_terms)
        for mass in ("lumped", "consistent"):
            grp = df[(df.dimension == dim) & (df.mass_type == mass)]
            assert int(grp.pauli_terms.min()) == int(
                grp[grp.ordering == "natural"].pauli_terms.iloc[0]), (dim, mass)

    # A random scramble saturates the real-symmetric maximum exactly (consistent mass).
    for dim in ("2D", "3D"):
        rnd = idx.loc[(dim, "consistent", "random")]
        assert int(rnd.pauli_terms) == int(rnd.pauli_bound_realsym), (dim, rnd.pauli_terms)

    # Lumped nnz(L) and density_A are permutation-invariant; only N_P moves.
    for dim in ("2D", "3D"):
        lump = df[(df.dimension == dim) & (df.mass_type == "lumped")]
        assert lump.nnz_L.nunique() == 1 and int(lump.nnz_L.iloc[0]) == int(lump.dof.iloc[0])
        assert lump.density_A.nunique() == 1, lump.density_A.tolist()
        assert lump.pauli_terms.nunique() == 3, lump.pauli_terms.tolist()

    # RCM reduces consistent-mass fill but raises N_P -- the trade, in both dimensions.
    for dim, fill_drop, np_rise in (("2D", 0.28, 1.30), ("3D", 0.33, 1.15)):
        nat = idx.loc[(dim, "consistent", "natural")]
        rcm = idx.loc[(dim, "consistent", "rcm")]
        assert rcm.nnz_L < nat.nnz_L * (1 - fill_drop), (dim, rcm.nnz_L, nat.nnz_L)
        assert rcm.pauli_terms > nat.pauli_terms * np_rise, (dim, rcm.pauli_terms)

    # Quoted ratios.
    close(float(idx.loc[("2D", "lumped", "rcm")].pauli_ratio_to_natural), 5.60, rtol=2e-3)
    close(float(idx.loc[("3D", "lumped", "rcm")].pauli_ratio_to_natural), 14.14, rtol=2e-3)
    close(float(idx.loc[("3D", "lumped", "random")].pauli_ratio_to_natural), 43.73, rtol=2e-3)
    assert int(idx.loc[("2D", "consistent", "natural")].pauli_terms) == 25082
    assert int(idx.loc[("3D", "consistent", "rcm")].nnz_L) == 19566


def validate_near_degenerate() -> None:
    """sec:neardegenerate -- the near-degenerate sweep and its two certificates."""
    df = pd.read_csv(ROOT / "dimension_generalization" / "data" / "near_degenerate_sweep.csv")
    assert len(df) == 72, len(df)

    # Both Davis-Kahan bounds must actually hold everywhere; this is the check that
    # would catch a sign or indexing error in the sweep.
    assert (df.sin_individual <= df.dk_individual + 1e-12).all()
    assert (df.sin_subspace <= df.dk_subspace + 1e-12).all()

    # Weyl's precondition and non-vacuity are the same inequality (claimed in the text).
    assert (df.weyl_ok_individual == (df.dk_individual < 1)).all()
    assert (df.weyl_ok_subspace == (df.dk_subspace < 1)).all()

    can = df[df.canonical].set_index("delta")
    assert len(can) == 6, len(can)
    # Quoted in tab:neardeg.
    close(float(can.loc[0.001, "mac_r"]), 0.026472, rtol=1e-3)
    close(float(can.loc[0.002, "mac_r"]), 0.147070, rtol=1e-3)
    close(float(can.loc[0.050, "mac_r"]), 0.999974, rtol=1e-5)
    close(float(can.loc[0.001, "dk_individual"]), 18.697146, rtol=1e-4)
    close(float(can.loc[0.050, "dk_individual"]), 0.392308, rtol=1e-4)
    # The per-mode diagnostic degrades monotonically as the gap closes.
    assert can.sort_index()["mac_r"].is_monotonic_increasing
    # eps_sub flat to three significant figures across the canonical sweep.
    assert can["eps_sub"].max() / can["eps_sub"].min() < 1.005, can["eps_sub"].describe()
    # The individual certificate lapses between delta = 1% and 2%.
    assert (can.loc[[0.001, 0.002, 0.005, 0.010], "dk_individual"] >= 1).all()
    assert (can.loc[[0.020, 0.050], "dk_individual"] < 1).all()
    # Perturbation held essentially constant along the canonical sweep, as the caption
    # states. Other probe edges give a different (but likewise constant) magnitude.
    assert 2.10e-3 <= can.norm_delta.min() and can.norm_delta.max() <= 2.16e-3

    # Robustness across all 12 probe locations, quoted in the text.
    lo, hi = df[df.delta == 0.001], df[df.delta == 0.05]
    assert (lo.dk_individual >= 1).all() and (hi.dk_individual < 1).all()
    close(float(lo.dk_individual.min()), 8.7625, rtol=1e-3)
    close(float(lo.dk_individual.max()), 27.667, rtol=1e-3)
    close(float(hi.dk_individual.min()), 0.17528, rtol=1e-3)
    close(float(hi.dk_individual.max()), 0.58052, rtol=1e-3)
    # The block certificate stays non-vacuous everywhere, by over an order of magnitude.
    assert df.dk_subspace.max() < 0.054, float(df.dk_subspace.max())
    assert 0.0021 <= df.eps_sub.min() and df.eps_sub.max() <= 0.0043
    assert (df.g_pair > 0.11).all()


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
    sys.path.insert(0, str(ROOT / "dimension_generalization"))
    import mesh_dimension_study as mds

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
        rev, fig, minor, dim = tmp / "rev", tmp / "fig", tmp / "minor", tmp / "dim"
        for d in (rev, fig, minor, dim):
            d.mkdir()
        saved = (rs.DATA, rs.FIG, mrs.DATA, mds.DATA)
        rs.DATA, rs.FIG, mrs.DATA, mds.DATA = rev, fig, minor, dim
        try:
            rs.theorem_penalty_study()
            diff("penalty_threshold.csv", rev, ROOT / "revision_data")
            diff("deflation_perturbation.csv", rev, ROOT / "revision_data")

            rs.scaling_and_pauli_study()
            diff("hamiltonian_scaling.csv", rev, ROOT / "revision_data")

            rs.generalized_quotient_comparison()
            diff("generalized_quotient_comparison.csv", rev, ROOT / "revision_data")

            mds.near_degenerate_sweep()
            diff("near_degenerate_sweep.csv", dim,
                 ROOT / "dimension_generalization" / "data")

            if full:                      # ~30 s: four n_q=8/9 Pauli decompositions
                mds.ordering_study()
                diff("ordering_study.csv", dim,
                     ROOT / "dimension_generalization" / "data")

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
            rs.DATA, rs.FIG, mrs.DATA, mds.DATA = saved


def main() -> None:
    full = "--full" in sys.argv
    validate_citations()
    validate_numbers()
    validate_complexity()
    validate_lumped_pauli_count()
    validate_near_degenerate()
    validate_ordering()
    validate_regeneration(full=full)
    scope = "including Monte Carlo" if full else "deterministic only; use --full for Monte Carlo"
    print(f"All citation, numerical, complexity, and regeneration checks passed ({scope}).")


if __name__ == "__main__":
    main()
