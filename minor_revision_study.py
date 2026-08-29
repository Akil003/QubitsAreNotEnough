"""Minor-revision calculations requested by the second reviewer.

This script adds two reproducible analyses to the major-revision study:

1. An optimistic qubit-wise-commuting (QWC) measurement floor for modal
   frequencies at the same total state-preparation budget as the independent
   Pauli baseline. Exact within-group covariances and continuous optimal shot
   allocation are used, so the result is a lower bound rather than a hardware
   forecast.
2. Delta-method shot requirements for detecting and localizing modal-strain-
   energy changes. The calculation reports both a single-element detection
   requirement and a more demanding pairwise localization requirement.

The model, Hamiltonian mapping, and Pauli utilities are imported from
``revision_study.py`` so the new calculations use exactly the same six-DOF
benchmark as the manuscript.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

import revision_study as rs

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "minor_revision_data"
DATA.mkdir(exist_ok=True)


def _is_identity(pauli: str) -> bool:
    return set(pauli) == {"I"}


def _qwc_compatible(pauli_a: str, pauli_b: str) -> bool:
    """Return True when two Pauli strings commute qubit-wise."""
    return all(a == "I" or b == "I" or a == b for a, b in zip(pauli_a, pauli_b))


def _greedy_qwc_groups(pauli_df: pd.DataFrame) -> list[list[int]]:
    """Return greedy QWC groups as dataframe-index lists.

    Identity terms are excluded because they require no measurements.
    """
    groups: list[list[int]] = []
    for idx, row in pauli_df.iterrows():
        pauli = str(row.pauli)
        if _is_identity(pauli):
            continue
        for group in groups:
            if all(_qwc_compatible(pauli, str(pauli_df.loc[j].pauli)) for j in group):
                group.append(int(idx))
                break
        else:
            groups.append([int(idx)])
    return groups


def _qwc_product(pauli_a: str, pauli_b: str) -> str:
    """Product of qubit-wise-compatible Pauli strings.

    QWC strings never contain different non-identity operators on the same
    qubit, so the product has no imaginary phase.
    """
    product: list[str] = []
    for a, b in zip(pauli_a, pauli_b):
        if a == "I":
            product.append(b)
        elif b == "I":
            product.append(a)
        elif a == b:
            product.append("I")
        else:
            raise ValueError(f"Pauli strings are not QWC: {pauli_a}, {pauli_b}")
    return "".join(product)


def operator_mean_variance_coefficient(
    pauli_df: pd.DataFrame, state: np.ndarray
) -> tuple[float, float, int]:
    """Return exact mean, equal-shots variance coefficient, and term count.

    If each non-identity Pauli term receives ``n`` shots independently, the
    estimator variance is ``variance_coefficient / n``.
    """
    mean = 0.0
    variance_coefficient = 0.0
    nonidentity_terms = 0
    for row in pauli_df.itertuples(index=False):
        coefficient = float(row.coefficient)
        pauli = str(row.pauli)
        identity = _is_identity(pauli)
        mu = 1.0 if identity else rs.pauli_expectation(state, pauli)
        mean += coefficient * mu
        if not identity:
            variance_coefficient += coefficient**2 * (1.0 - mu**2)
            nonidentity_terms += 1
    return float(mean), float(variance_coefficient), int(nonidentity_terms)


def qwc_optimal_variance_coefficient(
    pauli_df: pd.DataFrame, state: np.ndarray
) -> tuple[float, int]:
    """Return the optimal-QWC coefficient for a fixed total shot budget.

    For QWC group ``g``, let ``v_g`` be the exact single-shot variance of the
    weighted group observable. With continuous optimal allocation and total
    budget ``B``, the minimum variance is

        (sum_g sqrt(v_g))**2 / B.

    This assumes perfect simultaneous measurement, exact covariance knowledge,
    no integer-allocation loss, and no circuit or readout noise.
    """
    groups = _greedy_qwc_groups(pauli_df)
    mu: dict[int, float] = {}
    for idx, row in pauli_df.iterrows():
        pauli = str(row.pauli)
        mu[int(idx)] = 1.0 if _is_identity(pauli) else rs.pauli_expectation(state, pauli)

    group_variances: list[float] = []
    for group in groups:
        variance = 0.0
        for i in group:
            row_i = pauli_df.loc[i]
            h_i = float(row_i.coefficient)
            p_i = str(row_i.pauli)
            for j in group:
                row_j = pauli_df.loc[j]
                h_j = float(row_j.coefficient)
                p_j = str(row_j.pauli)
                product = _qwc_product(p_i, p_j)
                mu_product = 1.0 if _is_identity(product) else rs.pauli_expectation(state, product)
                variance += h_i * h_j * (mu_product - mu[i] * mu[j])
        # Numerical roundoff can produce tiny negative values.
        group_variances.append(max(float(variance), 0.0))

    coefficient = float(sum(math.sqrt(v) for v in group_variances) ** 2)
    return coefficient, len(groups)


def measurement_floor_study() -> pd.DataFrame:
    """Compute independent and optimistic grouped frequency RMSE estimates."""
    _, _, _, _, _, eigenvalues, modes, scale, hamiltonian = rs.exact_modal_data()
    pauli_df = rs.pauli_coefficients(hamiltonian, tol=1e-12)
    nonidentity_terms = sum(not _is_identity(str(p)) for p in pauli_df.pauli)

    rows: list[dict[str, float | int]] = []
    for mode_index in range(4):
        state = modes[:, mode_index]
        _, independent_coeff, _ = operator_mean_variance_coefficient(pauli_df, state)
        qwc_coeff, group_count = qwc_optimal_variance_coefficient(pauli_df, state)
        frequency = math.sqrt(float(eigenvalues[mode_index])) / (2.0 * math.pi)
        derivative = scale / (8.0 * math.pi**2 * frequency)

        for shots_per_term in (100, 1_000, 10_000, 100_000):
            total_budget = shots_per_term * nonidentity_terms
            independent_rmse_pct = (
                100.0
                * derivative
                * math.sqrt(independent_coeff / shots_per_term)
                / frequency
            )
            qwc_floor_rmse_pct = (
                100.0
                * derivative
                * math.sqrt(qwc_coeff / total_budget)
                / frequency
            )
            rows.append(
                {
                    "mode": mode_index + 1,
                    "frequency_hz": frequency,
                    "shots_per_term_independent": shots_per_term,
                    "matched_total_shots": total_budget,
                    "nonidentity_pauli_terms": nonidentity_terms,
                    "qwc_groups": group_count,
                    "independent_analytic_rmse_pct": independent_rmse_pct,
                    "qwc_optimal_floor_rmse_pct": qwc_floor_rmse_pct,
                    "rmse_improvement_factor": independent_rmse_pct / qwc_floor_rmse_pct,
                }
            )

    dataframe = pd.DataFrame(rows)
    dataframe.to_csv(DATA / "measurement_floor.csv", index=False)
    return dataframe


def _padded_element_operators(
    mass: np.ndarray,
    element_stiffness: list[np.ndarray],
    scale: float,
    dimension: int,
) -> list[np.ndarray]:
    operators: list[np.ndarray] = []
    for stiffness in element_stiffness:
        whitened, _, _, _ = rs.mass_whiten(mass, stiffness)
        operators.append(rs.pad(whitened / scale, dimension, 0.0))
    return operators


def damage_shot_requirement_study(convention: str = "baseline") -> pd.DataFrame:
    """Estimate shots needed for modal-energy detection and localization.

    ``convention`` selects the element stiffnesses used to build the numerator
    operators on the damaged structure: ``"baseline"`` (Ke0, the deployable and
    reported case) or ``"oracle"`` (Ke^d, retained only to reproduce the manuscript's
    Remark 3 sensitivity case). This calculation is fully deterministic, so both
    variants reproduce exactly.

    Assumptions
    -----------
    * First four exact modes of the six-DOF benchmark are used.
    * Element-energy numerators on both the baseline and the damaged structure are
      built from the BASELINE element stiffnesses, following the standard modal-
      strain-energy convention; the damaged element stiffness is the unknown and
      cannot be used to assemble a measurable operator.
    * Every non-identity Pauli term receives the same number ``n`` of shots.
    * Numerator and denominator operators, baseline and damaged states, and
      different element scores are treated as independent. This is conservative
      for shared-denominator implementations because beneficial covariance is
      ignored.
    * The absolute-value signs are linearized at their exact signs.
    * Detection uses a two-sided 5% test with 80% power.
    * Localization uses a one-sided Bonferroni family-wise 5% comparison against
      the five intact competitors, also at 80% power.
    """
    base = rs.exact_modal_data()
    mass_0, _, stiffness_0, _, _, _, modes_0, scale_0, hamiltonian_0 = base
    element_0 = _padded_element_operators(
        mass_0, stiffness_0, scale_0, hamiltonian_0.shape[0]
    )
    pauli_h0 = rs.pauli_coefficients(hamiltonian_0, tol=1e-12)
    pauli_e0 = [rs.pauli_coefficients(operator, tol=1e-12) for operator in element_0]

    z_detection = norm.ppf(0.975) + norm.ppf(0.80)
    z_localization = norm.ppf(1.0 - 0.05 / 5.0) + norm.ppf(0.80)

    rows: list[dict[str, float | int]] = []
    for severity in (0.02, 0.05, 0.10, 0.20):
        damaged = rs.exact_modal_data(damage_story=2, damage_fraction=severity)
        mass_d, _, stiffness_d, _, _, _, modes_d, scale_d, hamiltonian_d = damaged
        # Baseline element stiffnesses, not the damaged ones, build the operators
        # measured on the damaged structure -- see the note in
        # revision_study.damage_shot_study.
        stiffness_numerator = stiffness_0 if convention == "baseline" else stiffness_d
        element_d = _padded_element_operators(
            mass_d, stiffness_numerator, scale_d, hamiltonian_d.shape[0]
        )
        pauli_hd = rs.pauli_coefficients(hamiltonian_d, tol=1e-12)
        pauli_ed = [rs.pauli_coefficients(operator, tol=1e-12) for operator in element_d]

        scores = np.zeros(6)
        variance_coefficients = np.zeros(6)
        measurement_settings = 0

        for mode_index in range(4):
            b0, vb0, nb0 = operator_mean_variance_coefficient(
                pauli_h0, modes_0[:, mode_index]
            )
            bd, vbd, nbd = operator_mean_variance_coefficient(
                pauli_hd, modes_d[:, mode_index]
            )
            measurement_settings += nb0 + nbd

            for element_index in range(6):
                a0, va0, na0 = operator_mean_variance_coefficient(
                    pauli_e0[element_index], modes_0[:, mode_index]
                )
                ad, vad, nad = operator_mean_variance_coefficient(
                    pauli_ed[element_index], modes_d[:, mode_index]
                )
                measurement_settings += na0 + nad

                eta_0 = a0 / b0
                eta_d = ad / bd
                difference = eta_d - eta_0
                ratio_difference_variance_coefficient = (
                    va0 / b0**2
                    + (a0**2) * vb0 / b0**4
                    + vad / bd**2
                    + (ad**2) * vbd / bd**4
                )
                weight = 0.25
                scores[element_index] += weight * abs(difference)
                variance_coefficients[element_index] += (
                    weight**2 * ratio_difference_variance_coefficient
                )

        damaged_element = 2
        competitors = [index for index in range(6) if index != damaged_element]
        strongest_competitor = max(competitors, key=lambda index: scores[index])
        localization_margin = scores[damaged_element] - scores[strongest_competitor]
        localization_variance_coefficient = (
            variance_coefficients[damaged_element]
            + variance_coefficients[strongest_competitor]
        )

        detection_shots = (
            z_detection**2
            * variance_coefficients[damaged_element]
            / scores[damaged_element] ** 2
        )
        localization_shots = (
            z_localization**2
            * localization_variance_coefficient
            / localization_margin**2
        )
        snr_at_100k = localization_margin * math.sqrt(
            100_000.0 / localization_variance_coefficient
        )

        rows.append(
            {
                "damage_fraction": severity,
                "damaged_element_score": scores[damaged_element],
                "strongest_competitor_element": strongest_competitor + 1,
                "strongest_competitor_score": scores[strongest_competitor],
                "exact_localization_margin": localization_margin,
                "target_variance_coefficient": variance_coefficients[damaged_element],
                "margin_variance_coefficient": localization_variance_coefficient,
                "shots_per_term_detection_80pct_power": detection_shots,
                "shots_per_term_localization_80pct_power": localization_shots,
                "nonidentity_measurement_settings_per_trial": measurement_settings,
                "total_shots_detection": detection_shots * measurement_settings,
                "total_shots_localization": localization_shots * measurement_settings,
                "localization_snr_at_100k_shots_per_term": snr_at_100k,
            }
        )

    dataframe = pd.DataFrame(rows)
    suffix = "" if convention == "baseline" else "_oracle_convention"
    dataframe.to_csv(DATA / f"damage_shot_requirements{suffix}.csv", index=False)
    return dataframe


def main() -> None:
    measurement = measurement_floor_study()
    damage = damage_shot_requirement_study()
    # Remark 3 sensitivity case. Fully deterministic, so it reproduces exactly and
    # cannot perturb the reported baseline figures above.
    damage_shot_requirement_study("oracle")
    summary = {
        "qwc_floor_mode1_at_100k_pct": float(
            measurement[
                (measurement["mode"] == 1)
                & (measurement["shots_per_term_independent"] == 100_000)
            ]["qwc_optimal_floor_rmse_pct"].iloc[0]
        ),
        "damage_10pct_detection_shots_per_term": float(
            damage[damage["damage_fraction"] == 0.10][
                "shots_per_term_detection_80pct_power"
            ].iloc[0]
        ),
        "damage_10pct_localization_shots_per_term": float(
            damage[damage["damage_fraction"] == 0.10][
                "shots_per_term_localization_80pct_power"
            ].iloc[0]
        ),
    }
    (DATA / "minor_revision_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(measurement.to_string(index=False))
    print(damage.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
