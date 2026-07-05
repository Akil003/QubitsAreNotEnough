"""Lightweight integrity checks for the manuscript's generated outputs."""
from __future__ import annotations

import json
import re
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


def main() -> None:
    validate_citations()
    validate_numbers()
    print("All citation and numerical integrity checks passed.")


if __name__ == "__main__":
    main()
