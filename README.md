# TQE Submission-Ready Package

This package contains the manuscript converted to the official IEEE Transactions on Quantum Engineering LaTeX template supplied by the author.

## Main files

- `Qubits_Are_Not_Enough_TQE.pdf` - compiled journal manuscript
- `Qubits_Are_Not_Enough_TQE.tex` - editable TQE LaTeX source
- `references.bib` - bibliography
- `ieeeaccess.cls`, `IEEEtran.cls` - journal class files supplied in the uploaded TQE template
- `figures/` - manuscript figures
- `reproduce_study.py`, `revision_study.py`, `minor_revision_study.py` - numerical reproduction scripts
- `validate_outputs.py` - numerical and citation checks
- `data/`, `revision_data/`, `minor_revision_data/` - numerical outputs

## Important author check

The PDF lists Bikalpa Gautam as the sole author and ETH Zurich as the affiliation. The corresponding-author email is intentionally not invented. Add your official ETH email to the `\corresp{...}` line before the final upload if you want it printed in the manuscript.

## Compile

```bash
pdflatex -interaction=nonstopmode -halt-on-error Qubits_Are_Not_Enough_TQE.tex
bibtex Qubits_Are_Not_Enough_TQE
pdflatex -interaction=nonstopmode -halt-on-error Qubits_Are_Not_Enough_TQE.tex
pdflatex -interaction=nonstopmode -halt-on-error Qubits_Are_Not_Enough_TQE.tex
```

## Reproduce numerical outputs

```bash
python -m venv .venv
source .venv/bin/activate  # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -r requirements-lock.txt
python reproduce_study.py
python revision_study.py
python minor_revision_study.py
python validate_outputs.py
```

**Note on Python version:** The original outputs were generated with Python 3.13.5 (see `environment_versions.txt`). All deterministic results (Pauli counts, densities, Gershgorin ratios, analytic shot requirements) reproduce exactly on any Python 3.10+ with the pinned dependencies. However, Monte Carlo results (finite-shot frequency RMSE, damage localization accuracy) may differ by ~5–10% on other Python versions because `scipy.optimize.minimize` (L-BFGS-B) floating-point behavior is not identical across interpreters, which shifts the random number consumption order. The `validate_outputs.py` script checks tight tolerances calibrated to the original Python 3.13.5 run and may fail on stochastic checks if run under a different version. The scientific conclusions are unaffected.

## Submission status

Reviewer-response files and historical review reports have been removed. This folder contains only the manuscript, journal support files, reproducibility code, figures, data, and submission notes.
