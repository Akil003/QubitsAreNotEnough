# TQE Submission-Ready Package

This package contains the manuscript, its reproduction code, and the generated data.
The source builds against two journal classes -- see **Compile** below -- and the
destination template should be confirmed before submission.

## Main files

- `Qubits_Are_Not_Enough_TQE.pdf` - compiled journal manuscript
- `Qubits_Are_Not_Enough_TQE.tex` - LaTeX source, single file for both journal targets
- `references.bib` - bibliography
- `ieeeaccess.cls`, `IEEEtran.cls` - the two journal classes the source can build against
- `figures/` - manuscript figures
- `reproduce_study.py`, `revision_study.py`, `minor_revision_study.py` - numerical reproduction scripts
- `dimension_generalization/mesh_dimension_study.py` - 2D/3D scaling and degeneracy study
- `validate_outputs.py` - numerical, citation, and regeneration checks
- `data/`, `revision_data/`, `minor_revision_data/` - numerical outputs
- `notes/` - working notes, excluded from the manifest and from any submitted archive

Not every generated artefact is used by the manuscript: `reproduce_study.py` also writes
`figures/frequency_error`, `mac_heatmap`, `damage_index`, `sensitivity_heatmap`, and
`pauli_spectrum`, plus `figures/preprocessing_scaling` from `revision_study.py`. These are
retained as diagnostics and are deliberately not referenced in the text.

## Author metadata

The manuscript lists Bikalpa Gautam (ETH Zurich) and Akil Raj Subedi (Zynga Inc.), with Bikalpa Gautam as corresponding author at `bgautam@ethz.ch`. Confirm that address on the `\corresp{...}` line before final upload.

## Compile

The manuscript builds against either of two journal classes from a **single source
file**. The body is shared verbatim; only the front matter branches (class selection,
`\PARstart`/`keywords` naming, the author block, and where `\maketitle` sits).

```bash
make manuscript        # IEEE Access layout  -> Qubits_Are_Not_Enough_TQE.pdf
make manuscript-tqe    # IEEE TQE layout     -> Qubits_Are_Not_Enough_TQE_IEEEtran.pdf
```

Equivalently, by hand:

```bash
pdflatex Qubits_Are_Not_Enough_TQE.tex                                  # Access
pdflatex "\def\TQE{}\input{Qubits_Are_Not_Enough_TQE}"                   # TQE
```

In an editor that cannot pass command-line arguments (Overleaf, TeXShop, TeXstudio),
uncomment the `% \def\TQE{}` line at the top of the `.tex` instead.

`Qubits_Are_Not_Enough_TQE.pdf` (the Access build) is the version listed in
`MANIFEST.sha256`. The IEEEtran build is produced on demand and is not manifested.
**Confirm which template the destination journal requires before submitting** --
IEEE Access and IEEE Transactions on Quantum Engineering are different journals with
different classes, and TQE supplies its own IEEEtran-based template.

## Reproduce numerical outputs

```bash
python -m venv .venv
source .venv/bin/activate  # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -r requirements-lock.txt
python reproduce_study.py
python revision_study.py
python minor_revision_study.py
python dimension_generalization/mesh_dimension_study.py   # ~1 min; needs revision_study first
python validate_outputs.py            # add --full to include the Monte Carlo studies
```

**Note on reproducibility.** All outputs were produced under one environment (Python 3.12.3; see `environment_versions.txt`), which supersedes an earlier release that mixed Python 3.13.5 and 3.12.3. The `dimension_generalization/data/` files carry an older timestamp: they are deterministic and were re-verified column by column against a fresh run rather than rewritten, so they are current in content.

All deterministic results reproduce exactly: eigenvalues, Pauli coefficients and counts, transformed densities, Gershgorin ratios, the deflation-threshold sweep, the Davis–Kahan perturbation table, the QWC measurement floor, and the analytic (delta-method) shot requirements are bit-identical, or identical to within 1e-12, against an independently produced run.

Monte Carlo results depend on the **order** of the Pauli term list, not on the interpreter version, because each term draws one binomial from a single shared `numpy` `Generator`. `pauli_coefficients` sorts by `|coefficient|`, and several benchmarks contain exactly tied magnitudes — the six-DOF Hamiltonian has two groups of four. With pandas' default non-stable quicksort those ties broke differently across environments and shifted every Monte Carlo figure by up to ~12%. Both Pauli routines now sort with `kind="stable"` over the deterministic `itertools.product` generation order, which fixes the term list and removes this source of variation. The seed is never the variable: `numpy` guarantees `Generator` reproduces the same stream for a given seed across versions and platforms.

Two smaller sources remain and are not fixable by seeding: L-BFGS-B follows a slightly different path on different BLAS builds, so the *optimizer-limited* near-zero entries of the depth study can move by a factor of a few (the expressibility-limited entries carrying the scientific claim are stable to five significant figures); and wall-clock timings are machine-specific.

## Submission status

Reviewer-response files and historical review reports have been removed. This folder contains only the manuscript, journal support files, reproducibility code, figures, data, and submission notes.
