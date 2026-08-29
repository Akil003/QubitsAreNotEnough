# Verification Report

> **Scope.** This document records the original conversion of the manuscript to the TQE
> template. The package has since been revised; see `notes/CORRECTNESS_REVIEW.md` for the
> findings and the changes made in response, and `environment_versions.txt` for current
> provenance. Facts below have been updated where the conversion-time statement no longer
> holds.

## Source provenance

The source manuscript was taken from `Qubits_Are_Not_Enough_Minor_Revision_Final(1).zip`, which contained the complete 26-page minor-revision manuscript, all three reproduction scripts, generated data, and a passing SHA-256 manifest. Reviewer-response documents and review reports were intentionally excluded from this TQE package.

## Journal conversion

- Converted to the user-supplied IEEE Transactions on Quantum Engineering `ieeeaccess` LaTeX template.
- Article type prepared: Regular Article.
- Final title: **Qubits Are Not Enough: Deflation and Measurement Limits in Variational Quantum Modal Analysis**.
- Authors: Bikalpa Gautam (ETH Zurich), Akil Raj Subedi (Zynga Inc.); corresponding author Bikalpa Gautam, `bgautam@ethz.ch`.
- Affiliation: Department of Civil, Environmental and Geomatic Engineering, ETH Zurich, Zurich, Switzerland.
- Abstract length: 249 words, within the 150-250 word template requirement.
- Keywords are alphabetized.
- Bibliography uses `IEEEtran.bst` and numeric IEEE citations.
- Reviewer-response language was removed from the standalone manuscript.
- The broad title-level claim “Certified Deflation” was replaced by “Exact Deflation Conditions.” Narrow references to certified bounds remain only where an explicit mathematical bound is provided.

## Numerical integrity

`validate_outputs.py` completed successfully and reported:

```text
All citation, numerical, and complexity integrity checks passed.
```

Note that the tolerances on the two stochastic checks were subsequently widened; see the
reproducibility note in `README.md` for why four-digit agreement is not achievable on
Monte Carlo quantities.

## PDF checks

- Compiled with pdfLaTeX, BibTeX, and two final pdfLaTeX passes.
- 19 TQE-formatted pages.
- PDF is unencrypted, text-based, and openable with PyMuPDF.
- All pages were rendered and visually inspected.
- Wide measurement and shot-requirement tables were converted to full-width IEEE floats to prevent overlap.
- The original template’s hard-coded 2016 footer and appendix punctuation defect were corrected locally in `ieeeaccess.cls`.
- No undefined citations or references remain.

## Remaining author-controlled item

The corresponding-author address `bgautam@ethz.ch` is now printed on the `\corresp{...}` line. Confirm it is the address you want to appear before final upload.
