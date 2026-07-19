# Verification Report

## Source provenance

The source manuscript was taken from `Qubits_Are_Not_Enough_Minor_Revision_Final(1).zip`, which contained the complete 26-page minor-revision manuscript, all three reproduction scripts, generated data, and a passing SHA-256 manifest. Reviewer-response documents and review reports were intentionally excluded from this TQE package.

## Journal conversion

- Converted to the user-supplied IEEE Transactions on Quantum Engineering `ieeeaccess` LaTeX template.
- Article type prepared: Regular Article.
- Final title: **Qubits Are Not Enough: Exact Deflation Conditions and Measurement Limits in Variational Quantum Modal Analysis**.
- Author: Bikalpa Gautam.
- Affiliation: Department of Civil, Environmental and Geomatic Engineering, ETH Zurich, Zurich, Switzerland.
- Abstract length: 217 words, within the 150-250 word template requirement.
- Keywords are alphabetized.
- Bibliography uses `IEEEtran.bst` and numeric IEEE citations.
- Reviewer-response language was removed from the standalone manuscript.
- The broad title-level claim “Certified Deflation” was replaced by “Exact Deflation Conditions.” Narrow references to certified bounds remain only where an explicit mathematical bound is provided.

## Numerical integrity

`validate_outputs.py` completed successfully and reported:

```text
All citation and numerical integrity checks passed.
```

## PDF checks

- Compiled with pdfLaTeX, BibTeX, and two final pdfLaTeX passes.
- 13 TQE-formatted pages.
- PDF is unencrypted, text-based, and openable with PyMuPDF.
- All pages were rendered and visually inspected.
- Wide measurement and shot-requirement tables were converted to full-width IEEE floats to prevent overlap.
- The original template’s hard-coded 2016 footer and appendix punctuation defect were corrected locally in `ieeeaccess.cls`.
- No undefined citations or references remain.

## Remaining author-controlled item

The official ETH email address was not guessed. The manuscript identifies Bikalpa Gautam as corresponding author without printing an email. The submission portal will collect the email; add it to the `\\corresp{...}` line if you want it printed in the PDF.
