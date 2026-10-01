# Cover letter — draft (NOT part of the submission package)

Working draft for D4 of the revision brief and H of the master pre-submission checklist. This file is a note to the authors, not
manuscript text and not an artefact to be uploaded. Author names, ORCIDs, and the
corresponding-author details are left as placeholders pending F4.

---

Dear Editor,

We submit *"Qubits Are Not Enough: Deflation and Measurement Limits in Variational Quantum
Modal Analysis"* for consideration in *IEEE Transactions on Quantum Engineering*.

**What the paper is.** It is a quantum-engineering evaluation of variational modal analysis
carried out at the level an engineering decision requires, rather than at the level of qubit
count. It explains why logarithmic state encoding is insufficient, sharpens the exact
variational-deflation condition, quantifies the robustness lost to approximate projectors,
derives the Pauli representation structure of locality-preserving lumped-mass meshes, measures
how the mass formulation, whitening and degree-of-freedom ordering change the observable cost,
and propagates finite-shot uncertainty all the way through to a structural damage decision.

It is deliberately *not* a quantum-advantage paper, not a deployable structural-health-monitoring
method, and not "VQE applied to structural engineering". The numerical program assumes ideal
state-vector circuits wherever possible, so the negative resource conclusions are lower-bound
statements made under assumptions favourable to the quantum method.

**Prior art is separated from novelty explicitly.** The Cholesky mass-orthonormal mapping and
the two-scale normalization are adopted, not claimed. Higgott et al.\ established the
sufficiency of the deflation penalty condition, and we say so at every point the threshold
appears; our contribution is the converse and its extension to the padded spectrum, a
deterministic certified penalty rule, and the approximate-projector certificate together with
the penalty that minimizes it. Four contributions are claimed, each either proved or measured.

**Why TQE.** The venue has published both variational-algorithm optimizer studies and VQE
measurement-cost analyses — including Gokhale et al., *IEEE Trans. Quantum Eng.* **1**,
1–24 (2020), whose `O(N^4) -> O(N^3)` Pauli-partitioning result we cite and build on in the
measurement-accounting section. Both the method class and the negative-result form are
therefore within demonstrated scope for the journal.

**Reproducibility.** Every number in the manuscript is produced by a shipped generator. The
package includes a `MANIFEST.sha256` covering all artefacts and a validator that re-executes a
deterministic regression suite into temporary directories and diffs it against the committed
outputs at zero tolerance -- ten datasets by default and five more under `--full`, which also
reruns the Monte Carlo studies -- in addition to binding the manuscript-quoted values to their
generated summaries by assertion. The validator was itself fault-injection tested, catching 7
of 7 deliberately introduced defects with no false positive on the clean package.

**Recent alternative routes are discussed.** The introduction positions the work against the
fault-tolerant and Hamiltonian-simulation literature -- qubitization and phase estimation for
modal analysis, energy-conservation Hamiltonian simulation of elastodynamics with a hardware
demonstration, QSVT-based simple harmonic motion analysis, and a decomposition-free
potential-energy formulation. Our conclusions are therefore explicitly limited to the
variational setting studied, and the Conclusion says so in those terms.

**Scope and limitations, stated up front.** The numerical program is deliberately
idealized: exact circuit states, no gate or readout noise, and — for the damage-localization
study — oracle access to the damaged modes. We treat these as upper bounds on achievable
performance rather than as feasibility evidence, and say so in the abstract. The
dimensional and Pauli-count results are additionally conditioned on the locality-preserving
DOF ordering studied here; we show that a standard bandwidth-reducing reordering degrades the
Pauli count by up to an order of magnitude, and that a scrambled ordering saturates the
real-symmetric worst case exactly. Where measurement savings are reported they are savings in
distinct circuit settings, not a proportional reduction in total shots; the covariance-aware
power analysis that would license the latter is named as future work.

We confirm the manuscript is original, is not under consideration elsewhere, and that all
authors have approved the submission.

Sincerely,
[Author names]
[Corresponding author, affiliation, email, ORCID]

---

## Checklist before sending

- [ ] Fill author names, affiliations, corresponding-author email, **ORCIDs** (E6; still outstanding)
- [ ] Confirm the **APC** arrangement with the institution before submitting (E8; TQE lists US$1,995)
- [ ] Submit as a **Regular Article** (E7; Review/Tutorial are invitation-based)
- [ ] Decide the **archive URL** for the Data Availability statement, and confirm it does not
      expose `notes/` (E9 / the earlier brief's F2)
- [x] AI-use disclosure present in the Acknowledgment, per IEEE PSPB 8.2.1.B
- [x] Canonical build settled: IEEEtran is the default; `Qubits_Are_Not_Enough_TQE.pdf` is the
      submitted PDF and the manifested artefact (resolves the earlier brief's F1)
- [ ] Repeat the literature search in the upload week (checklist section I)
