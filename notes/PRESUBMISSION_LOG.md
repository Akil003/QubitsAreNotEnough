# Pre-submission log — IEEE_TQE_MASTER_PRE_SUBMISSION_CORRECTIONS.txt

Tracks the master checklist dated 14 Sep 2026. Separate from `REVISION_LOG.md`, which
tracks the earlier `QANE_TQE_revision_brief.md`.

Status: `todo` · `done` · `waived` (with reason) · `blocked` (needs author input)
· `optional` (Section G / L, not correctness blockers)

Baseline for this round: `b61d8b1`.

---

## A. Must-fix technical / mathematical

| ID | Item | Status | Notes |
|----|------|--------|-------|
| A1 | Generalized-quotient asymptotic scaling mis-stated as constant-factor | done | **My error, refuted by my own data.** Fitted exponents over `N>=16`: **1.475** for `N_P(A)` vs **1.052** for `N_P(K)+N_P(M)`; ratio falls monotonically 1.60 → 0.254, and `N_P(KM)/(N log2 N)` is slowly *declining*, consistent with `O(N log N)`. Replaced the closing paragraph of §VIII-E: the baseline changes from `O(kN^6 N_eval)` (whitened worst case) to **`O(kN^5 log N N_eval)`** (unwhitened pair) — an asymptotic improvement in the *representation factor*. Explicitly labelled a representation-plus-shot baseline, **not** an end-to-end complexity (excludes all four unmeasured costs; needs a solver + covariance-aware shot analysis). Practical conclusion restated as unchanged: still far above `O(kN T_it)`. Added a formulation-dependence pointer at the §VIII-D `O(kN^6)` statement. |
| A2 | Denominator claimed unbounded below; M is SPD | done | **My error.** `<psi\|M\|psi> >= lambda_min(M) > 0` for every normalized state, so positivity *is* guaranteed; the real issue is conditioning. Rewritten with the `1/kappa(M)` bound, and **quantified**: `kappa(M) <= 2.50` lumped, `<= 4.92` consistent across `N<=128`, both plateauing ⇒ normalized denominator bounded below by ~0.2. Turns the caveat into a measured, reassuring statement. |
| A3 | `Omega(N/log N)` depth stated without generic/universal qualifier | done | Complexity-table footnote now separates (i) generic real-state reachability of *this* ansatz by parameter counting from (ii) worst-case arbitrary amplitude-state preparation, and states explicitly it is not a lower bound for any individual structural eigenstate. |
| A4 | "Symmetry-protected" degeneracy is not supported | done | **My error, and my own measurement disproved it.** Verified computationally: of the 8 square-symmetry operations, only **identity + diagonal reflection** fix the pinned corner and preserve both `K` and `M` ⇒ group is `Z_2`, which has only 1-D irreps and therefore *cannot* force a twofold degeneracy. The `diag(-1,+1)` action I had cited as evidence *for* protection actually places the two modes in **different** 1-D irreps — evidence against. Multiplicity also grows as `n_side - 1` at a fixed `lambda = 6.6e4` (2,3,4,5 for n_side 3,4,5,6), which no point group produces. Reworded as an exact *algebraic* degeneracy that the anisotropy splits. Also reordered the practical criterion to lead with the rigorous `2||Delta|| < G_r` condition, demoting the order-of-magnitude heuristic to a labelled empirical feature of this sweep. |
| A5 | Categorical "consistent mass ... unusable" | done | Scoped to "the whitened, independently measured 1e5-shots-per-term configuration tested here", with an explicit pointer to §VIII-E showing 750 unwhitened terms vs 2950 whitened. |

## B. High-priority scientific / reviewer-proofing

| ID | Item | Status | Notes |
|----|------|--------|-------|
| B1 | Update literature through Sep 2026 (4 refs) | done | **All four verified via Crossref** before any citation, per the D2 protocol. Lee & Kanno **QCE 2025, pp. 2170–2180, DOI 10.1109/QCE65121.2025.00237** — this *contradicts* the earlier brief, which said the QCE 2025 claim "does not match the record"; that warning was stale, the conference record is real, so the `@misc` arXiv entry should be replaced by the peer-reviewed one. Xu et al., *JMPS* **213**, 106639 (2026), DOI 10.1016/j.jmps.2026.106639. Wagatsuma, Endo & Terada, *IJNME* **127**(16), e70407 (2026), DOI 10.1002/nme.70407. Xu & Hu, *Comput. Struct.* **330**, 108309 (2026), DOI 10.1016/j.compstruc.2026.108309. Applied: `LeeKanno2023` @misc replaced by `@inproceedings{LeeKanno2025}`; three new `@article` entries added; Introduction positioning paragraph added after the feasibility sentence, framing all four as evidence that representation/state-preparation/output precision are algorithm-dependent and explicitly scoping this paper's negative conclusions to the variational setting. Bibliography 36 → 39 entries, all cited, 0 undefined. Also braced `{H}amiltonian` in two titles, since `IEEEtran.bst` lowercases mid-title words and had rendered "hamiltonian". |
| B2 | State what the generalized-quotient study does/doesn't prove | done | Explicit scope paragraph added **early** in §VIII-E (before the results, not only after them): representation cost only; no end-to-end generalized-eigenvalue solver; ratio-estimator finite-shot bias/variance not measured; no `M`-metric deflation threshold or its measurement cost derived. Closes with "establishes that whitening can be avoided at the representation level; it does not establish a total runtime". Complements the A1 rewrite, which states the same limits after the asymptotic claim. |
| B3 | Trainability claim conservative everywhere | done | Audited: **zero** occurrences of the forbidden phrasings ("we prove a barren plateau", "the barren-plateau rate is reached"), and no claim that the `n_q=10` target requires `L>=101`. Contribution (iv) tightened to B3's preferred wording: "a **finite-size** decay rate comparable to `2^-n_q` **at L=8**", retaining "universal real-state expressibility". |
| B4 | Preserve VQD novelty scope (Higgott attribution) | done | **Verified, no edit needed.** `Higgott2019` cited at both threshold sites; the Introduction states the sufficiency condition "is theirs, not ours" and lists this paper's claim as the converse plus padded-spectrum extension, certified rule, and approximate-projector certificate — i.e. the stronger contribution is not sold as the converse alone. |
| B5 | beta-star always "certificate-minimizing" | done | **Verified, no edit needed.** Every `beta*` mention is certificate-attached: "the penalty that minimizes the approximate-projector certificate" (contributions), "minimized at" + "It locates the minimum of the *certificate*, not of the error" (§V), "certificate-optimal" (§XI-E), "locates the penalty that minimizes it [the certificate]" (Conclusions). No bare "optimal penalty" anywhere. |
| B6 | Damage-localization scope explicit everywhere | done | **Real gap found in the figure caption.** `fig:shotdamage` said only "Correct damage-localization probability from independently sampled element-energy observables" — no oracle-state caveat, despite carrying the 84.5% headline. Caption now names Benchmark~C and 200 trials/condition and states the modal states are the **exact** healthy and damaged eigenvectors, not variationally obtained, so the curves isolate measurement noise and are an **upper bound** on an end-to-end pipeline. Claim verified against the code: `damage_shot_study` calls `exact_modal_data()` for both epochs, which does `eigh(A)`. Conclusion quotes only the `~1e8` design estimate (no 84.5%), so needs no caveat; cover letter already states oracle access. |
| B7 | settings != total shots | done | **Verified, no edit needed.** Text reads "$8\times19=152$ measurement settings rather than 648, a factor of $4.3$ reduction in the number of **distinct circuits that must be prepared**" — explicitly a settings reduction, with no proportional shot claim anywhere. |
| B8 | RCM terminology | done | Both sites fixed. §XI-C now reads "a standard bandwidth- and profile-reducing sparse-matrix ordering", and adds explicitly that RCM **is not** a fill-minimizing ordering in the sense of AMD or nested dissection — it is used because it is the common default — while noting it does reduce the consistent-mass fill measurably. §VIII-A's forward pointer changed to "bandwidth-reducing reordering" and now attributes the no-fill-benefit claim specifically to the lumped case, which is where our data actually shows it. |

## C. Reproducibility / code / data

| ID | Item | Status | Notes |
|----|------|--------|-------|
| C1 | Validation wording must match what the validator actually reruns | done | Enumerated the actual coverage from the `diff()` call sites: **10 datasets regeneration-diffed by default, 5 more under `--full`** (15 total), everything else bound by assertion to generated summaries. Data Availability now states exactly that, replacing "re-executes the deterministic generators", which over-read as *all* of them. |
| C2 | Diff `trainability_depth_sweep_ci.csv` under `--full` | done | Added immediately after the `trainability_depth_sweep()` rerun, so the bootstrap CI artefact is itself regeneration-checked rather than only spot-asserted through the quoted bounds. `--full` passes. |
| C3 | Per-term RNG substreams for order-invariant Monte Carlo | optional | Disclosed already; not a blocker. |
| C4 | Cross-platform CI | optional | |
| C5 | Every headline number generator-backed | done | Swept 188 numeric tokens in prose against all code, CSV and JSON artefacts. 7 appeared unbacked; **all 7 resolved on inspection** — 127.99/130.20 are roundings of validator-asserted restart values, 16.7% is `1/6`, 2.15e-3 is the asserted `norm_delta` max, 2025 is a citation year, and 5.16%/9.82e-4 are roundings of `larger_vqe_benchmarks.csv` (5.158516) and `gradient_trainability.csv` (0.00098197). **Net: 188/188 trace to a generator.** Four of these had a generator but no assertion binding the printed rounding, so assertions were added for the 16-DOF depth-1/depth-2 errors and the 2- and 12-qubit gradient variances. |

## D. Manuscript claim / writing cleanup

| ID | Item | Status | Notes |
|----|------|--------|-------|
| D1 | Abstract final check | done | **Verified against the rendered PDF**, not the source: 244 words (within 150–250), 0 citations, 0 footnotes. Contains the oracle-state disclosure, "known deflation" attribution, the eighty-percent-power design-estimate framing, and "claim no speedup". The generalized quotient is not mentioned at all, so nothing implies it was solved end to end. |
| D2 | Contribution list: keep four; scope (iii) after A1 | done | Still four contributions. (iii) now reads "end-to-end resource accounting *for the whitened variational pipeline*", so the A1 correction cannot be misread as claiming the generalized-quotient route was accounted end to end. (ii) retains the locality-preserving nearest-neighbour scope; (iv) retains "universal real-state" and gained the B3 finite-size qualifier. |
| D3 | Soften categorical language | done | Swept: `unusable` **0**, `impossible` **0**, `proves` **0** (A5 removed the last "unusable"). The single `dominates` is a correct statement about which row of the asymptotic table is largest, not an empirical overreach, so it was kept deliberately. |
| D4 | Conclusion foregrounds three durable conclusions | done | Final paragraph restructured to state all three explicitly: (1) logarithmic qubit storage does not imply logarithmic circuit or measurement cost; (2) deflation is spectrally certifiable but approximate projectors and small gaps make penalty selection and **block treatment** of near-degenerate pairs consequential; (3) measurement/decision precision is the binding resource and no advantage is established — now explicitly "a statement about this algorithmic regime, not about the fault-tolerant, qubitization and Hamiltonian-simulation routes surveyed in §I", which is exactly D4's warning and ties to the B1 additions. The title's payoff sentence is retained. |
| D5 | Title unchanged | done | "Variational" retained. |
| D6 | Reduce repetition (1–2 pp, optional) | optional | Now 23 pp Access / 22 pp IEEEtran after A1 added a paragraph. |

## E. IEEE / TQE format and submission package

| ID | Item | Status | Notes |
|----|------|--------|-------|
| E1 | Make IEEEtran the default and manifested build | done | **Resolves the previously blocked F1** of the earlier brief. Implemented by *inverting the flag* rather than deleting the dual target: the guard is now `\ifdefined\ACCESS`, so building with no flags gives IEEEtran. This makes `Qubits_Are_Not_Enough_TQE.pdf` — already manifested and already named `_TQE` — actually *be* the TQE build (22 pp, US letter), so no manifest path surgery was needed. The Access layout writes to `..._IEEEAccess.pdf`, is gitignored, not manifested, not preflighted, and explicitly not part of the submission archive. Stale `_IEEEtran.*` artefacts deleted; Makefile target renamed `manuscript-access`; README build section rewritten. One defect found and fixed during the flip: the author block had an explicit `\else`, so flipping the guard mechanically produced `Extra \else` — the two branches had to be swapped instead. Front matter verified in the rendered PDF. |
| E2 | Acknowledgment after appendices, before References | done | Order is now COI → Appendices → Acknowledgment → References, matching TQE Template v4. Verified 0 undefined references after the move (the disclosure's `\cref` to §VIII-E/XI-C/XI-E still resolves). |
| E3 | Refine AI disclosure | done | Checked against IEEE PSPB 8.2.1.B: system named (Claude, Anthropic), affected sections identified (§VIII-E, §XI-C, §XI-E plus revised passages and the analysis/validation scripts), level of use stated, and it does **not** imply the AI independently verified correctness — the audit is described as AI-assisted with findings "assessed before adoption". Author responsibility explicit; "No AI system is an author". **Author decision on file:** no version string added — "Claude (Anthropic)" judged sufficient. |
| E4 | PDF_PREFLIGHT / VERIFICATION_REPORT consistency | done | `VERIFICATION_REPORT.md` corrected from "19 TQE-formatted pages" to 22 pp IEEEtran letter plus the 1-pp supplement, and now points at `PDF_PREFLIGHT.txt` as the generated source of truth. `make_preflight.py` retargeted from the retired `_IEEEtran` name to the two artefacts actually submitted (main + supplement); the dev-only Access build is deliberately excluded. Both files now agree at 22 pp. |
| E5 | SUBMISSION_CHECKLIST "sole-authorship" wording | done | Replaced with author names/order/affiliations/corresponding author/per-author approval. |
| E6 | Verify author metadata | blocked | Author input. |
| E7 | Submit as Regular Article | blocked | Author action at portal. |
| E8 | APC planning (US$1,995 listed) | blocked | Author action. |
| E9 | Submission archive hygiene | todo | |
| E10 | Regenerate manifest last | todo | Ordering discipline. |

## F. Figures / tables / PDF quality

| ID | Item | Status | Notes |
|----|------|--------|-------|
| F1 | Zero Type-3 / unembedded fonts | done | 0/0 in all three PDFs; recheck after every regeneration (matplotlib can revert). |
| F2 | Verify wide tables under IEEEtran | done | All 9 tables inspected in the **canonical IEEEtran** build (0 overfull boxes). The two widest single-column tables — 8-column `tab:genquotient` (p. 10) and 7-column `tab:neardeg` (p. 17) — sit comfortably within the narrower IEEEtran column at `\small`; `tab:ordering` (p. 16) and the two full-width `table*` floats (`tab:shotfreq`, `tab:dimscaling`) render correctly. Captions explain dashes and thresholds where relevant, and precision is justified (integer counts, ratios to 2 dp, densities to 3 dp). **The visual pass caught a cross-subsection contradiction:** §XI-D still said the diagonal reflection *causes* the degeneracy ("...is an exact symmetry, **so** a uniform square grid has genuinely degenerate mode pairs"), which A4 had just disproved two subsections later. Reworded to an algebraic repetition with a pointer to §XI-E. |
| F3 | Figure-caption self-containment | done | Audited all 7 surviving figures against F3's criteria (benchmark, sample count, interval meaning, coefficient threshold, and whether exact / Monte Carlo / analytic / hardware). Five were incomplete and were extended: `fig:penalty_threshold` (+Benchmark C, exact spectrum, no sampling); `fig:density` (+Test B family, exact, and the `tau=1e-12` *relative entry* cut); `fig:pauliscale` (+Test B family, exact, the `1e-10` *absolute coefficient* cut, and a pointer that the consistent-mass count is cut-sensitive); `fig:largervqe` (+Test B members, exact state-vector, so residual error is ansatz/optimizer only); `fig:gradients` (+Test B chains, **80 draws per qubit count**, and that the intervals reflect parameter-distribution sampling rather than shot noise). `fig:shotdamage` was already completed under B6; `fig:workflow` is a schematic, so the data criteria do not apply. |
| F4 | Fitted exponents not asymptotic laws | done | Already labelled illustrative; 3D three-size caveat in body and `tab:dimscaling` caption. |

## G / L. Optional experiments

| ID | Item | Status |
|----|------|--------|
| G1 | End-to-end generalized-quotient VQE benchmark | optional (highest value) |
| G2 | AMD / nested-dissection ordering | optional |
| G3 | Vector-DOF FE benchmark | optional |
| G4 | Real-hardware demo | optional (explicitly low priority) |
| G5 | Covariance-aware damage power calculation | optional |
| G6 | Per-location finite-shot damage sweep | optional |
| G7 | Second near-degenerate geometry | optional |

## H / I. Cover letter and final literature check

| ID | Item | Status | Notes |
|----|------|--------|-------|
| H | Cover-letter strategy | done | Rewritten to the checklist's pitch: leads with "evaluation at the level an engineering decision requires, rather than at the level of qubit count", then states explicitly what it is **not** (not a quantum-advantage paper, not deployable SHM, not "VQE applied to structural engineering"). Adds a dedicated prior-art/novelty paragraph naming Higgott for sufficiency and the converse + certified rule + approximate-projector certificate as ours; a reproducibility paragraph with the corrected 10/5 regeneration split and the 7-of-7 injection result; and an alternative-routes paragraph covering all four new references so the negative conclusions are visibly scoped. **Three stale items fixed:** the B8 "fill-reducing" wording had leaked into the letter (now bandwidth-reducing), the B7 settings-vs-shots distinction was missing, and the pre-send checklist still listed the AI disclosure and the build-default question as open — both now closed, with ORCIDs/APC/article-type/archive-URL remaining. |
| I | Final literature search immediately before upload | done (this round) | Four 2026/2025 refs verified via Crossref on 17 Sep 2026. Repeat in the upload week. |
