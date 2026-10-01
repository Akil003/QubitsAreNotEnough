# Final TQE Submission Checklist

- [ ] Confirm that `bgautam@ethz.ch` is the official corresponding-author email you want printed in the PDF.
- [ ] Confirm final author names, order, affiliations, corresponding author, and each
      author's approval of the submitted version.
- [ ] Confirm that every derivation, numerical result, and citation has been personally checked.
- [ ] Run `python make_submission_archive.py`; it performs the full validation, rebuilds the TQE PDFs, regenerates preflight, refreshes/verifies the manifest, and creates `tqe_submission_bundle.zip`.
- [ ] Visually inspect the freshly rebuilt `Qubits_Are_Not_Enough_TQE.pdf` page by page.
- [ ] Inspect `tqe_submission_bundle.zip` and confirm it contains no IEEE Access class/logo assets, working notes, or alternate Access PDF.
- [ ] Upload the freshly rebuilt manuscript PDF and the clean TQE source bundle through the TQE submission portal.
- [ ] Select article type: Regular Article.
- [ ] Upload code/data as supplementary material or provide a repository link, depending on the portal options.
- [ ] Retain the AI-use acknowledgment unless your actual use was limited only to grammar editing; IEEE requires disclosure of AI-generated article content.
- [ ] Do not submit the same manuscript simultaneously to the QCE workshop or another journal.
- [ ] If a conference/workshop version is later published, disclose and cite it and explain the journal extension.
