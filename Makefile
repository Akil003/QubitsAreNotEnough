PYTHON ?= python
PDFLATEX ?= pdflatex
BIBTEX ?= bibtex
MAIN = Qubits_Are_Not_Enough_TQE

.PHONY: all manuscript manuscript-access supplementary study validate clean

all: manuscript

# CANONICAL target: IEEE TQE layout (IEEEtran.cls), producing $(MAIN).pdf -- the exact
# PDF submitted to TQE, and the artefact recorded in MANIFEST.sha256 and PDF_PREFLIGHT.txt.
manuscript:
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error $(MAIN).tex
	$(BIBTEX) $(MAIN)
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error $(MAIN).tex
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error $(MAIN).tex

# IEEE Access layout from the same source, written to a separate file so the two never
# overwrite each other. DEVELOPMENT CONVENIENCE ONLY: not manifested, not preflighted,
# and not part of the submission archive.
ACCESSOUT = $(MAIN)_IEEEAccess
SUPP      = $(MAIN)_supplementary
manuscript-access:
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error -jobname=$(ACCESSOUT) "\def\ACCESS{}\input{$(MAIN)}"
	$(BIBTEX) $(ACCESSOUT)
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error -jobname=$(ACCESSOUT) "\def\ACCESS{}\input{$(MAIN)}"
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error -jobname=$(ACCESSOUT) "\def\ACCESS{}\input{$(MAIN)}"

supplementary:
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error $(SUPP).tex
	$(BIBTEX) $(SUPP)
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error $(SUPP).tex
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error $(SUPP).tex

study:
	$(PYTHON) reproduce_study.py
	$(PYTHON) revision_study.py
	$(PYTHON) minor_revision_study.py
	$(PYTHON) dimension_generalization/mesh_dimension_study.py

validate:
	$(PYTHON) validate_outputs.py

# NB: $(MAIN).bbl is a shipped artefact listed in MANIFEST.sha256 and is deliberately
# not removed here; regenerate it with `make manuscript` if you delete it by hand.
clean:
	rm -f $(MAIN).aux $(MAIN).blg $(MAIN).log $(MAIN).out
	rm -f $(TQEOUT).aux $(TQEOUT).bbl $(TQEOUT).blg $(TQEOUT).log $(TQEOUT).out
