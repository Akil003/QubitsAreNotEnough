PYTHON ?= python
PDFLATEX ?= pdflatex
BIBTEX ?= bibtex
MAIN = Qubits_Are_Not_Enough_TQE

.PHONY: all manuscript manuscript-tqe supplementary study validate clean

all: manuscript

# Default target: IEEE Access layout, producing $(MAIN).pdf -- the PDF shipped in
# MANIFEST.sha256. Switch the default here once the destination journal's template
# is confirmed; the source needs no other change.
manuscript:
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error $(MAIN).tex
	$(BIBTEX) $(MAIN)
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error $(MAIN).tex
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error $(MAIN).tex

# IEEE TQE layout (IEEEtran.cls) from the same source, written to a separate file so
# the two never overwrite each other. Not in the manifest; build on demand.
TQEOUT = $(MAIN)_IEEEtran
SUPP   = $(MAIN)_supplementary
manuscript-tqe:
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error -jobname=$(TQEOUT) "\def\TQE{}\input{$(MAIN)}"
	$(BIBTEX) $(TQEOUT)
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error -jobname=$(TQEOUT) "\def\TQE{}\input{$(MAIN)}"
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error -jobname=$(TQEOUT) "\def\TQE{}\input{$(MAIN)}"

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
