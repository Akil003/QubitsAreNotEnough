PYTHON ?= python
PDFLATEX ?= pdflatex
BIBTEX ?= bibtex
MAIN = Qubits_Are_Not_Enough_TQE

.PHONY: all manuscript study validate clean

all: manuscript

manuscript:
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error $(MAIN).tex
	$(BIBTEX) $(MAIN)
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error $(MAIN).tex
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error $(MAIN).tex

study:
	$(PYTHON) reproduce_study.py
	$(PYTHON) revision_study.py
	$(PYTHON) minor_revision_study.py

validate:
	$(PYTHON) validate_outputs.py

clean:
	rm -f $(MAIN).aux $(MAIN).bbl $(MAIN).blg $(MAIN).log $(MAIN).out
