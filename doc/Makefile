# All source files
# Figures:
ALLDOT := $(wildcard *.dot)
ALLFIG := $(subst .dot,,$(ALLDOT))
ALLFIGPNG := $(addsuffix .png,$(ALLFIG))

# Source files that only render properly as PNGs
FIGURES := $(ALLFIGPNG)

# LaTeX:
TEXINCLUDES := defs.tex notation.tex
ALLTEX := $(filter-out $(TEXINCLUDES),$(filter-out $(ALLFIGJUNK),$(wildcard *.tex)))
ALLTEXPDF := $(subst .tex,.pdf,$(ALLTEX))
ALLTEXJUNK := $(addsuffix .aux,$(ALLTEX)) $(addsuffix .log,$(ALLTEX))
ALLBIB := $(wildcard *.bib)

# Top-level rules
all: $(FIGURES) $(ALLTEXPDF)

open: $(ALLTEXPDF)
	open $(ALLTEXPDF)

tidy:
	rm -f $(ALLFIGJUNK) $(ALLTEXJUNK) *~

clean: tidy
	rm -f $(ALLTEXPDF) $(ALLFIGPDF) $(ALLFIGPNG)

# General rules:
%.open: %
	open $<

# Paper rules:
%.pdf: %.tex $(FIGURES) $(TEXINCLUDES) $(ALLBIB)
	test -e $*.aux && rm $*.aux || eval
	pdflatex -shell-escape $*
	bibtex $*
	pdflatex -shell-escape $*
	pdflatex -shell-escape $*

# Figure rules:
# PNGs (requires graphviz only)
%.png: %.dot
	dot -Tpng $< >$@

# hey make! don't delete all my stuff
.SECONDARY:
