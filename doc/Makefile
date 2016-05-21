# All source files
# Figures:
ALLDOT := $(wildcard *.dot)
ALLLADOT := $(wildcard *.ladot)
ALLFIG := $(subst .dot,,$(ALLDOT))
ALLLAFIG := $(subst .ladot,,$(ALLLADOT))
ALLFIGPNG := $(addsuffix .png,$(ALLFIG))
ALLFIGTEX := $(addsuffix .tex,$(ALLLAFIG))
ALLFIGJUNK := $(ALLFIGTEX)

FIGURES := $(ALLFIGPNG) $(ALLFIGTEX)

# LaTeX:
TEXINCLUDES := defs.tex notation.tex
ALLTEX := $(filter-out $(TEXINCLUDES),$(filter-out $(ALLFIGJUNK),$(wildcard *.tex)))
ALLTEXPS := $(subst .tex,.ps,$(ALLTEX))
ALLTEXJUNK := $(addsuffix .aux,$(ALLTEX)) $(addsuffix .log,$(ALLTEX))
ALLBIB := $(wildcard *.bib)

# Top-level rules
all: $(FIGURES) $(ALLTEXPS)

open: $(ALLTEXPS)
	open $(ALLTEXPS)

tidy:
	rm -f $(ALLFIGJUNK) $(ALLTEXJUNK) *~

clean: tidy
	rm -f $(ALLTEXPS) $(ALLFIGTEX) $(ALLFIGPNG) *.dvi

# General rules:
%.open: %
	open $<

# Paper rules:
%.ps: %.tex $(FIGURES) $(TEXINCLUDES) $(ALLBIB)
	test -e $*.aux && rm $*.aux || eval
	latex -shell-escape $*
	bibtex $*
	latex -shell-escape $*
	latex -shell-escape $*
	dvips $*

# Figure rules:
# PNGs (requires graphviz only)
%.png: %.dot
	dot -Tpng $< >$@

# ladot
%.tex: %.ladot
	./ladot $<

# hey make! don't delete all my stuff
.SECONDARY:
