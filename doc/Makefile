# All source files
# Figures:
ALLDOT := $(wildcard *.dot)
ALLLADOT := $(wildcard *.ladot)
ALLLACIRCO := $(wildcard *.lacirco)
ALLFIG := $(subst .dot,,$(ALLDOT))
ALLLAFIG := $(subst .ladot,,$(ALLLADOT)) $(subst .lacirco,,$(ALLLACIRCO))
ALLFIGPS := $(addsuffix .ps,$(ALLFIG))
ALLFIGTEX := $(addsuffix .tex,$(ALLLAFIG))
ALLFIGJUNK := $(ALLFIGTEX)

FIGURES := $(ALLFIGPS) $(ALLFIGTEX)

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
	rm -f $(ALLTEXPS) $(ALLFIGTEX) $(ALLFIGPS) *.dvi

# General rules:
%.open: %
	open $<

test:
	@echo $(FIGURES) $(TEXINCLUDES) $(ALLBIB)

# Paper rules:
%.ps: %.tex $(FIGURES) $(TEXINCLUDES) $(ALLBIB)
	test -e $*.aux && rm $*.aux || eval
	latex -shell-escape $*
	bibtex $*
	latex -shell-escape $*
	latex -shell-escape $*
	dvips $*

# Figure rules:
# dot
%.ps: %.dot
	dot -Tps $< >$@

# ladot
%.tex: %.ladot ladot
	./ladot $<

# lacirco
%.tex: %.lacirco lacirco
	./lacirco $<

# Hamming(7,4) transducer
hamming74.dot: hamming74.pl
	perl $< >$@

# DNAstore codes
dna2.lacirco: filter.pl
	../bin/dnastore -l2 -c1 -k --dot | perl filter.pl >$@

dna3.lacirco: filter.pl
	../bin/dnastore -l3 --dot | perl filter.pl >$@

# hey make! don't delete all my stuff
.SECONDARY:
