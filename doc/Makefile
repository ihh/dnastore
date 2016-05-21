# All source files
# Figures:
DNASTORE := dna1rep dna2rep dna2repnoend dna2noend dna2 dna3
DNASTOREPS := $(addsuffix .ps,$(DNASTORE))
DNASTORETEX := $(addsuffix .tex,$(DNASTORE))

ALLDOT := $(wildcard *.dot)
ALLLADOT := $(wildcard *.ladot)
ALLLACIRCO := $(wildcard *.lacirco)
ALLFIG := $(subst .dot,,$(ALLDOT))
ALLLAFIG := $(subst .ladot,,$(ALLLADOT)) $(subst .lacirco,,$(ALLLACIRCO))
ALLFIGPS := $(addsuffix .ps,$(ALLFIG))
ALLFIGTEX := $(addsuffix .tex,$(ALLLAFIG))
ALLFIGJUNK := $(ALLFIGTEX)

FIGURES := $(ALLFIGPS) $(ALLFIGTEX) $(DNASTOREPS) $(DNASTORETEX)

# LaTeX:
TEXINCLUDES := defs.tex notation.tex
ALLTEX := $(filter-out $(DNASTORETEX),$(filter-out $(TEXINCLUDES),$(filter-out $(ALLFIGJUNK),$(wildcard *.tex))))
ALLTEXPS := $(subst .tex,.ps,$(ALLTEX))
ALLTEXJUNK := $(addsuffix .aux,$(ALLTEX)) $(addsuffix .log,$(ALLTEX))
ALLBIB := $(wildcard *.bib)

# Top-level rules
all: figures texps

figures: $(FIGURES)

texps: $(ALLTEXPS)

bib: $(ALLBIB)

open: $(ALLTEXPS)
	open $(ALLTEXPS)

tidy:
	rm -f $(ALLFIGJUNK) $(ALLTEXJUNK) *~

clean: tidy
	rm -f $(ALLTEXPS) $(ALLFIGTEX) $(ALLFIGPS) *.dvi

# General rules:
%.open: %
	open $<

# Figure rules:
# dot
%.ps: %.dot
	dot -Tps $< >$@

# ladot
%.tex %.ps: %.ladot ladot
	./ladot $<

# lacirco
%.tex %.ps: %.lacirco lacirco
	./lacirco $<

# Hamming(7,4) transducer
hamming74.dot: hamming74.pl
	perl $< >$@

# DNAstore machines
dna1rep.lacirco: filter.pl Makefile
	../bin/dnastore -l1 -c0 -k -t0 --no-start --no-eof --dot | perl filter.pl >$@

dna2.lacirco: filter.pl Makefile
	../bin/dnastore -l2 -c1 -k --dot | perl filter.pl >$@

dna2noend.lacirco: filter.pl Makefile
	../bin/dnastore -l2 -c1 -k --no-eof --dot | perl filter.pl >$@

dna2repnoend.lacirco: filter.pl Makefile
	../bin/dnastore -l2 -c1 -k -t0 --no-eof --dot | perl filter.pl >$@

dna2rep.lacirco: filter.pl Makefile
	../bin/dnastore -l2 -c0 -k -t0 --no-start --no-eof --dot | perl filter.pl >$@

dna3.lacirco: filter.pl Makefile
	../bin/dnastore -l3 --dot | perl filter.pl >$@

# Paper rules:
%.ps: %.tex figures bib
	test -e $*.aux && rm $*.aux || eval
	latex -shell-escape $*
	bibtex $*
	latex -shell-escape $*
	latex -shell-escape $*
	dvips $*

# hey make! don't delete all my stuff
.SECONDARY:
