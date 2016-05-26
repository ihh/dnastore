# All source files
# Figures:
DNASTORE := dna1full dna2full dna2start dna2startend dna2norep dna2startnorep dna2startendnorep dna3
DNASTORELACIRCO := $(addsuffix .lacirco,$(DNASTORE))
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

FIGURES := $(ALLFIGPS) $(ALLFIGTEX) $(DNASTOREPS) $(DNASTORETEX) hamming74.ps eof2.ps

# LaTeX:
TEXINCLUDES := defs.tex notation.tex
ALLTEX := $(filter-out $(DNASTORETEX),$(filter-out $(TEXINCLUDES),$(filter-out $(ALLFIGJUNK),$(wildcard *.tex))))
ALLTEXPS := $(subst .tex,.ps,$(ALLTEX))
ALLTEXJUNK := $(addsuffix .aux,$(ALLTEX)) $(addsuffix .log,$(ALLTEX))
ALLBIB := $(wildcard *.bib)

# Top-level rules
all: figures texps

dnastore: $(DNASTOREPS)

figures: $(FIGURES)

texps: $(ALLTEXPS)

bib: $(ALLBIB)

open: $(ALLTEXPS)
	open $(ALLTEXPS)

tidy:
	rm -f $(ALLFIGJUNK) $(ALLTEXJUNK) *~

clean: tidy
	rm -f $(ALLTEXPS) $(ALLFIGTEX) *.ps *.dvi

cleanstore:
	rm -f $(DNASTORELACIRCO)

# General rules:
%.open: %
	open $<

# Figure rules:
# dot
%.ps: %.dot
	dot -Tps $< >$@

%.png: %.dot
	dot -Tpng $< >$@

# ladot
%.tex %.ps: %.ladot ladot
	./ladot $<

# lacirco
%.tex %.ps: %.lacirco lacirco
	./lacirco $<

# Hamming(7,4) transducer
hamming74.dot: hamming74.pl
	perl $< >$@

# Arithmetic code with wordlen=2
eof2.dot: mixradar.pl Makefile
	perl mixradar.pl 2 1/10 >$@

# Arithmetic code with wordlen=3
eof3.dot: mixradar.pl Makefile
	perl mixradar.pl 3 1/10 >$@

# DNAstore machines
dna1full.lacirco: filter.pl Makefile
	../bin/dnastore -l1 -c0 -t0 --no-start --no-eof --dot | perl filter.pl >$@

dna2full.lacirco: filter.pl Makefile
	../bin/dnastore -l2 -c0 -t0 --no-start --no-eof --dot | perl filter.pl >$@

dna2start.lacirco: filter.pl Makefile
	../bin/dnastore -l2 -c1 -t0 --no-eof --dot | perl filter.pl >$@

dna2startend.lacirco: filter.pl Makefile
	../bin/dnastore -l2 -c1 -t0 --dot | perl filter.pl >$@

dna2norep.lacirco: filter.pl Makefile
	../bin/dnastore -l2 -c0 --no-eof --dot | perl filter.pl >$@

dna2startnorep.lacirco: filter.pl Makefile
	../bin/dnastore -l2 -c1 --no-eof --dot | perl filter.pl >$@

dna2startendnorep.lacirco: filter.pl Makefile
	../bin/dnastore -l2 -c1 --dot | perl filter.pl >$@

dna3.lacirco: filter.pl Makefile
	../bin/dnastore -l3 --dot | perl filter.pl >$@

# Paper rules:
%.ps: %.tex figures bib $(TEXINCLUDES)
	test -e $*.aux && rm $*.aux || eval
	latex -shell-escape $*
	bibtex $*
	latex -shell-escape $*
	latex -shell-escape $*
	dvips $*

# hey make! don't delete all my stuff
.SECONDARY:
