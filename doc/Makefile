# Makefile dir
TOP := $(dir $(CURDIR)/$(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST)))

# LDPC
LDPCDIR := $(TOP)/../../LDPC-codes

# All source files
# Figures:
DNASTORE := dna1full dna2full dna2start dna2startend dna2norep dna2startnorep dna2startendnorep dna3
DNASTORECIRCO := $(addsuffix .circo,$(DNASTORE))
DNASTOREPS := $(addsuffix .ps,$(DNASTORE))

ALLDOT := $(wildcard *.dot)
ALLCIRCO := $(wildcard *.circo)
ALLFIG := $(subst .dot,,$(ALLDOT)) $(subst .circo,,$(ALLCIRCO))
ALLFIGPS := $(addsuffix .ps,$(ALLFIG))

FIGURES := $(ALLFIGPS) eof2.ps subplot.ps dupplot.ps delplot.ps watplot.ps

# LaTeX:
TEXINCLUDES := defs.tex notation.tex mixradtab.tex dnastoretab.tex
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

# circo
%.ps: %.circo
	circo -Tps $< >$@

# Hamming(7,4) transducer
hamming74.dot: ../bin/hamming74.pl
	perl $< >$@

# Arithmetic code with wordlen=2
eof2.dot: ../bin/mixradar.pl
	perl ../bin/mixradar.pl -k 2 1/100 >$@

# Arithmetic code with wordlen=3
eof3.dot: ../bin/mixradar.pl
	perl ../bin/mixradar.pl -k 3 1/100 >$@

# DNAstore machines
dna1full.circo: filter.pl
	../bin/dnastore -l1 -c0 -t0 --no-start --no-end --dot | perl filter.pl >$@

dna2full.circo: filter.pl
	../bin/dnastore -l2 -c0 -t0 --no-start --no-end --dot | perl filter.pl >$@

dna2start.circo: filter.pl
	../bin/dnastore -l2 -c1 -t0 --no-end --dot | perl filter.pl >$@

dna2startend.circo: filter.pl
	../bin/dnastore -l2 -c1 -t0 --dot | perl filter.pl >$@

dna2norep.circo: filter.pl
	../bin/dnastore -l2 -c0 --no-end --dot | perl filter.pl >$@

dna2startnorep.circo: filter.pl
	../bin/dnastore -l2 -c1 --no-end --dot | perl filter.pl >$@

dna2startendnorep.circo: filter.pl
	../bin/dnastore -l2 -c1 --dot | perl filter.pl >$@

dna3.circo: filter.pl
	../bin/dnastore -l3 --dot | perl filter.pl >$@

# Results data tables
# Code sizes:
codetabs: mixradtab.tex dnastoretab.tex

mixradtab.tex: mixradtab.pl
	perl $< >$@

dnastoretab.tex: dnastoretab.pl
	perl $< >$@

# Simulations of decoding with errors:
errtabs: deltabs duptabs subtabs

deltabs: len4.dels.tab len4.mix2.dels.tab len4.ham.dels.tab len4.wat64.1.ldpc.dels.tab len4.wat64.16.ldpc.dels.tab len4.wat64.ldpc.dels.tab len4.wat128.1.ldpc.dels.tab

duptabs: len4.mix2.dups.tab len8.mix2.dups.tab

subtabs: len4.subs.tab len4.mix2.subs.tab len4.ham.subs.tab len4.wat64.1.ldpc.subs.tab

w128deltabs: len4.wat128.ldpc.dels.tab len4.wat128.16.ldpc.dels.tab len4.wat128.16b.ldpc.dels.tab len4.wat128n.ldpc.dels.tab len4.wat128a2.ldpc.dels.tab

w64deltabs: len4.wat64.ldpc.dels.tab len4.wat64.16.ldpc.dels.tab len4.wat64.16b.ldpc.dels.tab len4.wat64n.ldpc.dels.tab len4.wat64a2.ldpc.dels.tab

LOGSCALE = .001,.0014,.002,.0028,.004,.0057,.008,.011,.016,.023,.032,.045,.064,.09,.128,.18,.26,.36,.5
ERRDECODE = perl errdecode.pl -bits 8192 -trainalign 10 -reps 20 -u 3

len4.subs.tab:
	$(ERRDECODE) -length 4 -subrate 0,$(LOGSCALE) -delrate 0 -duprate 0 >$@

len4.dels.tab:
	$(ERRDECODE) -length 4 -subrate 0 -delrate $(LOGSCALE) -duprate 0 >$@

len4.mix2.subs.tab:
	$(ERRDECODE) -length 4 -mix2 -subrate 0,$(LOGSCALE) -delrate 0 -duprate 0 >$@

len4.mix6.subs.tab:
	$(ERRDECODE) -length 4 -mix6 -subrate 0,$(LOGSCALE) -delrate 0 -duprate 0 >$@

len4.mix2.dels.tab:
	$(ERRDECODE) -length 4 -mix2 -subrate 0 -delrate $(LOGSCALE) -duprate 0 >$@

len4.mix2.dups.tab:
	$(ERRDECODE) -length 4 -mix2 -subrate 0 -delrate 0 -duprate $(LOGSCALE) -maxdupsize 4 >$@

len8.mix2.dups.tab:
	$(ERRDECODE) -length 8 -mix2 -subrate 0 -delrate 0 -duprate $(LOGSCALE) -maxdupsize 4 >$@

len4.ham.subs.tab:
	$(ERRDECODE) -length 4 -hamming -subrate 0,$(LOGSCALE) -delrate 0 -duprate 0 >$@

len4.ham.dels.tab:
	$(ERRDECODE) -length 4 -hamming -subrate 0 -delrate 0,$(LOGSCALE) -duprate 0 >$@

len4.sync64.mix2.dels.tab:
	$(ERRDECODE) -length 4 -syncfreq 64 -mix2 -subrate 0 -delrate 0,$(LOGSCALE) -duprate 0 >$@


len4.wat128.ldpc.subs.tab:
	$(ERRDECODE) -length 4 -watermark 128 -ldpc $(LDPCDIR) -subrate 0,$(LOGSCALE) -delrate 0 -duprate 0 >$@

len4.wat128.ldpc.dups.tab:
	$(ERRDECODE) -length 4 -watermark 128 -ldpc $(LDPCDIR) -subrate 0 -delrate 0 -duprate $(LOGSCALE) >$@

len4.wat128.16.ldpc.subs.tab:
	$(ERRDECODE) -length 4 -watermark 128.16 -ldpc $(LDPCDIR) -subrate 0,$(LOGSCALE) -delrate 0 -duprate 0 >$@

len4.wat128.16.ldpc.dups.tab:
	$(ERRDECODE) -length 4 -watermark 128.16 -ldpc $(LDPCDIR) -subrate 0 -delrate 0 -duprate $(LOGSCALE) >$@


len4.wat128.ldpc.dels.tab:
	$(ERRDECODE) -length 4 -watermark 128 -ldpc $(LDPCDIR) -subrate 0 -delrate $(LOGSCALE) -duprate 0 >$@

len4.wat128.16.ldpc.dels.tab:
	$(ERRDECODE) -length 4 -watermark 128.16 -ldpc $(LDPCDIR) -subrate 0 -delrate $(LOGSCALE) -duprate 0 >$@

len4.wat128.16b.ldpc.dels.tab:
	$(ERRDECODE) -length 4 -watermark 128.16b -ldpc $(LDPCDIR) -subrate 0 -delrate $(LOGSCALE) -duprate 0 >$@

len4.wat128n.ldpc.dels.tab:
	$(ERRDECODE) -length 4 -watermark 128n -ldpc $(LDPCDIR) -subrate 0 -delrate $(LOGSCALE) -duprate 0 >$@

len4.wat128a2.ldpc.dels.tab:
	$(ERRDECODE) -length 4 -watermark 128a2 -ldpc $(LDPCDIR) -subrate 0 -delrate $(LOGSCALE) -duprate 0 >$@

len4.wat128.1.ldpc.dels.tab:
	$(ERRDECODE) -length 4 -watermark 128.1 -ldpc $(LDPCDIR) -subrate 0 -delrate $(LOGSCALE) -duprate 0 >$@


len4.wat64.ldpc.dels.tab:
	$(ERRDECODE) -length 4 -watermark 64 -ldpc $(LDPCDIR) -subrate 0 -delrate $(LOGSCALE) -duprate 0 >$@

len4.wat64.16.ldpc.dels.tab:
	$(ERRDECODE) -length 4 -watermark 64.16 -ldpc $(LDPCDIR) -subrate 0 -delrate $(LOGSCALE) -duprate 0 >$@

len4.wat64.16b.ldpc.dels.tab:
	$(ERRDECODE) -length 4 -watermark 64.16b -ldpc $(LDPCDIR) -subrate 0 -delrate $(LOGSCALE) -duprate 0 >$@

len4.wat64n.ldpc.dels.tab:
	$(ERRDECODE) -length 4 -watermark 64n -ldpc $(LDPCDIR) -subrate 0 -delrate $(LOGSCALE) -duprate 0 >$@

len4.wat64a2.ldpc.dels.tab:
	$(ERRDECODE) -length 4 -watermark 64a2 -ldpc $(LDPCDIR) -subrate 0 -delrate $(LOGSCALE) -duprate 0 >$@

len4.wat64.1.ldpc.dels.tab:
	$(ERRDECODE) -length 4 -watermark 64.1 -ldpc $(LDPCDIR) -subrate 0 -delrate $(LOGSCALE) -duprate 0 >$@

len4.wat64.1.ldpc.subs.tab:
	$(ERRDECODE) -length 4 -watermark 64.1 -ldpc $(LDPCDIR) -subrate $(LOGSCALE) -delrate 0 -duprate 0 >$@

# R
RSCRIPT = /Library/Frameworks/R.framework/Versions/3.2/Resources/Rscript

%.ps: %.r
	$(RSCRIPT) $<

# Paper rules:
%.ps: %.tex figures bib $(TEXINCLUDES)
	test -e $*.aux && rm $*.aux || eval
	latex -shell-escape $*
	bibtex $*
	latex -shell-escape $*
	latex -shell-escape $*
	dvips $*

# arxiv
ARXIV = trans.tex trans.bbl defs.tex notation.tex hamming31.ps hamming74.ps b2t.ps eof2.ps dna2full.ps dna2start.ps dna2startend.ps dna2norep.ps dna2startnorep.ps dna2startendnorep.ps blocks.ps error.ps dupplot.ps subplot.ps delplot.ps watplot.ps mixradtab.tex dnastoretab.tex
arxiv: $(ARXIV)
	test -e $@ || mkdir $@
	cp $(ARXIV) $@

arxiv.tar: $(ARXIV)
	tar cvf $@ $(ARXIV)

# hey make! don't delete all my stuff
.SECONDARY:
