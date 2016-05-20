# All source files
# Figures:
ALLDOT := $(wildcard *.dot)
ALLBIGDOT := $(wildcard *.bigdot)
ALLFDPDOT := $(wildcard *.fdpdot)
ALLFIG := $(subst .dot,,$(ALLDOT)) $(subst .bigdot,,$(ALLBIGDOT)) $(subst .fdpdot,,$(ALLFDPDOT))

ALLFIGPNG := $(addsuffix .png,$(ALLFIG))
ALLFIGPDF := $(addsuffix -fig.pdf,$(ALLFIG))

ALLFIGJUNK := $(addsuffix .tex,$(ALLFIG)) $(addsuffix .aux,$(ALLFIG)) $(addsuffix .log,$(ALLFIG)) $(addsuffix -pics.pdf,$(ALLFIG))

# Source files that only render properly as PNGs
BAD := 
BADPNG := $(addsuffix .png,$(BAD))
GOODPDF := $(addsuffix -fig.pdf,$(filter-out $(BAD),$(ALLFIG)))

FIGURES := $(BADPNG) $(GOODPDF)

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
# PDFs (requires dot2tex + graphviz + pyparsing)
# Landscape orientation, letter paper
%.tex: %.dot
	dot2tex --crop --format pstricks --texmode raw $< | perl -ne 'if(/documentclass/){print"\\documentclass{article}\n\\usepackage[x11names,rgb]{xcolor}\n\\usepackage{auto-pst-pdf}\n\\usepackage[landscape]{geometry}\n"}else{print unless /usepackage.*xcolor/}' >$@

# Landscape orientation, A1 paper (workaround to avoid cropping big graphs)
%.tex: %.bigdot
	dot2tex --crop --format pstricks --texmode raw $< | perl -ne 'if(/documentclass/){print"\\documentclass[a1paper]{article}\n\\usepackage[x11names,rgb]{xcolor}\n\\usepackage{auto-pst-pdf}\n\\usepackage[landscape]{geometry}\n"}else{print unless /usepackage.*xcolor/}' >$@

# Using TikZ instead of PSTricks (some problems with this, e.g. TikZ doesn't have same node shapes as normal GraphViz...)
%.tikz.tex: %.dot
	cat $*.dot | perl -pe 's/invhouse/pentagon/g;s/house/triangle/g' >$*.tmpdot
	dot2tex --format tikz --texmode raw $*.tmpdot | perl -ne 'if(/documentclass/){print"\\documentclass{article}\n\\usepackage[x11names,rgb]{xcolor}\n\\usepackage{auto-pst-pdf}\n\\usepackage[landscape]{geometry}\n"}else{print unless /usepackage.*xcolor/}' >$@

# LaTeX -> PDF
%-fig.pdf: %.tex
	pdflatex -shell-escape $*.tex
	pdfcrop $*.pdf $*-fig.pdf
	rm $*.pdf

# PNGs (requires graphviz only)
%.png: %.dot
	dot -Tpng $< >$@

%.png: %.bigdot
	dot -Tpng $< >$@

%.png: %.fdpdot
	dot -Tpng -Kfdp $< >$@

# hey make! don't delete all my stuff
.SECONDARY:
