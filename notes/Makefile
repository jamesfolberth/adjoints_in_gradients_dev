
# Files
TEXSRC = $(shell egrep -l '^[^%]*\\begin\{document\}' *.tex)
PDF = $(TEXSRC:%.tex=%.pdf)
BIBSRC = $(shell *.bib)

# Commands
PDFLATEX = lualatex -halt-on-error
BIBTEX = bibtex

.PHONY: all pdflatex bibtex clean cleanall figures force
.SUFFIXES: .tex .pdf .aux .log .out .dvi .ps .bib

all: figures pdflatex force

#########################################
# This will recompile LaTeX every time (regardless of LaTeX code update or not)
#pdflatex:  $(TEXSRC) 
#	$(PDFLATEX) $(TEXSRC)
#########################################
# This rule works for 1 tex file projects
pdflatex: $(PDF)

%.pdf: %.tex
	$(PDFLATEX) $(TEXSRC)
#########################################

figures:
	make -C figures

#bibtex:
#	$(BIBTEX) $(BIBSRC)


clean:
	@rm -f \
	$(TEXSRC:.tex=.aux) \
	$(TEXSRC:.tex=.log) \
	$(TEXSRC:.tex=.out) \
	$(TEXSRC:.tex=.bbl) \
	$(TEXSRC:.tex=.toc) \
	$(TEXSRC:.tex=.blg) \
	$(TEXSRC:.tex=.dvi) \
	$(TEXSRC:.tex=.pdf) \
	$(TEXSRC:.tex=.ps)


cleanall: clean
	make -C code clean
	make -C figures clean
