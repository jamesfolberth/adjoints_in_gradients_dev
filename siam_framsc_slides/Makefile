
# Files
TEXSRC = $(shell egrep -l '^[^%]*\\begin\{document\}' *.tex)
PDF = $(TEXSRC:%.tex=%.pdf)
BIBSRC = $(shell *.bib)

# Commands
PDFLATEX = lualatex -halt-on-error -shell-escape -enable-write18
BIBTEX = bibtex

.PHONY: all pdflatex code bibtex clean cleanall figures force
.SUFFIXES: .tex .pdf .aux .log .out .dvi .ps .bib

all: code $(PDF)

$(PDF): $(TEXSRC) force
	$(PDFLATEX) $(TEXSRC)

#num_est_lyap_exp_pres.pdf: force num_est_lyap_exp_pres.tex
#	$(PDFLATEX) num_est_lyap_exp_pres.tex

#code:
#	make -C code

figures:
	make -C figures

#bibtex:
#	$(BIBTEX) $(BIBSRC)

force:

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
	$(TEXSRC:.tex=.ps) \
	$(TEXSRC:.tex=.auxlock) \
	$(TEXSRC:.tex=.nav) \
	$(TEXSRC:.tex=.snm)
	@rm -f tikz_cache/*


cleanall: clean
	#make -C code clean
	make -C figures clean
